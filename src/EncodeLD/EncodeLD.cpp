/*********************************************************************/
/* EncodeLD.cpp                                                      */
/* Author: Tim Borer                                                 */
/* This version 27th September 2013                                  */
/*                                                                   */
/* Reads image data in from a planar file.                           */
/* Compresses image using SMPTE VC-2 Low Delay profile.              */
/* Write compressed transform data out (not complete stream).        */
/* It is not necessarily complet nor korrect.                        */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

const char version[] = __DATE__ " @ " __TIME__ ;
const char summary[] = "Encodes an uncompressed planar video file with VC-2 Low Delay profile at constant bit rate";
const char description[] = "\
This program compresses an image sequence using SMPTE VC-2 LD profile.\n\
SMPTE VC-2 LD profile is strictly constant bit rate coding.\n\
The bit rate is specified by defining the number of compressed bytes per frame.\n\
Its primary output is the VC-2 stream. However it may produce alternative outputs which are:\n\
  1 the wavelet transform of the input\n\
  2 the quantised wavelet coefficients\n\
  3 the quantisation indices used for each slice\n\
  4 compressed bytes\n\
  5 VC2 bitstream (default output)\n\
  6 the decoded sequence\n\
  7 the PSNR for each frame\n\
Input and output (where appropriate) are in planar format (4:4:4, 4:2:2, 4:2:0 or RGB).\n\
There can be 1 to 4 bytes per sample and the data is left (MSB) justified.\n\
Data is assumed offset binary (which is fine for both YCbCr or RGB).\n\
\n\
Example: EncodeLD -v -x 1920 -y 1080 -f 4:2:2 -l 10 -k LeGall -d 3 -u 1 -a 2 -s 829440 -i inFileName outFileName";
const char* details[] = {version, summary, description};

#include <cstdlib> //for EXIT_SUCCESS, EXIT_FAILURE, atoi
#include <stdexcept> //For standard logic errors
#include <iostream> //For cin, cout, cerr
#include <string>
#include <fstream>
#include <cstdio> // for perror
#include <iomanip> // For reporting stats only
#include <algorithm>
#include <functional>
#include <cmath>

#include "EncodeParams.h"
#include "Arrays.h"
#include "Picture.h"
#include "Frame.h"
#include "WaveletTransform.h"
#include "Quantisation.h"
#include "Slices.h"
#include "Utils.h"

using std::cout;
using std::cin;
using std::cerr;
using std::clog;
using std::endl;
using std::string;
using std::filebuf;
using std::streambuf;
using std::ios_base;
using std::istream;
using std::ostream;

// Declare algorithm to calculate an array of quantisation indices
const Array2D quantIndices(const Picture& coefficients,
                           const Array1D& qMatrix,
                           const Array2D& sliceBytes);

int main(int argc, char * argv[]) {
  
try { //Giant try block around all code to get error messages

  // get command line parameters
  ProgramParams params = getCommandLineParams(argc, argv, details);
  if (!params.error.empty()) {
    cerr << "Command line error: " << params.error << endl;
    return EXIT_FAILURE;
  }

  // Create convenient aliases for program parameters
  const string inFileName = params.inFileName;
  const string outFileName = params.outFileName;
  const bool verbose = params.verbose;
  const int height = params.height;
  const int width = params.width;
  const ColourFormat chromaFormat = params.chromaFormat;
  const int bytes = params.bytes;
  const int lumaDepth = params.lumaDepth;
  const int chromaDepth = params.chromaDepth;
  const bool interlaced = params.interlaced;
  const bool topFieldFirst = params.topFieldFirst;
  const WaveletKernel kernel = params.kernel;
  const int waveletDepth = params.waveletDepth;
  const int ySize = params.ySize;
  const int xSize = params.xSize;
  const int compressedBytes = params.compressedBytes;
  const Output output = params.output;
  const FrameRate frame_rate = params.frame_rate;
  const int fragmentLength = params.fragment_length;

  if (verbose) {
    clog << endl;
    for (int arg=0; arg<argc; ++arg) { //Output command line
      if (arg) clog << " ";
      clog << argv[arg];
    }
    clog << endl;
    clog << "input file = " << inFileName << endl;
    clog << "output file = " << outFileName << endl;
  }

  // Open input file or use standard input
  // Input stream is read only binary mode.
  // No point in continuing if can't open input file.
  filebuf inFileBuffer; // For file input. Needs to be defined here to remain in scope
  streambuf *pInBuffer; // Either standard input buffer or a file buffer
  if (inFileName=="-") { // Use standard in
    //Set standard input to binary mode.
    //Only relevant for Windows (*nix is always binary)
    if ( utils::setstdinmode(std::ios_base::binary) == -1 ) {
        cerr << "Error: could not set standard input to binary mode" << endl;
        return EXIT_FAILURE;
    }
    pInBuffer = cin.rdbuf();
  }
  else { // Open file inFileName and use it for input
    pInBuffer = inFileBuffer.open(inFileName.c_str(), ios_base::in|ios_base::binary);
    if (!pInBuffer) {
      perror((string("Failed to open input file \"")+inFileName+"\"").c_str());
      return EXIT_FAILURE;
    }
  }
  istream inStream(pInBuffer);

  // Configure input stream to read the required picture format
  inStream >> pictureio::wordWidth(bytes); // Set number of bytes per value in file
  inStream >> pictureio::left_justified;
  inStream >> pictureio::offset_binary;
  inStream >> pictureio::bitDepth(lumaDepth, chromaDepth); // Set luma and chroma bit depths

  // Open output file or use standard output.
  // Output stream is write only binary mode
  // No point in continuing if can't open output file.
  filebuf outFileBuffer; // For file output. Needs to be defined here to remain in scope
  streambuf *pOutBuffer; // Either standard output buffer or a file buffer
  if (outFileName=="-") { // Use standard out
    //Set standard input to binary mode.
    //Only relevant for Windows (*nix is always binary)
    if ( utils::setstdoutmode(std::ios_base::binary) == -1 ) {
        cerr << "Error: could not set standard output to binary mode" << endl;
        return EXIT_FAILURE;
    }
    pOutBuffer = cout.rdbuf();
  }
  else { // Open file outFileName and use it for output
    pOutBuffer = outFileBuffer.open(outFileName.c_str(), ios_base::out|ios_base::binary);
    if (!pOutBuffer) {
      perror((string("Failed to open output file \"")+outFileName+"\"").c_str());
      return EXIT_FAILURE;
    }
  }
  ostream outStream(pOutBuffer);

  PictureFormat format(height, width, chromaFormat);

  if (verbose) {
    clog << "bytes per sample= " << bytes << endl;
    clog << "luma depth (bits) = " << lumaDepth << endl;
    clog << "chroma depth (bits) = " << chromaDepth << endl;
    clog << "height = " << format.lumaHeight() << endl;
    clog << "width = " << format.lumaWidth() << endl;
    clog << "chroma format = " << format.chromaFormat() << endl;
    clog << "interlaced = " << std::boolalpha << interlaced << endl;
    if (interlaced) clog << "top field first = " << std::boolalpha << topFieldFirst << endl;
    clog << "wavelet kernel = " << kernel << endl;
    clog << "wavelet depth = " << waveletDepth << endl;
    clog << "vertical slice size (in units of 2**(wavelet depth)) = " << ySize << endl;
    clog << "horizontal slice size (in units of 2**(wavelet depth)) = " << xSize << endl;
    clog << "compressed bytes = " << compressedBytes << endl;
    clog << "output = " << output << endl;
  }

  // Calculate number of slices per picture
  const int yTransformSize = ySize*utils::pow(2,waveletDepth);
  const int xTransformSize = xSize*utils::pow(2,waveletDepth);
  const int pictureHeight = ( (interlaced) ? height/2 : height);
  const int paddedPictureHeight = paddedSize(pictureHeight, waveletDepth);
  const int paddedWidth = paddedSize(width, waveletDepth);
  const int ySlices = paddedPictureHeight/yTransformSize;
  const int xSlices = paddedWidth/xTransformSize;
  if (paddedPictureHeight != (ySlices*yTransformSize) ) {
    throw std::logic_error("Padded picture height is not divisible by slice height");
	  return EXIT_FAILURE;
  } 
  if (paddedWidth != (xSlices*xTransformSize) ) {
    throw std::logic_error("Padded width is not divisible by slice width");
	  return EXIT_FAILURE;
  }

  if (verbose) {
    clog << "Vertical slices per picture          = " << ySlices << endl;
    clog << "Horizontal slices per picture        = " << xSlices << endl;
    // Calculate slice bytes numerator and denominator
    const utils::Rational sliceBytesNandD =
      utils::rationalise((interlaced ? compressedBytes/2 : compressedBytes), (ySlices*xSlices));
    const int SliceBytesNum = sliceBytesNandD.numerator;
    const int SliceBytesDenom = sliceBytesNandD.denominator;
    clog << "Slice bytes numerator                = " << SliceBytesNum << endl;
    clog << "Slice bytes denominator              = " << SliceBytesDenom << endl;
  }

  // Calculate the quantisation matrix
  const Array1D qMatrix = quantMatrix(kernel, waveletDepth);
  if (verbose) {
    clog << "Quantisation matrix = " << qMatrix[0];
    for (unsigned int i=1; i<qMatrix.size(); ++i) {
      clog << ", " << qMatrix[i];
    }
    clog << endl;
  }

  const int pictureSlices = ySlices*xSlices;
  const int framePics = (interlaced ? 2 : 1);
  const int totalSlices = framePics*pictureSlices;

  //Create input & output frames
  Frame inFrame(format, interlaced, topFieldFirst);
  Frame outFrame(format, interlaced, topFieldFirst); //used for decoded picture only
  Array2D yDiff(format.lumaShape()); // used for PSNR calculation only
  Array2D uDiff(format.chromaShape()); // used for PSNR calculation only
  Array2D vDiff(format.chromaShape()); // used for PSNR calculation only

  int frame = 0;
  if (output==STREAM) {
    if (verbose) clog << endl << "Writing Sequence Header" << endl << endl;
    outStream << dataunitio::start_sequence;
    if (fragmentLength > 0)
      outStream << dataunitio::fragmentedPictures(fragmentLength);
    outStream << SequenceHeader(PROFILE_LD, format.lumaHeight(), format.lumaWidth(), format.chromaFormat(), interlaced, frame_rate, topFieldFirst, lumaDepth);
  }
  while (true) {

    // Read input from planar file
    if (verbose) clog << "Reading input frame number " << frame;
    inStream >> inFrame; // Read the input frame
    // Check frame was read OK
    if (!inStream) {
      if (frame==0) {
        cerr << "\rFailed to read input frame number " << frame << endl;
	      return EXIT_FAILURE;
      }
      else {
        if (verbose) clog << "\rEnd of input reached after " << frame << " frames" << endl;
        break;
      }
    }
    else if (verbose) clog << endl;

    int stats[128] = {0}; //Define and initialise array to hold quantiser stats

    Picture picture; // define a picture, either a field or frame
    for (int pic=0 ; pic<framePics; ++pic) { // loop over fields if interlaced input

      if (interlaced) { // assign picture to be either a field or the whole frame
        picture = (pic==0 ? inFrame.firstField() : inFrame.secondField());
      }
      else { //progressive
        picture = inFrame;
      }

      //Forward wavelet transform
      if (verbose) clog << "Forward transform" << endl;
      Picture transform = waveletTransform(picture, kernel, waveletDepth);

      if (output==TRANSFORM) {
        //Write transform output as 4 byte 2's comp values
        clog << "Writing transform coefficients to output file" << endl;
        outStream << pictureio::wordWidth(4); //4 bytes per sample
        outStream << pictureio::signed_binary;
        outStream << transform;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }
      
      // Choose quantisation indices to achieve a size of compressedBytes for the frame
      if (verbose) clog << "Determine quantisation indices" << endl;
      const int pictureBytes = (interlaced ? compressedBytes/2 : compressedBytes);
      // Calculate number of bytes for each slice
      const Array2D bytes = slice_bytes(ySlices, xSlices, pictureBytes, 1);
      Array2D qIndices = quantIndices(transform, qMatrix, bytes);
    
      // Analyse quantiser index stats
      for (int v=0; v<ySlices; ++v) {
        for (int h=0; h<xSlices; ++h) {
          ++stats[qIndices[v][h]];
        }
      }

      if (output==INDICES) {
        //Write quantisation indices as 1 byte unsigned values
        clog << "Writing quantisation indices to output file" << endl;
        outStream << arrayio::wordWidth(1); //1 byte per sample
        outStream << arrayio::unsigned_binary; // unsigned output
        outStream << qIndices;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }

      if (verbose) clog << "Quantise transform coefficients" << endl;
      const Picture quantisedSlices = quantise_transform(transform, qIndices, qMatrix);
      
      if (output==QUANTISED) {
        //Write quantised transform output as 4 byte 2's comp values
        clog << "Writing quantised transform coefficients to output file" << endl;
        outStream << pictureio::wordWidth(4); // 4 bytes per sample
        outStream << pictureio::signed_binary; // 2's comp output
        outStream << quantisedSlices;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }

      // Split transform into slices
      if (verbose) clog << "Split quantised coefficients into slices" << endl;
      const PictureArray slices = split_into_blocks(quantisedSlices, ySlices, xSlices);
      
      if (output==PACKAGED) {
        // Package up data for output
        const Slices outSlices(slices, waveletDepth, qIndices);

        //Write packaged output
        if (verbose) clog << "Writing compressed output to file" << endl;
        outStream << sliceio::lowDelay(bytes); // Write output in Low Delay mode
        outStream << outSlices;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE;
        }
        continue; // omit rest of processing for this picture
      }

      if (output==STREAM) { 
        const utils::Rational rationalBytes = utils::rationalise(pictureBytes, (ySlices*xSlices));
        // Package up data for output
        const Slices outSlices(slices, waveletDepth, qIndices);
        const WrappedPicture outWrapped(frame,
                                        kernel,
                                        waveletDepth,
                                        xSlices,
                                        ySlices,
                                        rationalBytes,
                                        outSlices);

        //Write packaged output
        if (verbose) clog << "Writing compressed picture to file" << endl;
        outStream << sliceio::lowDelay(bytes); // Write output in Low Delay mode
        outStream << outWrapped;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE;
        }
        continue; // omit rest of processing for this picture
      }
    
      // Inverse quantise in transform order
      if (verbose) clog << "Inverse quantise" << endl;
      const Picture yuvTransform = inverse_quantise_transform(quantisedSlices, qIndices, qMatrix);

      // Inverse wavelet transform
      if (verbose) clog << "Inverse transform" << endl;
      picture = inverseWaveletTransform(yuvTransform, kernel, waveletDepth, picture.format());

      if (verbose) clog << "Clip decoded picture" << endl;
      {
        const int yMin = -utils::pow(2, lumaDepth-1);
        const int yMax = utils::pow(2, lumaDepth-1)-1;
        const int uvMin = -utils::pow(2, chromaDepth-1);
        const int uvMax = utils::pow(2, chromaDepth-1)-1;
        picture = (clip(picture, yMin, yMax, uvMin, uvMax));
      }

      // Assign either a field or the whole frame to outFrame
      if (interlaced) {
        (pic==0) ? outFrame.firstField(picture) : outFrame.secondField(picture);
      }
      else { //progressive
        outFrame.frame(picture);
      }

    } // end picture loop

    // Calculate mean and standard deviation of quantiser index
    float mean, meanSquare, stdDev;
    if (output!=TRANSFORM) { // No quant stats if output is TRANSFORM
      int topIndex = -1;
      for (int i=0; i<128; ++i) if (stats[i]>0) topIndex=i;
      mean = 0.0;
      meanSquare = 0.0;
      for (int z=0; z<=topIndex; ++z) {
        mean += (z*stats[z]);
        meanSquare += (z*z*stats[z]);
      }
      mean /= totalSlices;
      meanSquare /= totalSlices;
      stdDev = sqrt(meanSquare - (mean*mean));
      if (verbose) {
        clog << endl;
        clog << std::fixed << std::setprecision(2);
        clog << "Mean, Standard Deviation of quantiser index = " << mean << ", " << stdDev << endl;
      }
    }

    // Calculate PSNR of decoded frame
    float YPSNR, UPSNR, VPSNR;
    if ((output==DECODED) || (output==PSNR)) {
      // Calculate difference and difference square picture, and sum of squares
      // Y or R component
      std::transform(inFrame.y().data(), inFrame.y().data()+inFrame.y().num_elements(),
                     outFrame.y().data(),
                     yDiff.data(),
                     std::minus<int>() );
      std::transform(yDiff.data(), yDiff.data()+yDiff.num_elements(),
                     yDiff.data(),
                     yDiff.data(),
                     std::multiplies<int>() );
      const long long YSS = std::accumulate(yDiff.data(), yDiff.data()+yDiff.num_elements(), 0LL);
      const int yPixels = width*height;
      const float YRMS = sqrt(float(YSS)/float(yPixels))/utils::pow(2, lumaDepth);
      YPSNR = -20*log10(YRMS);
      // U or G component
      const int uvPixels = format.chromaWidth()*format.chromaHeight();
      std::transform(inFrame.c1().data(), inFrame.c1().data()+inFrame.c1().num_elements(),
                     outFrame.c1().data(),
                     uDiff.data(),
                     std::minus<int>() );
      std::transform(uDiff.data(), uDiff.data()+uDiff.num_elements(),
                     uDiff.data(),
                     uDiff.data(),
                     std::multiplies<int>() );
      const long long USS = std::accumulate(uDiff.data(), uDiff.data()+uDiff.num_elements(), 0LL);
      const float URMS = sqrt(float(USS)/float(uvPixels))/utils::pow(2, chromaDepth);
      UPSNR = -20*log10(URMS);
      // V or B component
      std::transform(inFrame.c2().data(), inFrame.c2().data()+inFrame.c2().num_elements(),
                     outFrame.c2().data(),
                     vDiff.data(),
                     std::minus<int>() );
      std::transform(vDiff.data(), vDiff.data()+vDiff.num_elements(),
                     vDiff.data(),
                     vDiff.data(),
                     std::multiplies<int>() );
      const long long VSS = std::accumulate(vDiff.data(), vDiff.data()+vDiff.num_elements(), 0LL);
      const float VRMS = sqrt(float(VSS)/float(uvPixels))/utils::pow(2, chromaDepth);
      VPSNR = -20*log10(VRMS);
      if (verbose) {
        clog << std::fixed << std::setprecision(4);
        clog << "PSNR for Y/R, U/G, V/B = " << YPSNR << ", " << UPSNR << ", " << VPSNR  << endl;
      }
    }

    if (output==DECODED) {
      if (verbose) clog << "Writing decoded output frame " << frame << endl;
      outStream << pictureio::wordWidth(bytes); // Define output word width
      outStream << pictureio::offset_binary; // Write output as offset binary
      outStream << pictureio::bitDepth(lumaDepth, chromaDepth); // Set bit depths
      outStream << outFrame;
      if (!outStream) {
        cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
        return EXIT_FAILURE;
      }
    }

    if (output==PSNR) {
        outStream << "Frame " << frame << endl;
        outStream << std::fixed << std::setprecision(2);
        outStream << mean << " " << stdDev << endl;
        outStream << std::fixed << std::setprecision(4);
        outStream << YPSNR << " " << UPSNR << " " << VPSNR  << endl;
    }

    ++frame;
  } //End frame loop
  if (output==STREAM) {
    outStream << dataunitio::end_sequence;
  }
  
  if (inFileName!="-") inFileBuffer.close();
  if (outFileName!="-") outFileBuffer.close();

} // end of try block

// Report error messages from try block
catch (const std::exception& ex) {
    cout << "Error: " << ex.what() << endl;
    return EXIT_FAILURE;
}

  return EXIT_SUCCESS;
}

class SliceQuantiserRef: public SliceQuantiser {
  public:
    SliceQuantiserRef(const Array2D& coefficients,
                      int vSlices, int hSlices,
                      const Array1D& quantMatrix);
    virtual const Array2D& quantise_slice(int qIndex);
  private:
    const Array2D& coeffs; //Keep a reference to transform coefficients
    Array2D decodedLLCoeffs; //Locally decoded LL subband coefficients
    Array2D qMatrix; //quantMatrix values in slice format
};

SliceQuantiserRef::SliceQuantiserRef(const Array2D& coefficients,
                                     int vSlices, int hSlices,
                                     const Array1D& quantMatrix):
  SliceQuantiser(coefficients, vSlices, hSlices, quantMatrix),
  coeffs(coefficients)
{
  const int LLHeight = coeffsHeight/transformSize;
  const int LLWidth = coeffsWidth/transformSize;
  decodedLLCoeffs.resize(extents[LLHeight][LLWidth]);
  qMatrix.resize(extents[sliceHeight][sliceWidth]);
  // Fill qMatrix with the appropriate values
  BlockVector xQMatrix = split_into_subbands(qMatrix, waveletDepth);
  for (int band=0; band<numberOfSubbands; ++band) {
    std::fill(xQMatrix[band].data(),
              xQMatrix[band].data()+xQMatrix[band].num_elements(),
              quantMatrix[band]);
  }
  qMatrix = merge_subbands(xQMatrix);
}

const Array2D& SliceQuantiserRef::quantise_slice(int qIndex) {
  for (int y=0, yPos=v*sliceHeight; y<sliceHeight; ++y, ++yPos) {
    for (int x=0, xPos=h*sliceWidth; x<sliceWidth; ++x, ++xPos) {
      const int adjustedQ = adjust_quant_index(qIndex, qMatrix[y][x]);
      // Note: for efficiency test could be done by bit compare of LSBs
      if ( ((y%transformSize)==0)&&((x%transformSize)==0) ) { // LL Subband
        // Note: For efficiency division could be done by a bit shift.
        const int yLL = yPos/transformSize; // index to LL subband
        const int xLL = xPos/transformSize; // index to LL subband
        const int prediction = predictDC(decodedLLCoeffs, yLL, xLL);
        qSlice[y][x] = quant(coeffs[yPos][xPos]-prediction, adjustedQ);
        decodedLLCoeffs[yLL][xLL] = scale(qSlice[y][x], adjustedQ)+prediction;
      }
      else {
        qSlice[y][x] = quant(coeffs[yPos][xPos], adjustedQ);
      }
    }
  }
  return qSlice;
}

const Array2D quantIndices(const Picture& coefficients,
                           const Array1D& qMatrix,
                           const Array2D& sliceBytes) {
  const int ySlices = sliceBytes.shape()[0];
  const int xSlices = sliceBytes.shape()[1];
  // Create an empty array of indices to fill and return
  Array2D indices(extents[ySlices][xSlices]); 
  // Wavelet depth & number of subbands derived from dimensions of qMatrix
  const int numberOfSubbands = qMatrix.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  // Create state machines to quantise slices, with trial quantisers, in raster order
  SliceQuantiserRef ySliceQuantiser(coefficients.y(), ySlices, xSlices, qMatrix);
  SliceQuantiserRef uSliceQuantiser(coefficients.c1(), ySlices, xSlices, qMatrix);
  SliceQuantiserRef vSliceQuantiser(coefficients.c2(), ySlices, xSlices, qMatrix);
  bool sliceAvailable = true;
  while (sliceAvailable) {
    const int bytes = sliceBytes[ySliceQuantiser.row()][ySliceQuantiser.column()];
    const int length_bits = utils::intlog2(8*bytes-7);
    const int bitsAvailable = 8*bytes - 7 - length_bits;

    int trialQ = 63;
    int q = 127;
    int delta = 64;
    while (delta>0) {
      delta >>= 1;
      const Array2D yTrialSlice = ySliceQuantiser.quantise_slice(trialQ);
      const Array2D uTrialSlice = uSliceQuantiser.quantise_slice(trialQ);
      const Array2D vTrialSlice = vSliceQuantiser.quantise_slice(trialQ);
      int bitsRequired = luma_slice_bits(yTrialSlice, waveletDepth);
      bitsRequired += chroma_slice_bits(uTrialSlice, vTrialSlice, waveletDepth);
      if (bitsRequired<=bitsAvailable) {
        if (trialQ<q) q=trialQ;
        trialQ -= delta;
      }
      else {
        trialQ +=delta;
      }
    }
    // The slices must be requantised with the correct q to ensure correct DC prediction
    ySliceQuantiser.quantise_slice(q);
    uSliceQuantiser.quantise_slice(q);
    vSliceQuantiser.quantise_slice(q);
    indices[ySliceQuantiser.row()][ySliceQuantiser.column()] = q;

    // Go to next slice
    sliceAvailable = ySliceQuantiser.next_slice();
    uSliceQuantiser.next_slice();
    vSliceQuantiser.next_slice();
  }
  return indices;
}
