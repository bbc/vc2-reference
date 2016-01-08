/***********************************************************************/
/* DecodeHQ.cpp                                                        */
/* Author: Tim Borer                                                   */
/* This version 4th November 2013                                      */
/*                                                                     */
/* Reads compressed transform data in                                  */
/* Decompresses image using VC-2 High Quality profile                  */
/* Writes image data out to a planar file.                             */
/* It is not necessarily complet nor korrect.                          */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file     */
/***********************************************************************/

const char version[] = __DATE__ " @ " __TIME__ ;
const char summary[] = "Decodes the compressed bytes of a VC-2 High Quality profile to an uncompressed planar file";
const char description[] = "\
This program decodes SMPTE VC-2 HQ profile compressed transform data to regenerate an image sequence.\n\
Its primary output is the decoded image sequence. However it may produce alternative outputs which are:\n\
  1 the wavelet transform of the decoded output (inverse quantised wavelet coefficients)\n\
  2 the quantised wavelet coefficients\n\
  3 the quantisation indices used for each slice\n\
  4 the decoded sequence\n\
Input is just a sequence of compressed bytes.\n\
Output (where appropriate) are in planar format (4:4:4, 4:2:2, 4:2:0 or RGB).\n\
There can be 1 to 4 bytes per sample and the data is left (MSB) justified.\n\
Data is assumed offset binary (which is fine for both YCbCr or RGB).\n\
\n\
Example: DecodeHQ -v -x 1920 -y 1080 -f 4:2:2 -i -l 10 -k LeGall -d 3 -u 1 -a 2 inFileName outFileName";
const char* details[] = {version, summary, description};

#include <cstdlib> //for EXIT_SUCCESS, EXIT_FAILURE, atoi
#include <stdexcept> //For standard logic errors
#include <iostream> //For cin, cout, cerr
#include <string>
#include <fstream>
#include <cstdio> // for perror

#include "DecodeParams.h"
#include "Arrays.h"
#include "Slices.h"
#include "Picture.h"
#include "Frame.h"
#include "Quantisation.h"
#include "WaveletTransform.h"
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

int main(int argc, char * argv[]) {
  
try { //Giant try block around all code to get error messages

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
  const Output output = params.output;
  const int sliceScalar = params.slice_scalar;
  const int slicePrefix = params.slice_prefix;

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

  if (verbose) {
    clog << "bytes per sample= " << bytes << endl;
    clog << "luma depth (bits) = " << lumaDepth << endl;
    clog << "chroma depth (bits) = " << chromaDepth << endl;
    clog << "height = " << height << endl;
    clog << "width = " << width << endl;
    clog << "chroma format = " << chromaFormat << endl;
    clog << "interlaced = " << std::boolalpha << interlaced << endl;
    if (interlaced) clog << "top field first = " << std::boolalpha << topFieldFirst << endl;
    clog << "wavelet kernel = " << kernel << endl;
    clog << "wavelet depth = " << waveletDepth << endl;
    clog << "vertical slice size (in units of 2**(wavelet depth)) = " << ySize << endl;
    clog << "horizontal slice size (in units of 2**(wavelet depth)) = " << xSize << endl;
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

  const int framePics = (interlaced ? 2 : 1);

  // Construct an container to read the compressed data into.
  const PictureFormat transformFormat(paddedPictureHeight, paddedWidth, chromaFormat);
  Slices inSlices(transformFormat, waveletDepth, ySlices, xSlices);

  // Define picture format (field or frame)
  const PictureFormat picFormat(pictureHeight, width, chromaFormat);

  // Create Frame to hold output data
  const PictureFormat frameFormat(height, width, chromaFormat);
  Frame outFrame(frameFormat, interlaced, topFieldFirst);
  
  int frame = 0;
  while (true) {

    for (int pic=0; pic<framePics; ++pic) {

      // Read input from planar file
      if (verbose) {
        if (interlaced)
          clog << "Reading compressed input field " << pic << " of frame " << frame;
        else
          clog << "Reading compressed input frame number " << frame;
      }
      clog.flush(); // Make sure comments written to log file.
      inStream >> sliceio::highQualityVBR(slicePrefix, sliceScalar); // Read input in HQ VBR mode
      inStream >> inSlices; // Read the compressed input picture
      // Check picture was read OK
      if (!inStream) {
        if (frame==0) {
          cerr << "\rFailed to read the first compressed frame" << endl;
	        return EXIT_FAILURE;
        }
        else {
          if (verbose) clog << "\rEnd of input reached after " << frame << " frames     " << endl;
          if (inFileName!="-") inFileBuffer.close();
          if (outFileName!="-") outFileBuffer.close();
          return EXIT_SUCCESS;
        }
      }
      else clog << endl;
    
      // Reorder quantised coefficients from slice order to transform order
      if (verbose) clog << "Merge slices into full picture" << endl;
      const Picture yuvQCoeffs = merge_blocks(inSlices.yuvSlices);

      if (output==INDICES) {
        //Write quantisation indices as 1 byte unsigned values
        clog << "Writing quantisation indices to output file" << endl;
        outStream << arrayio::wordWidth(1); //1 byte per sample
        outStream << arrayio::unsigned_binary; // unsigned output
        outStream << inSlices.qIndices;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }
      
      if (output==QUANTISED) {
        //Write quantised transform output as 4 byte 2's comp values
        clog << "Writing quantised transform coefficients to output file" << endl;
        outStream << pictureio::wordWidth(4); // 4 bytes per sample
        outStream << pictureio::signed_binary; // 2's comp output
        outStream << yuvQCoeffs;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }
    
      // Inverse quantise in transform order
      if (verbose) clog << "Inverse quantise" << endl;
      const Picture yuvTransform = inverse_quantise_transform_np(yuvQCoeffs, inSlices.qIndices, qMatrix);

      if (output==TRANSFORM) {
        //Write transform output as 4 byte 2's comp values
        clog << "Writing transform coefficients to output file" << endl;
        outStream << pictureio::wordWidth(4); //4 bytes per sample
        outStream << pictureio::signed_binary;
        outStream << yuvTransform;
        if (!outStream) {
          cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
	        return EXIT_FAILURE; }
        continue; // omit rest of processing for this picture
      }

      // Inverse wavelet transform
      if (verbose) clog << "Inverse transform" << endl;
      const Picture outPicture = inverseWaveletTransform(yuvTransform, kernel, waveletDepth, picFormat);
  
      // Copy picture to output frame
      if (verbose) clog << "Copy picture to output frame" << endl;
      if (interlaced) {
        (pic==0) ? outFrame.firstField(outPicture) : outFrame.secondField(outPicture);
      }
      else { //progressive
        outFrame.frame(outPicture);
      }

    }

    if (verbose) clog << "Clipping output" << endl;
    {
      const int yMin = -utils::pow(2, lumaDepth-1);
      const int yMax = utils::pow(2, lumaDepth-1)-1;
      const int uvMin = -utils::pow(2, chromaDepth-1);
      const int uvMax = utils::pow(2, chromaDepth-1)-1;
      outFrame.frame(clip(outFrame, yMin, yMax, uvMin, uvMax));
    }

    if (verbose) clog << "Writing decoded output file" << endl;
    outStream << pictureio::wordWidth(bytes); // Set number of bytes per value in file
    outStream << pictureio::left_justified;
    outStream << pictureio::offset_binary;
    outStream << pictureio::bitDepth(lumaDepth, chromaDepth); // Set luma and chroma bit depths
    outStream << outFrame;
    if (!outStream) {
      cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
      return EXIT_FAILURE;
    }

    ++frame;

  } //End frame loop

} // end of try block

// Report error messages from try block
catch (const std::exception& ex) {
    cout << "Error: " << ex.what() << endl;
    return EXIT_FAILURE;
}

  // Program should never reach here!
  return EXIT_FAILURE;
}
