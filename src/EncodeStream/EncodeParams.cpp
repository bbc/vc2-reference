/*********************************************************************/
/* EncodeParams.cpp                                                  */
/* Author: Tim Borer and Galen Reich,  BBC Research and Development  */
/* This version July 2020                                            */
/*                                                                   */
/* Defines getting program parameters from command line.             */
/* Copyright (c) BBC 2011-2020 -- For license see the LICENSE file   */
/*********************************************************************/

#include "EncodeParams.h"
#include "Picture.h"
#include "WaveletTransform.h"

#include <iostream> //For cin, cout, cerr
#include <stdexcept> // For invalid_argument
#include <string>

using std::clog;
using std::cerr;
using std::endl;
using std::string;
using std::invalid_argument;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;

// Tell tclap that various enums are to be treated as tclap values 
namespace TCLAP {
  template <>
  struct ArgTraits<ColourFormat> { // Let TCLAP parse ColourFormat objects
    typedef ValueLike ValueCategory;
  };

  template <>
  struct ArgTraits<WaveletKernel> { // Let TCLAP parse WaveletKernel objects
    typedef ValueLike ValueCategory;
  };

  template <>
  struct ArgTraits<Output> { // Let TCLAP parse Output objects
    typedef ValueLike ValueCategory;
  };

  template <>
  struct ArgTraits<Mode> { // Let TCLAP parse Mode objects
    typedef ValueLike ValueCategory;
  };

}

ProgramParams getCommandLineParams(int argc, char * argv[], const char * details[]) {

  const char * const version = details[0];
  std::string description;
  {
    std::stringstream ss;
    ss << details[1] << "\n\n" << details[2] << "\n";
    description = ss.str();
  }

  if (argc<2) {
    clog << "Version: " << version << endl;
    clog << description << endl;
    clog << "\nFor more details and useage use -h or --help" << endl;
    exit(EXIT_SUCCESS);
  }

  ProgramParams params;

  try {

    // Define tclap command line object
    CmdLine cmd(description, ' ', version);

    // Define tclap command line parameters (and add them to tclap command line)
    UnlabeledValueArg<string> inFile("inFile", "Input file name (use \"-\" for standard input)", true, "-", "string", cmd);
    UnlabeledValueArg<string> outFile("outFile", "Output file name (use \"-\" for standard output)", true, "-", "string", cmd);
    SwitchArg verbosity("v", "verbose", "Output extra information to standard log", cmd);
    // "cla" prefix == command line argument
    ValueArg<Mode> cla_mode("m", "mode", "Encoding mode (HQ_ConstQ, HQ_CBR, LD [depreciated])", true, HQ_ConstQ, "string", cmd);
    ValueArg<Output> cla_output("o", "output", "Program output (Transform, Quantised, Indices, Packaged, Stream, Decoded, PSNR)", false, STREAM, "string", cmd);
    ValueArg<int> cla_hSliceSize("a", "hSlice", "Horizontical slice size (in units of 2**(wavelet depth))", true, 0, "integer", cmd);
    ValueArg<int> cla_vSliceSize("u", "vSlice", "Vertical slice size (in units of 2**(wavelet depth))", true, 0, "integer", cmd);
    ValueArg<int> cla_waveletDepth("d", "waveletDepth", "Wavelet transform depth", true, 0, "integer", cmd);
    ValueArg<WaveletKernel> cla_kernel("k", "kernel", "Wavelet kernel (DD97, LeGall, DD137, Haar0, Haar1, Fidelity, Daub97)", true, NullKernel, "string", cmd);
    SwitchArg cla_bottomFieldFirst("b", "bottomFieldFirst", "Bottom field is earliest (defaults to top field first))", cmd, false);
    SwitchArg cla_topFieldFirst("t", "topFieldFirst", "Top field is earliest (defaults to top field first))", cmd, true);
    SwitchArg cla_interlace("i", "interlace", "Use interlace coding (defaults to progressive coding))", cmd, false);
    SwitchArg cla_progressive("p", "progressive", "Use progressive coding (defaults to progressive coding))", cmd, true);
    ValueArg<int> cla_chromaDepth("c", "chromaDepth", "Bit depth for chroma (defaults to luma_depth), for RGB use -z)", false, 0, "integer", cmd);
    ValueArg<int> cla_lumaDepth("l", "lumaDepth", "Bit depth for luma (defaults to bits per input sample), for RGB use -z", false, 0, "integer", cmd);
    ValueArg<int> cla_bitDepth("z", "bitDepth", "Common bit depth for all components (defaults to bits per input sample)", false, 0, "integer", cmd);
    ValueArg<int> cla_bytes("n", "bytes", "Number of bytes per sample in image file (default 2)", false, 2, "integer", cmd);
    ValueArg<ColourFormat> cla_format("f", "format", "Colour format (4:4:4, 4:2:2, 4:2:0)", true, CF_UNSET, "string", cmd);
    ValueArg<int> cla_width("x", "width", "Picture width", true, 0, "integer", cmd);
    ValueArg<int> cla_height("y", "height", "Picture height", true, 0, "integer", cmd);
    ValueArg<int> cla_framerate("r", "framerate", "Frame Rate ( 1 = 24/1.001, 2 = 24, 3 = 25, 4 = 30/1.001, 5 = 30, 6 = 50, 7 = 60/1.001, 8 = 60, 9 = 15/1.001, 10 = 25/2, 11 = 48, 12=48/1.001, 13=96, 14=100, 15=120/1.001, 16=120 (default 3)", false, 3, "integer", cmd);
    
    ValueArg<int> cla_sliceScalar("S", "scalar", "Slice Size Scalar (default 1) [HQ_CBR and HQ_ConstQ modes only]", false, 1, "integer", cmd);
    ValueArg<int> cla_slicePrefix("P", "prefix", "Slice Prefix Bytes (default 0) [HQ_CBR and HQ_ConstQ modes only]", false, 0, "integer", cmd);
    ValueArg<int> cla_fragmentLength("F", "fragmentLength", "Maximum length in bytes for picture fragments (default = 0 = don't fragment) [HQ_CBR and LD modes only]", false, 0, "integer", cmd);
    ValueArg<int> cla_compressedBytes("s", "compressedBytes", "compressed bytes (size in bytes) [HQ_CBR and LD modes only]", false, 0, "integer", cmd);
    ValueArg<int> cla_quantIndex("q", "quantIndex", "Quantiser index (0 to 119) [HQ_ConstQ mode only]", false, 0, "integer", cmd);

    // Parse the argv array
    cmd.parse(argc, argv);

    // Initialise program parameters
    const string inFileName = inFile.getValue();
    const string outFileName = outFile.getValue();
    const bool verbose = verbosity.getValue();
    const int height = cla_height.getValue();
    const int width = cla_width.getValue();
    const ColourFormat chromaFormat = cla_format.getValue();
    const int bytes = cla_bytes.getValue();
    int bitDepth = cla_bitDepth.getValue();
    int lumaDepth = cla_lumaDepth.getValue();
    int chromaDepth = cla_chromaDepth.getValue();
    bool interlaced = cla_interlace.isSet();
    bool topFieldFirst = !cla_bottomFieldFirst.isSet();
    const WaveletKernel kernel = cla_kernel.getValue();
    const int waveletDepth = cla_waveletDepth.getValue();
    const int ySize = cla_vSliceSize.getValue();
    const int xSize = cla_hSliceSize.getValue();
    const Output output = cla_output.getValue();
    const Mode mode = cla_mode.getValue();
    const int frame_rate = cla_framerate.getValue();

    const int slice_scalar = cla_sliceScalar.getValue();
    const int slice_prefix = cla_slicePrefix.getValue();
    const int fragment_length = cla_fragmentLength.getValue();
    const int compressedBytes = cla_compressedBytes.getValue();
    const int qIndex = cla_quantIndex.getValue();

    // Check for valid combinations of parameters and options
    if (cla_bitDepth.isSet() && (cla_lumaDepth.isSet() || cla_chromaDepth.isSet()))
      throw invalid_argument("bitDepth is incompatible with luma depth (and/or chroma depth): use one or the other");
    if (cla_progressive.isSet() && cla_interlace.isSet())
      throw invalid_argument("image can't be both interlaced and progressive: specify one or the other");
    if (cla_progressive.isSet() && (cla_topFieldFirst.isSet() || cla_bottomFieldFirst.isSet()))
      throw invalid_argument("field parity is incompatible with progressive image");
    if (cla_topFieldFirst.isSet() && cla_bottomFieldFirst.isSet())
      throw invalid_argument("image can't be both top field first and bottom field first: specify one or the other");

    // Set default values
    if (!cla_bitDepth.isSet()) bitDepth = 8*bytes;
    if (!cla_lumaDepth.isSet()) lumaDepth = bitDepth;
    if (!cla_chromaDepth.isSet()) chromaDepth = lumaDepth;

    // Check parameter values
    if (height<1) throw invalid_argument("picture height must be > 0");
    if (width<1) throw invalid_argument("picture width must be > 0");
    if (chromaFormat==CF_UNSET)
      throw std::invalid_argument("unknown colour format");
    if ( (1>bytes) || (bytes>4) )
      throw std::invalid_argument("bytes must be in range 1 to 4");
    if (cla_bitDepth.isSet()) {
      if ( (1>bitDepth) || (bitDepth>(8*bytes)) )
        throw std::invalid_argument("bit depth must be in range 1 to 8*(bytes per sample)");
    }
    else {
      if ( (1>lumaDepth) || (lumaDepth>(8*bytes)) )
        throw std::invalid_argument("luma bit depth must be in range 1 to 8*(bytes per sample)");
      if ( (1>chromaDepth) || (chromaDepth>(8*bytes)) )
        throw std::invalid_argument("chroma bit depth must be in range 1 to 8*(bytes per sample)");
    }
    if (kernel==NullKernel)
      throw std::invalid_argument("invalid wavelet kernel");
    if (waveletDepth<1)
      throw std::invalid_argument("wavelet depth must be 1 or more");

    // Multiple modes option logic - commandline args are given correctly
    if (!(mode == HQ_CBR || mode == HQ_ConstQ) && cla_sliceScalar.isSet())
        throw std::invalid_argument("Slice Scalar is only used in HQ_CBR and HQ_ConstQ modes");
    if (!(mode == HQ_CBR || mode == HQ_ConstQ) && cla_slicePrefix.isSet())
        throw std::invalid_argument("Slice Prefix is only used in HQ_CBR and HQ_ConstQ modes");
    if (!(mode == HQ_CBR || mode == LD) && cla_fragmentLength.isSet())
        throw std::invalid_argument("Fragment length is only used in HQ_CBR and LD modes");
    if (!(mode == HQ_CBR || mode == LD) && cla_compressedBytes.isSet())
        throw std::invalid_argument("Compressed bytes is only used in HQ_CBR and LD modes");
    if (!(mode == HQ_ConstQ) && cla_quantIndex.isSet())
        throw std::invalid_argument("Quantisation index is only used in HQ_ConstQ mode");

    if ((mode == HQ_CBR || mode == LD) && !cla_compressedBytes.isSet())
        throw std::invalid_argument("Compressed bytes must be set in HQ_CBR and LD modes");
    if ((mode == HQ_ConstQ) && !cla_quantIndex.isSet())
        throw std::invalid_argument("Quantisation index must be set in HQ_ConstQ mode");


    // Multiple modes option logic - args are valid
    if ((mode == HQ_CBR || mode == HQ_ConstQ) && slice_scalar < 1) {
      throw std::invalid_argument("slice scalar must be >=1");
    }
    if ((mode == HQ_CBR || mode == HQ_ConstQ) && slice_prefix < 0) {
      throw std::invalid_argument("slice prefix must be >=0");
    }
    if ((mode == HQ_CBR || mode == LD) && compressedBytes<1)
      throw std::invalid_argument("number of compressed bytes must be >0");
    if (mode == HQ_ConstQ && ((qIndex<0) || (qIndex>119)))
      throw std::invalid_argument("quantisation index must be in the range 0 to 119");

    params.inFileName = inFileName;
    params.outFileName = outFileName;
    params.verbose = verbose;
    params.height = height;
    params.width = width;
    params.chromaFormat = chromaFormat;
    params.bytes = bytes;
    params.lumaDepth = lumaDepth;
    params.chromaDepth = chromaDepth;
    params.interlaced = interlaced;
    params.topFieldFirst = topFieldFirst;
    params.kernel = kernel;
    params.waveletDepth = waveletDepth;
    params.ySize = ySize;
    params.xSize = xSize;
    params.output = output;
    params.mode = mode;

    params.slice_scalar = slice_scalar;
    params.slice_prefix = slice_prefix;
    params.fragment_length = fragment_length;
    params.compressedBytes = compressedBytes;
    params.qIndex = qIndex;


    try {
      params.frame_rate = (FrameRate)frame_rate;
    } catch(const std::exception& e) {
      params.error = string("Invalid Frame Rate: ") + e.what();
    }
  }

  // catch any TCLAP exceptions
  catch (TCLAP::ArgException &e) {
    params.error = string("Command line error: ") + e.error() + " for arg " + e.argId();
  }

  // catch other exceptions
  catch(const std::exception& ex) {
    params.error = string(ex.what());
  }

  return params;
}

std::ostream& operator<<(std::ostream& os, Output output) {
  const char* s;
  switch (output) {
    case TRANSFORM:
      s = "Transform";
      break;
    case QUANTISED:
      s = "Quantised";
      break;
    case INDICES:
      s = "Indices";
      break;
    case PACKAGED:
      s = "Packaged";
      break;
    case STREAM:
      s = "Stream";
      break;
    case DECODED:
      s = "Decoded";
      break;
    case PSNR:
      s = "PSNR";
      break;
    default:
      s = "Unknown output!";
      break;
  }
  return os<<s;
}

std::istream& operator>>(std::istream& is, Output& output) {
        std::string text;
        is >> text;
        if (text == "Transform") output = TRANSFORM;
        else if (text == "Quantised") output = QUANTISED;
        else if (text == "Indices") output = INDICES;
        else if (text == "Packaged") output = PACKAGED;
        else if (text == "Stream") output = STREAM;
        else if (text == "Decoded") output = DECODED;
        else if (text == "PSNR") output = PSNR;
        else is.setstate(std::ios_base::badbit|std::ios_base::failbit);
        // Alternatively
        // else throw std::invalid_argument("invalid input");
        return is;
}

std::ostream& operator<<(std::ostream& os, Mode mode) {
  const char* s;
  switch (mode) {
    case HQ_CBR:
      s = "HQ_CBR";
      break;
    case HQ_ConstQ:
      s = "HQ_ConstQ";
      break;
    case LD:
      s = "LD";
      break;
    default:
      s = "Unknown mode!";
      break;
  }
  return os<<s;
}

std::istream& operator>>(std::istream& is, Mode& mode) {
        std::string text;
        is >> text;
        if (text == "HQ_ConstQ") mode = HQ_ConstQ;
        else if (text == "HQ_CBR") mode = HQ_CBR;
        else if (text == "LD") mode = LD;
        else is.setstate(std::ios_base::badbit|std::ios_base::failbit);
        // Alternatively
        // else throw std::invalid_argument("invalid input");
        return is;
}
