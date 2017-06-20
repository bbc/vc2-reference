/*********************************************************************/
/* EncodeParams.cpp                                                  */
/* Author: Tim Borer,  BBC Research                                  */
/* This version 19th September 2013                                  */
/*                                                                   */
/* Defines getting program parameters from command line.             */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
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
}

ProgramParams getCommandLineParams(int argc, char * argv[], const char * details[]) {

  const char * const version = details[0];
  std::string description;
  {
    std::stringstream ss;
    ss <<  details[1] << "\n\n" << details[2] << "\n";
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
    ValueArg<Output> cla_output("o", "output", "Program output (Transform, Quantised, Indices, Packaged, Stream, Decoded, PSNR)", false, STREAM, "string", cmd);
    ValueArg<int> cla_quantIndex("q", "quantIndex", "Quantiser index (0 to 63)", true, 0, "integer", cmd);
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
    ValueArg<ColourFormat> cla_format("f", "format", "Colour format (4:4:4, 4:2:2, 4:2:0 or RGB)", true, UNKNOWN, "string", cmd);
    ValueArg<int> cla_width("x", "width", "Picture width", true, 0, "integer", cmd);
    ValueArg<int> cla_height("y", "height", "Picture height", true, 0, "integer", cmd);
    ValueArg<int> cla_framerate("r", "framerate", "Frame Rate ( 1 = 24/1.001, 2 = 24, 3 = 25, 4 = 30/1.001, 5 = 30, 6 = 50, 7 = 60/1.001, 8 = 60, 9 = 15/1.001, 10 = 25/2, 11 = 48, 12=48/1.001, 13=96, 14=100, 15=120/1.001, 16=120(default 3)", false, 3, "integer", cmd);
    ValueArg<int> cla_sliceScalar("S", "scalar", "Slice size Scalar (default 1)", false, 1, "integer", cmd);
    ValueArg<int> cla_slicePrefix("P", "prefix", "Slice Prefix Bytes (default 0)", false, 0, "integer", cmd);

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
    int bitDepth = cla_bitDepth.getValue();;
    int lumaDepth = cla_lumaDepth.getValue();;
    int chromaDepth = cla_chromaDepth.getValue();;
    bool interlaced = cla_interlace.isSet();
    bool topFieldFirst = !cla_bottomFieldFirst.isSet();
    const WaveletKernel kernel = cla_kernel.getValue();
    const int waveletDepth = cla_waveletDepth.getValue();
    const int ySize = cla_vSliceSize.getValue();
    const int xSize = cla_hSliceSize.getValue();
    const int qIndex = cla_quantIndex.getValue();
    const Output output = cla_output.getValue();
    const int frame_rate = cla_framerate.getValue();
    const int sliceScalar = cla_sliceScalar.getValue();
    const int slicePrefix = cla_slicePrefix.getValue();

    // Check for valid combinations of parameters and options
    if ((chromaFormat==RGB) && (cla_lumaDepth.isSet() || cla_chromaDepth.isSet()))
      throw invalid_argument("luma/chroma depth is not appropriate for RGB (use -z or --bitDepth)");
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
    if (chromaFormat==UNKNOWN)
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
    if ((qIndex<0) || (qIndex>63))
      throw std::invalid_argument("quantisation index must be in the range 0 to 119");

    if (sliceScalar < 1)
      throw std::invalid_argument("slice size scalar must be at least 1");

    if (slicePrefix < 0)
      throw std::invalid_argument("slice prefix bytes must be at least 0");

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
    params.qIndex = qIndex;
    params.output = output;
    params.slice_scalar = sliceScalar;
    params.slice_prefix = slicePrefix;

    switch (frame_rate) {
    case 1:
      params.frame_rate = FR24000_1001;
      break;
    case 2:
      params.frame_rate = FR24;
      break;
    case 3:
      params.frame_rate = FR25;
      break;
    case 4:
      params.frame_rate = FR30000_1001;
      break;
    case 5:
      params.frame_rate = FR30;
      break;
    case 6:
      params.frame_rate = FR50;
      break;
    case 7:
      params.frame_rate = FR60000_1001;
      break;
    case 8:
      params.frame_rate = FR60;
      break;
    case 9:
      params.frame_rate = FR15000_1001;
      break;
    case 10:
      params.frame_rate = FR25_2;
      break;
    case 11:
      params.frame_rate = FR48;
      break;
    case 12:
      params.frame_rate = FR48_1001;
      break;
    case 13:
      params.frame_rate = FR96;
      break;
    case 14:
      params.frame_rate = FR100;
      break;
    case 15:
      params.frame_rate = FR120_1001;
      break;
    case 16:
      params.frame_rate = FR120;
      break;
    default:
      params.frame_rate = FR0;
      break;
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
        else if (text == "Packaged") output = PACKAGED;
        else if (text == "Stream") output = STREAM;
        else if (text == "Decoded") output = DECODED;
        else if (text == "PSNR") output = PSNR;
        else is.setstate(std::ios_base::badbit|std::ios_base::failbit);
        // Alternatively
        // else throw std::invalid_argument("invalid input");
        return is;
}
