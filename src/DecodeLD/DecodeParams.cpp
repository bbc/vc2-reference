/*********************************************************************/
/* DecodeParams.cpp                                                  */
/* Author: Tim Borer                                                 */
/* This version 4th November 2013                                    */
/*                                                                   */
/* Defines getting program parameters from command line.             */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include "DecodeParams.h"
#include "Picture.h"
#include "WaveletTransform.h"

#include <iostream> //For cin, cout, cerr, clog
#include <stdexcept> // For invalid_argument
#include <string>

using std::clog;
using std::endl;
using std::string;
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

ProgramParams getCommandLineParams(int argc, char * argv[], const char* details[]) {

  const char * const version = details[0];
  const char * const description = details[2];

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
    UnlabeledValueArg<string> inFile("inFile", "Input file name", true, "-", "string", cmd);
    UnlabeledValueArg<string> outFile("outFile", "Output file name", true, "-", "string", cmd);
    SwitchArg verbosity("v", "verbose", "Output extra information to standard log", cmd);
    // "cla" prefix == command line argument
    ValueArg<Output> cla_output("o", "output", "Program output (Transform, Quantised, Indices, Decoded)", false, DECODED, "string", cmd);
    ValueArg<int> cla_compressedBytes("s", "compressedBytes", "compressed bytes (size in bytes)", true, 0, "integer", cmd);
    ValueArg<int> cla_hSliceSize("a", "hSlice", "Horizontical slice size (in units of 2**(wavelet depth))", true, 0, "integer", cmd);
    ValueArg<int> cla_vSliceSize("u", "vSlice", "Vertical slice size (in units of 2**(wavelet depth))", true, 0, "integer", cmd);
    ValueArg<int> cla_waveletDepth("d", "waveletDepth", "Wavelet transform depth", true, 0, "integer", cmd);
    ValueArg<WaveletKernel> cla_kernel("k", "kernel", "Wavelet kernel (DD97, LeGall, DD137, Haar0, Haar1, Fidelity, Daub97)", true, NullKernel, "string", cmd);
    SwitchArg cla_bottomFieldFirst("b", "bottomFieldFirst", "Bottom field is earliest (defaults to top field first))", cmd, false);
    SwitchArg cla_topFieldFirst("t", "topFieldFirst", "Top field is earliest (defaults to top field first))", cmd, true);
    SwitchArg cla_interlace("i", "interlace", "Coded using interlace coding (defaults to progressive coding))", cmd, false);
    SwitchArg cla_progressive("p", "progressive", "Coded using progressive coding (defaults to progressive coding))", cmd, true);
    ValueArg<int> cla_chromaDepth("c", "chromaDepth", "Bit depth for chroma (defaults to luma_depth), for RGB use -z)", false, 0, "integer", cmd);
    ValueArg<int> cla_lumaDepth("l", "lumaDepth", "Bit depth for luma (defaults to bits per input sample), for RGB use -z", false, 0, "integer", cmd);
    ValueArg<int> cla_bitDepth("z", "bitDepth", "Common bit depth for all components (defaults to bits per input sample)", false, 0, "integer", cmd);
    ValueArg<int> cla_bytes("n", "bytes", "Number of bytes per sample in image file (default 2)", false, 2, "integer", cmd);
    ValueArg<ColourFormat> cla_format("f", "format", "Colour format (4:4:4, 4:2:2, 4:2:0 or RGB)", true, UNKNOWN, "string", cmd);
    ValueArg<int> cla_width("x", "width", "Picture width", true, 0, "integer", cmd);
    ValueArg<int> cla_height("y", "height", "Picture height", true, 0, "integer", cmd);

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
    const int compressedBytes = cla_compressedBytes.getValue();
    const Output output = cla_output.getValue();

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
    interlaced = cla_interlace.isSet();
    topFieldFirst = !cla_bottomFieldFirst.isSet();

    // Check parameter values
    if (height<1) throw invalid_argument("picture height must be > 0");
    if (width<1) throw invalid_argument("picture width must be > 0");
    if (chromaFormat==UNKNOWN)
      throw std::invalid_argument("unknown colour format");
    if ( (1>bytes) | (bytes>4) )
      throw std::invalid_argument("bytes must be in range 1 to 4");
    if (cla_bitDepth.isSet()) {
      if ( (1>bitDepth) | (bitDepth>(8*bytes)) )
        throw std::invalid_argument("bit depth must be in range 1 to 8*(bytes per sample)");
    }
    else {
      if ( (1>lumaDepth) | (lumaDepth>(8*bytes)) )
        throw std::invalid_argument("luma bit depth must be in range 1 to 8*(bytes per sample)");
      if ( (1>chromaDepth) | (chromaDepth>(8*bytes)) )
        throw std::invalid_argument("chroma bit depth must be in range 1 to 8*(bytes per sample)");
    }
    if (kernel==NullKernel)
      throw std::invalid_argument("invalid wavelet kernel");
    if (waveletDepth<1)
      throw std::invalid_argument("wavelet depth must be 1 or more");
    if (compressedBytes<1)
      throw std::invalid_argument("number of compressed bytes must be >0");

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
    params.compressedBytes = compressedBytes;
    params.output = output;

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
    case DECODED:
      s = "Decoded";
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
        else if (text == "Decoded") output = DECODED;
        else is.setstate(std::ios_base::badbit|std::ios_base::failbit);
        // Alternatively
        // else throw std::invalid_argument("invalid input");
        return is;
}
