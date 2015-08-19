/*********************************************************************/
/* EncodeParams.h                                                    */
/* Author: Tim Borer,  BBC Research                                  */
/* This version 17th September 2013                                  */
/*                                                                   */
/* Declares getting program parameters from command line.            */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef ENCODERPARAMS_17SEPT13
#define ENCODERPARAMS_17SEPT13

#include <string>

#include "Picture.h"
#include "WaveletTransform.h"
#include "DataUnit.h"

enum Output {TRANSFORM, QUANTISED, INDICES, PACKAGED, STREAM, DECODED, PSNR};

std::ostream& operator<<(std::ostream&, Output value);

std::istream& operator>>(std::istream&, Output& value);

struct ProgramParams {
  std::string inFileName;
  std::string outFileName;
  bool verbose;
  int height;
  int width;
  enum ColourFormat chromaFormat;
  int bytes;
  int lumaDepth;
  int chromaDepth;
  bool interlaced;
  bool topFieldFirst;
  enum WaveletKernel kernel;
  int waveletDepth;
  int ySize;
  int xSize;
  int compressedBytes;
  enum Output output;
  FrameRate frame_rate;
  int slice_scalar;
  std::string error;
};

ProgramParams getCommandLineParams(int argc, char * argv[], const char* details[]);

#endif // ENCODERPARAMS_17SEPT13
