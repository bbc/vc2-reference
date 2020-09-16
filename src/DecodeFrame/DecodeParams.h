/*********************************************************************/
/* DecodeParams.h                                                    */
/* Author: Tim Borer and Galen Reich                                 */
/* This version September 2020                                       */
/*                                                                   */
/* Declares getting program parameters from command line.            */
/* Copyright (c) BBC 2011-2020 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef DECODEPARAMS_SEPT20
#define DECODEPARAMS_SEPT20

#include <string>
#include "Picture.h"
#include "WaveletTransform.h"

enum Output {TRANSFORM, QUANTISED, INDICES, DECODED};
enum Mode {HQ, LD};

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
  enum Output output;
  enum Mode mode;
  std::string error;

  int slice_scalar;
  int slice_prefix;
  int compressedBytes;
};

ProgramParams getCommandLineParams(int argc, char * argv[], const char* details[]);

#endif // DECODEPARAMS_SEPT20
