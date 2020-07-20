/*********************************************************************/
/* DecodeParams.h                                                    */
/* Author: Tim Borer                                                 */
/* This version 20th September 2013                                  */
/*                                                                   */
/* Declares getting program parameters from command line.            */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef DECODEPARAMS_20SEPT13
#define DECODEPARAMS_20SEPT13

#include <string>
#include "Picture.h"
#include "WaveletTransform.h"

enum Output {TRANSFORM, QUANTISED, INDICES, DECODED};

std::ostream& operator<<(std::ostream&, Output value);

std::istream& operator>>(std::istream&, Output& value);

struct ProgramParams {
  std::string inFileName;
  std::string outFileName;
  bool verbose;
  enum Output output;
  std::string error;
};

ProgramParams getCommandLineParams(int argc, char * argv[], const char* details[]);

#endif // DECODEPARAMS_20SEPT13
