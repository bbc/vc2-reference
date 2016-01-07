/*********************************************************************/
/* DecodeParams.cpp                                                  */
/* Author: Tim Borer                                                 */
/* This version 20th September 2013                                  */
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

    // Parse the argv array
    cmd.parse(argc, argv);

    // Initialise program parameters
    const string inFileName = inFile.getValue();
    const string outFileName = outFile.getValue();
    const bool verbose = verbosity.getValue();
    const Output output = cla_output.getValue();

    params.inFileName = inFileName;
    params.outFileName = outFileName;
    params.verbose = verbose;
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
