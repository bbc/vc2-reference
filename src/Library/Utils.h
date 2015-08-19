/*********************************************************************/
/* Utils.h                                                           */
/* Author: Tim Borer                                                 */
/* This version 7th July 2011                                        */
/*                                                                   */
/* Declares utility functions listed below (in namespace utils)      */
/*   fileSize: returns the size of a file in bytes                   */
/*   pow: raises an integer to an integer power                      */
/*   intlog2: the number of bits needed to express a number          */
/*   rationalise: the simplest form of a rational number             */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef UTILS_25FEB10
#define UTILS_25FEB10

#include <fstream>

namespace utils {

// Does what it says on the tin. Returns file size in bytes.
// Leaves the read pointer in the same position
// TO DO: change return type to size_t?
const int fileSize(const std::ifstream& infile);

// A little utility function, raises base to power exp
const int pow(int base, int exp);

// Another utility function which counts the bits needed to express a number
const int intlog2(int value);

struct Rational {
  int numerator;
  int denominator;
};

// Returns a rational number in which both numerator & denominator have been
// divided by their greatest common divisor.  
const Rational rationalise(const int numerator,
                           const int denominator);

// Two functions to set the mode of stdio
// Utility for setting the mode of stdin/stdout and cin/cout to either
// binary or text mode.
// The function actually changes the mode of stdin/out but since these
// use the same file id as cin/cout it changes the mode of those as well.
// This function is only really relevant to Windows OS. *nixes use binary
// IO mode all the time (there is no distinction beween binary and text mode).
// The function does nothing under *nixes.
// An argument is needed to control the mode. This is should be a
// platform independent type. I have used std::ios_base::openmode for this
// purpose. When a value of std::ios_base::binary is passed as a parameter
// then the stdio and cin/out streams are set to binary mode (on Windows OS).
// Return value: as _setmode function for Windows (-1 indicates error)
//              0 for *nix (always succeeds)

int setstdinmode(std::ios_base::openmode);
int setstdoutmode(std::ios_base::openmode);

}  // End namespace utils

// Output formating for rational on streams
std::ostream &operator << (std::ostream &stream, utils::Rational &r);

#endif //UTILS_25FEB10
