/*********************************************************************/
/* Utils.cpp                                                         */
/* Author: Tim Borer and Galen Reich                                 */
/* This version July 2020                                            */
/*                                                                   */
/* Defines utility functions listed below (in namespace utils)       */
/*   fileSize: returns the size of a file in bytes                   */
/*   pow: raises an integer to an integer power                      */
/*   intlog2: the number of bits needed to express a number          */
/*   rationalise: the simplest form of a rational number             */
/* Copyright (c) BBC 2011-2020 -- For license see the LICENSE file   */
/*********************************************************************/

#include "Utils.h"

#include "boost/math/common_factor_rt.hpp"  // For gcd function

// Does what it says on the tin. Returns file size in bytes.
// Leaves the read pointer in the same position
// TO DO: change return type to size_t?
const int utils::fileSize(const std::ifstream& infile) {
  // const cast is becuase the function is logically const but
  // needs to modify the fstream then restore it.
  std::ifstream& file = const_cast<std::ifstream&>(infile);
  std::ifstream::pos_type filePosn = file.tellg();
  file.seekg(0, std::ios::end);
  int size = static_cast<int>(std::streamoff(file.tellg()));
  file.seekg(filePosn);
  return size;
}

// A little utility function, raises base to power exp
const int utils::pow(int base, int exp) {
  int value = 1;
  if (exp>0) while (exp--) value *= base;
  return value;
}

// Another utility function which counts the bits needed to express a number
const int utils::intlog2(int value) {
  int log = 0;
  --value;
  while (value>0) {
    value >>= 1;
    ++log;
  }
  return log;
}

// Calculate the picture number given a frame number, picture number, and number 
// of pictures per frame - 1 or 2 (interlaced)
unsigned long utils::getPictureNumber(int fieldNumber, unsigned long long frameNumber, const int fieldsPerFrame){

  // Input checking
  if (fieldNumber < 0) throw std::logic_error("field number should be positive");
  if (fieldNumber > fieldsPerFrame) throw std::logic_error("field number exceeds number of fields per frame");
  if (fieldsPerFrame!=1 && fieldsPerFrame!=2) throw std::logic_error("number of fields per frame should be 1 (progressive) or 2 (interlaced)");

  const unsigned long long bigPictureNumber = fieldNumber + frameNumber*fieldsPerFrame;
  const unsigned long long bitLimit32 = (1ULL) << 32;
  // Wrap round to 0 after 2^32 -1 bits
  return bigPictureNumber % bitLimit32;
}

const utils::Rational utils::rationalise(const int numerator,
                                         const int denominator) {
  const int gcd = boost::math::gcd(numerator, denominator);
  utils::Rational result;
  result.numerator = numerator/gcd;
  result.denominator = denominator/gcd;
  return result;
}

// Utilities to set iomode for stdio
#ifdef _WIN32

#include <stdio.h>  //Defines _fileno (needed to set standard i/o to binary mode)
#include <io.h>     //Defines _setmode (needed to set standard input to binary mode)
#include <fcntl.h>  //Contains definition of _O_BINARY

int utils::setstdinmode(std::ios_base::openmode mode) {
	int winMode;
	if ((mode&std::ios_base::binary)==std::ios_base::binary) winMode=_O_BINARY;
	else winMode=_O_TEXT;
    //Set standard input and standard output to binary mode.
	return _setmode(_fileno( stdin ), winMode );
}

int utils::setstdoutmode(std::ios_base::openmode mode) {
	int winMode;
	if ((mode&std::ios_base::binary)==std::ios_base::binary) winMode=_O_BINARY;
	else winMode=_O_TEXT;
    //Set standard input and standard output to binary mode.
	return _setmode(_fileno( stdout ), winMode );
}

#else

int utils::setstdinmode(std::ios_base::openmode) {
	return 0;
}

int utils::setstdoutmode(std::ios_base::openmode) {
	return 0;
}

std::ostream &operator << (std::ostream &stream, utils::Rational &r) {
  stream << r.numerator << "/" << r.denominator;
  return stream;
}

#endif

