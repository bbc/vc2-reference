/*********************************************************************/
/* VLC.h                                                             */
/* Author: Tim Borer                                                 */
/* This version 26th June 2011                                       */
/*                                                                   */
/* Declares stuff related to variable length coding (exp-Golomb)     */
/* 21st June 2011: Copied version from 30th May 2010                 */
/* 26th June 2011: Made number of bytes in "Bytes" constant          */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef VLC_14MAY10
#define VLC_14MAY10

#include <iosfwd>

class UnsignedVLC {
  public:
    UnsignedVLC() {}
    UnsignedVLC(unsigned int value);
    UnsignedVLC(unsigned int numOfBits, unsigned int code):
      nBits(numOfBits), bits(code) {};
    const unsigned int numOfBits() const {return nBits;}
    const unsigned int code() const {return bits;}
    operator const unsigned int() const;
  public:
    unsigned int nBits;
    unsigned int bits;
};

class SignedVLC {
  public:
    SignedVLC() {}
    SignedVLC(int value);
    SignedVLC(unsigned int numOfBits, unsigned int code):
      nBits(numOfBits), bits(code) {};
    const unsigned int numOfBits() const {return nBits;}
    const unsigned int code() const {return bits;}
    operator const int() const;
  public:
    unsigned int nBits;
    unsigned int bits;
};

// ostream inserter for unsigned interleaved expGolomb
std::ostream& operator << (std::ostream& stream, UnsignedVLC);

// istream extractor for unsigned interleaved expGolomb
std::istream& operator >> (std::istream& stream, UnsignedVLC&);

// ostream inserter for signed interleaved expGolomb
std::ostream& operator << (std::ostream& stream, SignedVLC);

// istream extractor for signed interleaved expGolomb
std::istream& operator >> (std::istream& stream, SignedVLC&);

class Boolean {
  public:
    Boolean() {}
    Boolean(bool arg): bit(arg) {}
    operator bool() const {return bit;}
  private:
    bool bit;
};

std::ostream& operator << (std::ostream& stream, Boolean);

std::istream& operator >> (std::istream& stream, Boolean&);

class Bits {
  public:
    explicit Bits(unsigned int no_of_bits): nBits(no_of_bits) {}
    Bits(unsigned int no_of_bits, unsigned int value);
    Bits& operator=(unsigned int arg) {
      bits=arg;
      return *this; }
    operator unsigned int() const {return bits;}
    const unsigned int bitCount() const {return nBits;}
  private:
    const unsigned int nBits;
    unsigned int bits;
};

std::ostream& operator << (std::ostream& stream, Bits b);

std::istream& operator >> (std::istream& stream, Bits& b);

class Bytes {
  public:
    explicit Bytes(unsigned int no_of_bytes): nBytes(no_of_bytes) {}
    Bytes(unsigned int no_of_bytes, unsigned long value);
    Bytes& operator=(unsigned long arg) {
      bytes=arg;
      return *this; }
    operator unsigned long() const {return bytes;}
    const unsigned int byteCount() const {return nBytes;}
  private:
    const unsigned int nBytes;
    unsigned long bytes;
};

std::ostream& operator << (std::ostream& stream, Bytes b);

std::istream& operator >> (std::istream& stream, Bytes& b);

namespace vlc {

  class bounded {
    public:
      bounded(int maxBits): bits(maxBits) {}; 
      // Set maximum bits more to be written to stream
      void operator () (std::ios_base& stream) const;
    private:
      int bits;
  };

  // ostream unbounded format manipulator
  std::ostream& unbounded(std::ostream& stream);

  // ostream unbounded format manipulator
  std::istream& unbounded(std::istream& stream);

  // ostream flush writes zero bits to end of bounded stream
  // Note: end of bounded stream may not be byte aligned
  std::ostream& flush(std::ostream& stream);

  // istream flush moves read pointer to end of bounded stream
  // Note: end of bounded stream may not be byte aligned
  std::istream& flush(std::istream& stream);

  // ostream align moves read pinter to start of next byte
  std::ostream& align(std::ostream& stream);

  // istream align moves write pinter to start of next byte
  std::istream& align(std::istream& stream);

} // end namespace vlc

// ostream bounded format manipulator
std::ostream& operator << (std::ostream& stream, vlc::bounded b);

// istream bounded format manipulator
std::istream& operator >> (std::istream& stream, vlc::bounded b);

#endif //VLC_14MAY10
