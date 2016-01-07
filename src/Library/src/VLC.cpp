/*********************************************************************/
/* VLC.cpp                                                           */
/* Author: Tim Borer                                                 */
/* This version 9th August 2011                                      */
/*                                                                   */
/* Defines stuff related to variable length coding                   */
/* The variable length coding is interleaved exp-Golomb              */
/* 21st June 2011: Copied file from 30th May 2010                    */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include <ostream>
#include <istream>
#include <stdexcept>
#include <cstdlib> // for abs()

#include "VLC.h"

namespace {

  void encodeUnsignedVLC(unsigned int value,
                         unsigned int& nBits,
                         unsigned int& bits) {
    if (value==0) {
      nBits = 1;
      bits = 1;
    }
    else {
      value += 1;

      //Find top_bit
      int top_bit = 1;
      { unsigned int max_value = 1;
        while (value>max_value) {
          top_bit <<= 1;
          max_value <<= 1;
          max_value |= 1;
        }
      }

      nBits = 0;
      bits = 0;
      while ( top_bit>>=1 ) {
        bits <<= 2;
        if (value&top_bit) bits |= 0x1;
        nBits += 2;
      }
      bits <<= 1;
      bits |= 0x1;
      nBits += 1;
    }
  }

  const unsigned int decodeUnsignedVLC(unsigned int nBits, unsigned int bits) {
    // TO DO: Refactor with deterministic loop of (nBits-1)/2?
    unsigned int value = 1;
    unsigned int top_bit = (1<<(nBits-1));
    while ( (bits & top_bit)==0 ) {
      value <<= 1;
      top_bit >>= 1;
      if ( (bits & top_bit) ) value |= 0x1;
      top_bit >>= 1;
    }
    value -= 1;
    return value;
  }

} // end unnamed namespace

UnsignedVLC::UnsignedVLC(unsigned int value) {
  encodeUnsignedVLC(value, nBits, bits);
}

UnsignedVLC::operator const unsigned int() const {
  return decodeUnsignedVLC(nBits, bits);
}

SignedVLC::SignedVLC(int value): nBits(1), bits(1) {
  if (value) {
    encodeUnsignedVLC(abs(value), nBits, bits);
    bits <<= 1;
    if (value<0) bits |= 0x1;
    ++nBits;
  }
}

SignedVLC::operator const int() const {
  int result = 0;
  if (nBits>1) {
    result = decodeUnsignedVLC( (nBits-1), (bits>>1) );
    if ( result && (bits & 0x1) ) result *= -1;
  }
  return result;
}

namespace {

  long& cachedBits(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& cache(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& isBounded(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& bitsLeft(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }
} // end unnamed namespace

// bounded manipulator to limit number of (VLC) bits written to/read from stream
void vlc::bounded::operator()(std::ios_base& stream) const {
  isBounded(stream) = static_cast<long>(true);
  bitsLeft(stream) = static_cast<long>(bits);
}

// ostream bounded format manipulator
std::ostream& operator << (std::ostream& stream, vlc::bounded b) {
  b(stream);
  return stream;
}

//istream bounded format manipulator
std::istream& operator >> (std::istream& stream, vlc::bounded b) {
  b(stream);
  return stream;
}

// ostream unbounded format manipulator
std::istream& vlc::unbounded(std::istream& stream) {
  isBounded(stream) = static_cast<long>(false);
  return stream;
}

// ostream unbounded format manipulator
std::ostream& vlc::unbounded(std::ostream& stream) {
  isBounded(stream) = static_cast<long>(false);
  return stream;
}

namespace {

  void putBit(std::ostream& stream, bool bit) {
    if (isBounded(stream) && (bitsLeft(stream)<1)) {
        if (bit) return;
        else throw std::length_error("Attempt to write beyond end of bounded write");
    }
    unsigned char cache_byte = static_cast<unsigned char>(cache(stream));
    int cached_bits = static_cast<int>(cachedBits(stream));
    int bits_left = static_cast<int>(bitsLeft(stream));

    cache_byte <<= 1;
    if (bit) cache_byte |= 0x1;
    ++cached_bits;
    --bits_left;
    if (cached_bits==8) {
      stream.put(cache_byte);
      cached_bits=0;
    }

    cache(stream) = cache_byte;
    cachedBits(stream) = cached_bits;
    bitsLeft(stream) = bits_left;
  };

  // Put multiple bits to stream (n is number of bits to write)
  void putBits(std::ostream& stream, unsigned int n, unsigned int value) {
    while (n>0) {
      --n;
      putBit(stream, (value>>n)&0x1 );
    }
  }

  bool getBit(std::istream& stream) {
    if (isBounded(stream) && (bitsLeft(stream)<1)) {
      return true;
    }
    unsigned char cache_byte = static_cast<unsigned char>(cache(stream));
    int cached_bits = static_cast<int>(cachedBits(stream));
    int bits_left = static_cast<int>(bitsLeft(stream));
    
    if (cached_bits == 0) {
      cache_byte = stream.get();
      cached_bits = 8;
    }
    --cached_bits;
    --bits_left;

    cache(stream) = cache_byte;
    cachedBits(stream) = cached_bits;
    bitsLeft(stream) = bits_left;

    return ((cache_byte>>cached_bits) & 0x1);
  }

  // Get multiple bits from stream (n is number of bits to read)
  unsigned int getBits(std::istream& stream, unsigned int n) {
    unsigned int value = 0;
    while (n>0) {
      value <<= 1;
      if (getBit(stream)) value |= 0x1;
      --n;
    }
    return value;
  }

}

std::ostream& operator << (std::ostream& stream, Boolean bit) {
  putBit(stream, bit);
  return stream;
}

std::istream& operator >> (std::istream& stream, Boolean& bit) {
  bit = getBit(stream);
  return stream;
}

// ostream flush writes zero bits to end of bounded stream
// Note: end of bounded stream may not be byte aligned
std::ostream& vlc::flush(std::ostream& stream) {
  if (isBounded(stream)) {
    while (bitsLeft(stream)>0) putBit(stream, false);
  }
  return stream;
}

// istream flush moves read pointer to end of bounded stream
// Note: end of bounded stream may not be byte aligned
std::istream& vlc::flush(std::istream& stream) {
  if (isBounded(stream)) {
    while (bitsLeft(stream)>0) getBit(stream);
  }
  return stream;
}

// ostream align writes zero bits to start of next byte
std::ostream& vlc::align(std::ostream& stream) {
  isBounded(stream) = false;
  while (cachedBits(stream)) putBit(stream, false);
  return stream;
}

// istream align moves read pointer to start of next byte
std::istream& vlc::align(std::istream& stream) {
  isBounded(stream) = false;
  while (cachedBits(stream)) getBit(stream);
  return stream;
}

Bits::Bits(unsigned int no_of_bits, unsigned int value):
    nBits(no_of_bits), bits(value) {
  if ((value & ((1 << no_of_bits) - 1)) != value) {
    throw std::logic_error("Bits: value bigger than specified number of bits");
  }
}

std::ostream& operator << (std::ostream& stream, Bits b) {
  putBits(stream, b.bitCount(), b);
  return stream;
}

std::istream& operator >> (std::istream& stream, Bits& b) {
  b = getBits(stream, b.bitCount());
  return stream;
}

// ostream inserter for unsigned interleaved expGolomb
std::ostream& operator << (std::ostream& stream, UnsignedVLC g) {
  putBits(stream, g.numOfBits(), g.code());
  return stream;
}

// istream extractor for unsigned interleaved expGolomb
std::istream& operator >> (std::istream& stream, UnsignedVLC& g) {
  unsigned int numOfBits=0, code=0;
  while (!getBit(stream)) {
    code <<= 2;
    if (getBit(stream)) code |= 0x1;
    numOfBits += 2;
  }
  code <<= 1;
  code |= 0x1;
  numOfBits += 1;
  g = UnsignedVLC(numOfBits, code);
  return stream;
}

// ostream inserter for signed interleaved expGolomb
std::ostream& operator << (std::ostream& stream, SignedVLC g) {
  putBits(stream, g.numOfBits(), g.code());
  return stream;
}

// istream extractor for signed interleaved expGolomb
std::istream& operator >> (std::istream& stream, SignedVLC& g) {
  UnsignedVLC unsignedPart;
  stream >> unsignedPart;
  if (unsignedPart.numOfBits()==1) {
    g = SignedVLC(unsignedPart.numOfBits(), unsignedPart.code());
  }
  else {
    unsigned int signBit=0;
    if (getBit(stream)) signBit = 0x1;
    g = SignedVLC(unsignedPart.numOfBits()+1,
                 (unsignedPart.code()<<1) | signBit);
  }
  return stream;
}

Bytes::Bytes(unsigned int no_of_bytes, unsigned long value):
    nBytes(no_of_bytes), bytes(value) {
  if ((value & (0xFFFFFFFF >> ((4 - no_of_bytes) << 3))) != value) {
      throw std::logic_error("Bytes: value bigger than specified number of bytes");
  }
}

std::ostream& operator << (std::ostream& stream, Bytes b) {
  stream << vlc::align;
  int bytesLeft = b.byteCount();
  const unsigned long value = b;
  while (bytesLeft>0) {
    --bytesLeft;
    stream.put(static_cast<unsigned char>(value>>(bytesLeft<<3)));
  }
  return stream;
}

std::istream& operator >> (std::istream& stream, Bytes& b) {
  stream >> vlc::align;
  int bytesLeft = b.byteCount();
  unsigned long value = 0;
  while (bytesLeft>0) {
      value <<= 8;
      value |= stream.get();
      --bytesLeft;
  }
  b = value;
  return stream;
}
