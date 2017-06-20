/*********************************************************************/
/* Slices.cpp                                                        */
/* Author: Tim Borer                                                 */
/* This version 8th August 2011                                       */
/*                                                                   */
/* Defines stuff relating to slices in high quality profile          */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include <iostream> //For cin, cout, cerr

#include "Slices.h"
#include "WaveletTransform.h"
#include "VLC.h"
#include "Utils.h"
#include "DataUnit.h"

const int slice_bytes(int v, int h, // Slice co-ordinates
                     const int ySlices, const int xSlices, // Number of slices
                     const int sliceBytesNumerator, const int sliceBytesDenominator) { //Slice size
  const long long sliceNumber = v*xSlices+h;
  long long bytes;
  bytes = ((sliceNumber+1)*sliceBytesNumerator)/sliceBytesDenominator;
  bytes -= (sliceNumber*sliceBytesNumerator)/sliceBytesDenominator;
  return static_cast<const int>(bytes);
}

const Array2D slice_bytes(const int ySlices, const int xSlices, const int totalBytes, const int scalar) {
  const utils::Rational rationalBytes = utils::rationalise(totalBytes/scalar - 4*(ySlices*xSlices), (ySlices*xSlices));
  const int sliceBytesNumerator = rationalBytes.numerator;
  const int sliceBytesDenominator = rationalBytes.denominator;
  const int ratio = sliceBytesNumerator/sliceBytesDenominator;
  const int remainder = sliceBytesNumerator - (ratio*sliceBytesDenominator);
  int residue = 0;
  Array2D bytes(extents[ySlices][xSlices]);
  for (int v=0; v<ySlices; ++v) {
    for (int h=0; h<xSlices; ++h) {
      residue += remainder;
      if (residue<sliceBytesDenominator) {
        bytes[v][h] = ratio*scalar + 4;
      }
      else {
        bytes[v][h] = ((ratio+1)*scalar) + 4;
        residue -= sliceBytesDenominator;
      }
    }
  }
  return bytes;
}

const int luma_slice_bits(const Array2D& lumaSlice, const char waveletDepth) {
  const int numberOfSubbands = 3*waveletDepth+1;
  const BlockVector subbands = split_into_subbands(lumaSlice, waveletDepth);
  int count = 0;
  int gross = 0;
  for (int band=0; band<numberOfSubbands; ++band) {
    const Array2D& subband = subbands[band];
    const int height = subband.shape()[0];
    const int width = subband.shape()[1];
    for (int y=0; y<height; ++y) {
      for (int x=0; x<width; ++x) {
        const int numBits = SignedVLC(subband[y][x]).numOfBits();
        gross += numBits;
        if (numBits>1) count=gross;
      }
    }
  }
  return count;
}

const int chroma_slice_bits(const Array2D& uSlice, const Array2D& vSlice, const char waveletDepth) {
  const int numberOfSubbands = 3*waveletDepth+1;
  const BlockVector uSubbands = split_into_subbands(uSlice, waveletDepth);
  const BlockVector vSubbands = split_into_subbands(vSlice, waveletDepth);
  int count = 0;
  int gross = 0;
  for (int band=0; band<numberOfSubbands; ++band) {
    const Array2D& uSubband = uSubbands[band];
    const Array2D& vSubband = vSubbands[band];
    // TO DO: Check uSubband & vSubband have the same shape?
    const int height = uSubband.shape()[0];
    const int width = uSubband.shape()[1];
    for (int y=0; y<height; ++y) {
      for (int x=0; x<width; ++x) {
        int numBits = SignedVLC(uSubband[y][x]).numOfBits();
        gross += numBits;
        if (numBits>1) count=gross;
        numBits = SignedVLC(vSubband[y][x]).numOfBits();
        gross += numBits;
        if (numBits>1) count=gross;
      }
    }
  }
  return count;
}

const int component_slice_bytes(const Array2D& slice, const char waveletDepth, const int scalar) {
  const int numberOfSubbands = 3*waveletDepth+1;
  const BlockVector subbands = split_into_subbands(slice, waveletDepth);
  int count = 0;
  int gross = 0;
  for (int band=0; band<numberOfSubbands; ++band) {
    const Array2D& subband = subbands[band];
    const int height = subband.shape()[0];
    const int width = subband.shape()[1];
    for (int y=0; y<height; ++y) {
      for (int x=0; x<width; ++x) {
        const int numBits = SignedVLC(subband[y][x]).numOfBits();
        gross += numBits;
        if (numBits>1) count=gross;
      }
    }
  }
  return (((count+7)/8 + scalar - 1)/scalar)*scalar; // return whole number of scalar byte units
}

SliceQuantiser::SliceQuantiser(const Array2D& coefficients,
                               int vSlices, int hSlices,
                               const Array1D& quantMatrix):
  ySlices(vSlices),
  xSlices(hSlices),
  coeffsHeight(coefficients.shape()[0]),
  coeffsWidth(coefficients.shape()[1]),
  sliceHeight(coeffsHeight/vSlices),
  sliceWidth(coeffsWidth/hSlices),
  numberOfSubbands(quantMatrix.size()),
  waveletDepth((numberOfSubbands-1)/3),
  transformSize(utils::pow(2, waveletDepth)),
  v(0), h(0)
{
  qSlice.resize(extents[sliceHeight][sliceWidth]);
}

const bool SliceQuantiser::next_slice() {
  if (h<(xSlices-1)) { //Next column
    ++h;
    return true; }
  else if (v<(ySlices-1)){ //Next row
    h=0;
    ++v;
    return true; }
  return false; //No more slices
}

//**** IO functions ****//

namespace {

  // Choice of UNKNOWN, LD, HQVBR, HQCBR (see Slices.h)
  // Note: returns true if IO format has been set
  long& slice_IO_format(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& slice_sizes(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& slice_scalar(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

   long& slice_prefix(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& single_slice_size(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& num_slices(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& slice_offset_x(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& slice_offset_y(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  std::ostream& LDSliceIO(std::ostream& stream, const Slice& s) {
    //Get slice size from the stream
    const int sliceSize = static_cast<int>(single_slice_size(stream));

    const BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    const BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    const BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);

    stream << Bits(7, s.qIndex);

    const int yBits = luma_slice_bits(s.yuvSlice.y(), s.waveletDepth);
    const int uvSplitBits = utils::intlog2(8*sliceSize-7);
    const int uvBits = 8*sliceSize - 7 - uvSplitBits - yBits;
    stream << Bits(uvSplitBits, yBits);

    const int numberOfSubbands = 3*s.waveletDepth+1;
    stream << vlc::bounded(yBits);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& ySubband = ySliceSubbands[band];
      const int height = ySubband.shape()[0];
      const int width = ySubband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(ySubband[y][x]);
        }
      }
    }
    stream << vlc::flush;

    stream << vlc::bounded(uvBits);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& uSubband = uSliceSubbands[band];
      const Array2D& vSubband = vSliceSubbands[band];
      const int height = uSubband.shape()[0];
      const int width = uSubband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(uSubband[y][x]);
          stream << SignedVLC(vSubband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;
    return stream;
  }

  std::istream& LDSliceIO(std::istream& stream, Slice& s) {
    const int sliceSize = static_cast<int>(single_slice_size(stream));

    BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);

    Bits q(7);
    stream >> q;
    s.qIndex = q;
    const int uvSplitBits = utils::intlog2(8*sliceSize-7);
    int yBits;
    Bits yb(uvSplitBits);
    stream >> yb;
    yBits = yb;
    const int uvBits = 8*sliceSize - 7 - uvSplitBits - yBits;

    const int numberOfSubbands = 3*s.waveletDepth+1;
    stream >> vlc::bounded(yBits);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& ySubband = ySliceSubbands[band];
      const int height = ySubband.shape()[0];
      const int width = ySubband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          ySubband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush;

    stream >> vlc::bounded(uvBits);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& uSubband = uSliceSubbands[band];
      Array2D& vSubband = vSliceSubbands[band];
      const int height = uSubband.shape()[0];
      const int width = uSubband.shape()[1];
      // TO DO: Check u and v subbands have the same shape?
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          uSubband[y][x] = inVLC;
          stream >> inVLC;
          vSubband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;

    s.yuvSlice.y(merge_subbands(ySliceSubbands));
    s.yuvSlice.c1(merge_subbands(uSliceSubbands));
    s.yuvSlice.c2(merge_subbands(vSliceSubbands));

    return stream;
  }

  std::ostream& HQSliceIO_CBR(std::ostream& stream, const Slice& s) {
    const int sliceSize = static_cast<int>(single_slice_size(stream));

    const BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    const BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    const BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);
    const int numberOfSubbands = 3*s.waveletDepth+1;
    const int scalar = slice_scalar(stream);
    const int prefix = slice_prefix(stream);

    for (int n = 0; n < prefix; n++) {
      stream << Bytes(1, 0x00);
    }
    
    stream << Bytes(1, s.qIndex);

    // Output first (y/luma) component
    const int yBytes = component_slice_bytes(s.yuvSlice.y(), s.waveletDepth, scalar);
    stream << Bytes(1, yBytes/scalar);
    stream << vlc::bounded(8*yBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = ySliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;

    // Output secomd (u/c1/chroma) component
    const int uBytes = component_slice_bytes(s.yuvSlice.c1(), s.waveletDepth, scalar);
    stream << Bytes(1, uBytes/scalar);
    stream << vlc::bounded(8*uBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = uSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;
    
    // Output third (v/c2/chroma) component
    // Calculate bytes left for u, and throw if too few bytes avaiable
    const int vBytes = sliceSize - 4 - yBytes - uBytes;
    if (vBytes < component_slice_bytes(s.yuvSlice.c2(), s.waveletDepth, scalar) ) {
      throw std::logic_error("SliceIO, HQ CBR mode: Too many bytes for the slice");
    }
    stream << Bytes(1, vBytes/scalar);
    stream << vlc::bounded(8*vBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = vSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;

    return stream;
  }

  std::istream& HQSliceIO_CBR(std::istream& stream, Slice& s) {
    const int sliceSize = static_cast<int>(single_slice_size(stream));

    BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);
    const int numberOfSubbands = 3*s.waveletDepth+1;
    const int scalar = slice_scalar(stream);
    const int prefix = slice_prefix(stream);

    Bytes bytes(1);

    for (int n = 0; n < prefix; n++) {
      Bytes z(1);
      stream >> z;
    }

    Bytes q(1);
    stream >> q;
    s.qIndex = q;

    // Input first (y/luma) component
    stream >> bytes;
    const int yBytes = ((int)bytes)*scalar;
    stream >> vlc::bounded(8*yBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = ySliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;

    // Input second (u/c1/chroma) component
    stream >> bytes;
    const int uBytes = ((int)bytes)*scalar;
    stream >> vlc::bounded(8*uBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = uSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;
    
    // Input third (v/c2/chroma) component
    stream >> bytes;
    // Calculate bytes left for u, and throw if number of bytes read from stream disagrees
    const int vBytes = sliceSize - 4 - yBytes - uBytes;
    if (vBytes != static_cast<const int>(bytes) )
      throw std::logic_error("SliceIO, HQ CBR mode: Wrong number of bytes for a slice");
    stream >> vlc::bounded(8*vBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = vSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;

    s.yuvSlice.y(merge_subbands(ySliceSubbands));
    s.yuvSlice.c1(merge_subbands(uSliceSubbands));
    s.yuvSlice.c2(merge_subbands(vSliceSubbands));

    return stream;
  }

  std::ostream& HQSliceIO_VBR(std::ostream& stream, const Slice& s) {
    const BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    const BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    const BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);
    const int numberOfSubbands = 3*s.waveletDepth+1;

    const int scalar = slice_scalar(stream);
    const int prefix = slice_prefix(stream);

    for (int n = 0; n < prefix; n++) {
      stream << Bytes(1, 0x00);
    }

    stream << Bytes(1, s.qIndex);

    // Output first (y/luma) component
    const int yBytes = component_slice_bytes(s.yuvSlice.y(), s.waveletDepth, scalar);
    stream << Bytes(1, yBytes/scalar);
    stream << vlc::bounded(8*yBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = ySliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;

    // Output secomd (u/c1/chroma) component
    const int uBytes = component_slice_bytes(s.yuvSlice.c1(), s.waveletDepth, scalar);
    stream << Bytes(1, uBytes/scalar);
    stream << vlc::bounded(8*uBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = uSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;
    
    // Output third (v/c2/chroma) component
    const int vBytes = component_slice_bytes(s.yuvSlice.c2(), s.waveletDepth, scalar);
    stream << Bytes(1, vBytes/scalar);
    stream << vlc::bounded(8*vBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      const Array2D& subband = vSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream << SignedVLC(subband[y][x]);
        }
      }
    }
    stream << vlc::flush << vlc::align;

    return stream;
  }

  std::istream& HQSliceIO_VBR(std::istream& stream, Slice& s) {
    BlockVector ySliceSubbands = split_into_subbands(s.yuvSlice.y(), s.waveletDepth);
    BlockVector uSliceSubbands = split_into_subbands(s.yuvSlice.c1(), s.waveletDepth);
    BlockVector vSliceSubbands = split_into_subbands(s.yuvSlice.c2(), s.waveletDepth);
    const int numberOfSubbands = 3*s.waveletDepth+1;
    const int scalar = slice_scalar(stream);
    const int prefix = slice_prefix(stream);
    Bytes bytes(1);

    for (int n = 0; n < prefix; n++) {
      Bytes z(1);
      stream >> z;
    }

    Bytes q(1);
    stream >> q;
    s.qIndex = q;

    // Input first (y/luma) component
    stream >> bytes;
    const int yBytes = ((int)bytes)*scalar;
    stream >> vlc::bounded(8*yBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = ySliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;

    // Input second (u/c1/chroma) component
    stream >> bytes;
    const int uBytes = ((int)bytes)*scalar;
    stream >> vlc::bounded(8*uBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = uSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;
    
    // Input third (v/c2/chroma) component
    stream >> bytes;
    const int vBytes = ((int)bytes)*scalar;
    stream >> vlc::bounded(8*vBytes);
    for (int band=0; band<numberOfSubbands; ++band) {
      Array2D& subband = vSliceSubbands[band];
      const int height = subband.shape()[0];
      const int width = subband.shape()[1];
      SignedVLC inVLC;
      for (int y=0; y<height; ++y) {
        for (int x=0; x<width; ++x) {
          stream >> inVLC;
          subband[y][x] = inVLC;
        }
      }
    }
    stream >> vlc::flush >> vlc::align;

    s.yuvSlice.y(merge_subbands(ySliceSubbands));
    s.yuvSlice.c1(merge_subbands(uSliceSubbands));
    s.yuvSlice.c2(merge_subbands(vSliceSubbands));

    return stream;
  }

} // End unnamed namespace

sliceio::SliceIOMode &sliceio::sliceIOMode(std::ios_base& stream) {
  return reinterpret_cast<sliceio::SliceIOMode &>(slice_IO_format(stream));
}

Slices::Slices(const PictureArray& s, const int d, const Array2D& i):
  yuvSlices(s), waveletDepth(d), qIndices(i) {
};

Slices::Slices(const PictureFormat& picFormat, int d,int ySlices, int xSlices):
    waveletDepth(d) {
  const int lumaSliceHeight = picFormat.lumaHeight()/ySlices;
  const int lumaSliceWidth = picFormat.lumaWidth()/xSlices;
  const int chromaSliceHeight = picFormat.chromaHeight()/ySlices;
  const int chromaSliceWidth = picFormat.chromaWidth()/xSlices;
  const PictureFormat sliceFormat(lumaSliceHeight, lumaSliceWidth,
                                  chromaSliceHeight, chromaSliceWidth,
                                  picFormat.chromaFormat());
  const Shape2D shape = {{ySlices, xSlices}};
  yuvSlices = PictureArray(shape);
  for (int v=0; v<ySlices; ++v) {
    for (int h=0; h<xSlices; ++h) {
      yuvSlices[v][h]=Picture(sliceFormat);
    }
  }
  qIndices = Array2D(shape);
};

#include <iostream>

std::ostream& operator << (std::ostream& stream, const Slices& s) {
  const Array2D& bytes = *reinterpret_cast<const Array2D *>(slice_sizes(stream));
  const bool bytes_valid = (slice_sizes(stream)!=0);
  const PictureArray& yuvSlices = s.yuvSlices;
  const Array2D& qIndices = s.qIndices;
  const int waveletDepth = s.waveletDepth;
  const int ySlices = yuvSlices.shape()[0];
  const int xSlices = yuvSlices.shape()[1];
  for (int v=0; v<ySlices; ++v) {
    for (int h=0; h<xSlices; ++h) {
      if (bytes_valid) stream << sliceio::setBytes(bytes[v][h]);
      stream << Slice(yuvSlices[v][h], waveletDepth, qIndices[v][h]);
    }
  }
  return stream;
}

std::istream& operator >> (std::istream& stream, Slices& s) {
  Array2D& bytes = *reinterpret_cast<Array2D *>(slice_sizes(stream));
  const bool bytes_valid = (slice_sizes(stream)!=0);
  PictureArray& yuvSlices = s.yuvSlices;
  Array2D& qIndices = s.qIndices;
  const int waveletDepth = s.waveletDepth;
  const int ySlices = yuvSlices.shape()[0];
  const int xSlices = yuvSlices.shape()[1];
  const int n_slices = num_slices(stream);
  int n = 0;
  int h = slice_offset_x(stream);
  int v = slice_offset_y(stream);
  while (n_slices == 0 || n < n_slices) {
    Slice inSlice(yuvSlices[v][h].format(), waveletDepth);
    if (bytes_valid) stream >> sliceio::setBytes(bytes[v][h]);
    stream >> inSlice;
    yuvSlices[v][h].y(inSlice.yuvSlice.y());
    yuvSlices[v][h].c1(inSlice.yuvSlice.c1());
    yuvSlices[v][h].c2(inSlice.yuvSlice.c2());
    qIndices[v][h] = inSlice.qIndex;

    n++;
    h++;
    if (h == xSlices) {
      h = 0;
      v++;
      if (v == ySlices) {
        break;
      }
    }
  }
  return stream;
}

std::ostream& operator << (std::ostream& stream, const Slice& s) {
  if (!slice_IO_format(stream))
    throw std::logic_error("SliceIO: Output Format not set");
  switch (static_cast<sliceio::SliceIOMode>(slice_IO_format(stream))) {
    case sliceio::LD:
      return LDSliceIO(stream, s);
      break;
    case sliceio::HQVBR:
      return HQSliceIO_VBR(stream, s);
      break;
    case sliceio::HQCBR:
      return HQSliceIO_CBR(stream, s);
      break;
    default:
      throw std::logic_error("SliceIO: Unknown Output Format");
  }
}

std::istream& operator >> (std::istream& stream, Slice& s) {
  if (!slice_IO_format(stream))
    throw std::logic_error("SliceIO: Input Format not set");
  switch (static_cast<sliceio::SliceIOMode>(slice_IO_format(stream))) {
    case sliceio::LD:
      return LDSliceIO(stream, s);
      break;
    case sliceio::HQCBR:
      return HQSliceIO_CBR(stream, s);
      break;
    case sliceio::HQVBR:
      return HQSliceIO_VBR(stream, s);
      break;
    default:
      throw std::logic_error("SliceIO: Unknown Input Format");
  }
}

const Array2D *sliceio::SliceSizes(std::ios_base& stream) {
  return reinterpret_cast<const Array2D *>(slice_sizes(stream));
}

sliceio::ExpectedSlicesForFragment::ExpectedSlicesForFragment(Fragment &frag)
  : n_slices(frag.n_slices())
  , slice_offset_x(frag.slice_offset_x())
  , slice_offset_y(frag.slice_offset_y()) {}

void sliceio::ExpectedSlicesForFragment::operator()(std::ios_base &stream) const {
  num_slices(stream) = this->n_slices;
  ::slice_offset_x(stream) = this->slice_offset_x;
  ::slice_offset_y(stream) = this->slice_offset_y;
}

std::istream& operator >> (std::istream& stream, sliceio::ExpectedSlicesForFragment esf) {
  esf(stream);
  return stream;
}

// IO format manipulator to set the low delay IO format
void sliceio::lowDelay::operator()(std::ios_base& stream) const {
  slice_IO_format(stream) = static_cast<long>(LD);
  slice_sizes(stream) = reinterpret_cast<long>(&bytes);
}

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::lowDelay arg) {
  arg(stream);
  return stream;
}

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::lowDelay arg) {
  arg(stream);
  return stream;
}

// IO format manipulator to set the high quality CBR IO format
void sliceio::highQualityCBR::operator()(std::ios_base& stream) const {
  slice_IO_format(stream) = static_cast<long>(HQCBR);
  slice_sizes(stream) = reinterpret_cast<long>(&bytes);
  slice_prefix(stream) = static_cast<long>(prefix);
  slice_scalar(stream) = static_cast<long>(scalar);
}

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::highQualityCBR arg) {
  arg(stream);
  return stream;
}

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::highQualityCBR arg) {
  arg(stream);
  return stream;
}

// IO format manipulator to set the high quality CBR IO format
void sliceio::highQualityVBR::operator()(std::ios_base& stream) const {
  slice_IO_format(stream) = static_cast<long>(HQVBR);
  slice_prefix(stream) = static_cast<long>(prefix);
  slice_scalar(stream) = static_cast<long>(scalar);
}

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::highQualityVBR arg) {
  arg(stream);
  return stream;
}

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::highQualityVBR arg) {
  arg(stream);
  return stream;
}

// IO format manipulator to set the size of a single slice
void sliceio::setBytes::operator()(std::ios_base& stream) const {
  single_slice_size(stream) = bytes;
}

// ostream format manipulator to set the size of a single slice
std::ostream& operator << (std::ostream& stream, sliceio::setBytes arg) {
  arg(stream);
  return stream;
}

// istream format manipulator to set the size of a single slice
std::istream& operator >> (std::istream& stream, sliceio::setBytes arg) {
  arg(stream);
  return stream;
}
