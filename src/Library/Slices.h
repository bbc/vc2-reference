/*********************************************************************/
/* Slices.h                                                          */
/* Author: Tim Borer                                                 */
/* This version 8th August 2011                                      */
/*                                                                   */
/* Defines stuff releating to slices in both low delay and           */
/* high quality profile                                              */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef SLICES_24JUNE11
#define SLICES_24JUNE11

#include "Arrays.h"
#include "Picture.h"

// This slice_bytes returns the actual number of bytes for a slice at specific co-ordinates
const int slice_bytes(int v, int h, // Slice co-ordinates
                     const int ySlices, const int xSlices, // Number of slices
                     const int sliceBytesNumerator, const int sliceBytesDenominator);

// And this one fills an arrays with the number of bytes for corresponding slices (Scalar should be 1 for LD)
const Array2D slice_bytes(const int ySlices, const int xSlices, const int totalBytes, const int scalar);

// luma_slice_bits returns the number of bits in an LD luma slice after VLC coding
const int luma_slice_bits(const Array2D& lumaSlice, const char waveletDepth);

// chroma_slice_bits returns the number of bits in an LD luma slice after VLC coding
const int chroma_slice_bits(const Array2D& uSlice, const Array2D& vSlice, const char waveletDepth);

// Returns the number of bytes in a component HQ slice after VLC coding (no DC prediction)
const int component_slice_bytes(const Array2D& componentSlice, const char waveletDepth, const int scalar);

// Define a state machine for quantising slices
class SliceQuantiser {
  public:
    SliceQuantiser(const Array2D& coefficients,
                   int vSlices, int hSlices,
                   const Array1D& quantMatrix);
    const int row() const {return v;} //Get row number for slice
    const int column() const {return h;} //Get column number for slice
    const bool next_slice(); //Go to next slice in raster order (returns false if no more slices)
    virtual const Array2D& quantise_slice(int qIndex) = 0;
    virtual ~SliceQuantiser() {}
  protected:
    const int ySlices; //Number of rows
    const int xSlices; //Number of columns
    const int coeffsHeight;
    const int coeffsWidth;
    const int sliceHeight;
    const int sliceWidth;
    const int numberOfSubbands;
    const int waveletDepth;
    const int transformSize; // 2**waveletDepth
    int v; //Current row
    int h; //Current column
    Array2D qSlice; // Slice array of quantised coeffs to be returned by "quantise_slice"
  private:
    SliceQuantiser(const SliceQuantiser&); //No copying
    SliceQuantiser& operator=(const SliceQuantiser&); //No assignment
};

//**** Slice IO declarations ****//

struct Slices { 
  Slices(const PictureArray& yuvSlices, const int waveletDepth, const Array2D& qIndices);
  Slices(const PictureFormat& pictureFormat, int waveletDepth,
         int ySlices, int xSlices);
  PictureArray yuvSlices;
  const int waveletDepth;
  Array2D qIndices;

  int nSlices() { return yuvSlices.shape()[0]*yuvSlices.shape()[1]; }
};

std::ostream& operator << (std::ostream& stream, const Slices& s);

std::istream& operator >> (std::istream& stream, Slices& s);

struct Slice {
    Slice(const Picture& p, int d, int i):
      yuvSlice(p), waveletDepth(d), qIndex(i) {};
    Slice(const PictureFormat& f, int d):
      yuvSlice(f), waveletDepth(d) {};
    Picture yuvSlice;
    const int waveletDepth;
    int qIndex;
};

std::ostream& operator << (std::ostream& stream, const Slice& s);

std::istream& operator >> (std::istream& stream, Slice& s);

class Fragment;

namespace sliceio {

  enum SliceIOMode {UNKNOWN, LD, HQVBR, HQCBR};

  SliceIOMode &sliceIOMode(std::ios_base& stream);

  class lowDelay {
    public:
      lowDelay(const Array2D& b): bytes(b) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const Array2D& bytes;
  };

  class highQualityCBR {
    public:
      highQualityCBR(const Array2D& b, const int p, const int s): bytes(b), prefix(p), scalar(s) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const Array2D& bytes;
      const int prefix;
      const int scalar;
  };

  class highQualityVBR {
    public:
      highQualityVBR(const int p, const int s): prefix(p), scalar(s) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const int prefix;
      const int scalar;
  };

  // ostream format manipulator to set the size of a single slice
  class setBytes {
  public:
    setBytes(const int b): bytes(b) {}; 
    void operator () (std::ios_base& stream) const;
  private:
    const int bytes;
  };

  const Array2D *SliceSizes(std::ios_base& stream);

  class ExpectedSlicesForFragment {
  public:
    ExpectedSlicesForFragment(Fragment &frag);
    void operator () (std::ios_base& stream) const;
  private:
    const int n_slices;
    const int slice_offset_x;
    const int slice_offset_y;
  };
} // end namespace sliceio

// istream format manipulator to set slice offsets when reading fragments
std::istream& operator >> (std::istream& stream, sliceio::ExpectedSlicesForFragment esf);

// ostream format manipulator to set the size of a single slice
std::ostream& operator << (std::ostream& stream, sliceio::setBytes arg);

// istream format manipulator to set the size of a single slice
std::istream& operator >> (std::istream& stream, sliceio::setBytes arg);

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::lowDelay arg);

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::lowDelay arg);

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::highQualityCBR arg);

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::highQualityCBR arg);

// ostream low delay format manipulator
std::ostream& operator << (std::ostream& stream, sliceio::highQualityVBR arg);

// istream low delay format manipulator
std::istream& operator >> (std::istream& stream, sliceio::highQualityVBR arg);

#endif //SLICES_24JUNE11
