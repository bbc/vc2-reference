/*********************************************************************/
/* EncodeParams.h                                                    */
/* Author: Tim Borer and Galen Reich,  BBC Research & Development    */
/* This version 15th September 2020                                  */
/*                                                                   */
/* Declares getting program parameters from command line.            */
/* Copyright (c) BBC 2011-2020 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef ENCODERPARAMS_SEPT20
#define ENCODERPARAMS_SEPT20

#include <string>

#include "Picture.h"
#include "WaveletTransform.h"
#include "DataUnit.h"

enum Output {TRANSFORM, QUANTISED, INDICES, PACKAGED, STREAM, DECODED, PSNR};
enum Mode {HQ_ConstQ, HQ_CBR, LD};

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
    enum Mode mode;
    enum Output output;
    FrameRate frame_rate;
    std::string error;

    // Mode dependent
    int slice_scalar;
    int slice_prefix;
    int fragment_length;
    int compressedBytes;
    int qIndex;
};

ProgramParams getCommandLineParams(int argc, char * argv[], const char* details[]);

#endif // ENCODERPARAMS_SEPT20
