/*********************************************************************/
/* FrameResolutions.h                                                */
/* Author: Tim Borer                                                 */
/* This version 8th August 2010                                      */
/*                                                                   */
/* Defines a table of frame resolutions                              */
/*                                                                   */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef FRAMERESOLUTIONS_25FEB10
#define FRAMERESOLUTIONS_25FEB10

const int frameResolutions[][2] = {
//  height, width
  { 1080,   1920 },
  { 1080,   1440 },
  { 1080,   960  },
  { 720,    1280 },
  { 720,    960  },
  { 720,    640  },
  { 576,    720  },
  { 576,    704  },
  { 576,    540  },
  { 576,    360  },
  { 486,    720  },
  { 486,    704  },
  { 486,    540  },
  { 486,    360  },
  { 480,    720  },
  { 480,    704  },
  { 480,    540  },
  { 480,    360  }
};

#endif //FRAMERESOLUTIONS_25FEB10
