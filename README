VC-2 Reference Encoder and Decoder
----------------------------------

Copyright (C) Tim Borer and James Weaver 2010-2015, British
Broadcasting Corporation.
< james.barrett@bbc.co.uk >


This repository contains a SMPTE 2042-1 VC-2 reference encoder and
decoder. It can be compiled using autotools on Linux or Windows and
includes the following executables once compiled:

 o EncodeLD -- an encoder for the Low Delay (LD) profile of VC-2.
 o EncodeHQ-CBR -- an encoder for the High Quality (HQ) profile of
   VC-2 which encodes at a constant bit rate.
 o EncodeHQ-ConstQ -- an encoder for the High Quality (HQ) profile of
   VC-2 which encodes with a constant quantiser value.
 o DecodeStream -- a decoder which will decode a VC-2 compliant stream
   which complies with the LD or HQ profiles.

In addition code is included for decoders which take in the compressed
bytes of a VC-2 frame without any surrounding headers. These are not
compiled by default but can be enabled with the
--enable-other-decoders flag to configure.

Additional help on each executable will be printed if it is run with
the --help parameter.
