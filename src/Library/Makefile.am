noinst_LTLIBRARIES = libVC2.la

libVC2_la_LIBADD = \
	$(BOOST_LDFLAGS) \
	$(BOOST_SYSTEM_LIB)\
	$(BOOST_THREAD_LIB)

libVC2_la_LDFLAGS = \
	-no-undefined \
	-version-info $(VC2REFERENCE_LIBVERSION) \
	$(VC2REFERENCE_LDFLAGS) \
	-lpthread

libVC2_la_CPPFLAGS = \
	-I$(top_srcdir)/src/ \
	$(VC2REFERENCE_CFLAGS) \
	$(BOOST_CFLAGS)

libVC2_la_SOURCES = src/DataUnit.cpp src/Arrays.cpp  src/Frame.cpp  src/Picture.cpp  src/Quantisation.cpp  src/Slices.cpp  src/Utils.cpp  src/VLC.cpp  src/WaveletTransform.cpp

pkginclude_HEADERS = 

noinst_HEADERS = DataUnit.h Arrays.h Frame.h FrameResolutions.h Picture.h Quantisation.h Slices.h Utils.h VLC.h WaveletTransform.h

EXTRA_DIST =
