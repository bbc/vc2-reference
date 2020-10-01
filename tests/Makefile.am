check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = ../googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS =  -I$(top_srcdir)/googletest/googletest/include \
                        -I$(top_srcdir)/googletest/googletest -pthread 
libgtest_la_LDFLAGS = -pthread

MAIN_LDADD =   libgtest.la \
                $(top_srcdir)/src/Library/libVC2.la \
                $(BOOST_LDFLAGS) \
                $(BOOST_SYSTEM_LIB)\
                $(BOOST_THREAD_LIB)

MAIN_CPPFLAGS =    -I$(top_srcdir)/src/ \
                    -I$(top_srcdir)/src/Library \
                    -I$(top_srcdir)/googletest/googletest/include \
                    -I$(top_srcdir)/googletest/googletest -pthread  \
                    $(VC2REFERENCE_CFLAGS) \
	                $(BOOST_CFLAGS)

MAIN_LDFLAGS = -pthread

check_PROGRAMS = arrays quantisation dataunit utils

arrays_SOURCES = Arrays.cpp
arrays_LDADD = $(MAIN_LDADD)
arrays_LDFLAGS = $(MAIN_LDFLAGS)
arrays_CPPFLAGS = $(MAIN_CPPFLAGS)

dataunit_SOURCES = DataUnit.cpp
dataunit_LDADD = $(MAIN_LDADD)
dataunit_LDFLAGS = $(MAIN_LDFLAGS)
dataunit_CPPFLAGS = $(MAIN_CPPFLAGS)

quantisation_SOURCES = Quantisation.cpp
quantisation_LDADD = $(MAIN_LDADD)
quantisation_LDFLAGS = $(MAIN_LDFLAGS)
quantisation_CPPFLAGS = $(MAIN_CPPFLAGS)

utils_SOURCES = Utils.cpp
utils_LDADD = $(MAIN_LDADD)
utils_LDFLAGS = $(MAIN_LDFLAGS)
utils_CPPFLAGS = $(MAIN_CPPFLAGS)

TESTS = arrays quantisation dataunit utils
