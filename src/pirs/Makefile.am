AM_CPPFLAGS = @RNG_CPPFLAGS@ -DPKGDATADIR='"$(pkgdatadir)"'
AM_CXXFLAGS = -Wall -Wextra @RNG_CXXFLAGS@  @OPENMP_CXXFLAGS@
AM_LDFLAGS  = @OPENMP_LDFLAGS@

profile1dir   = $(pkgdatadir)/InDel_Profiles
profile2dir   = $(pkgdatadir)/Base-Calling_Profiles
profile3dir   = $(pkgdatadir)/GC-depth_Profiles

dist_profile1_DATA = $(srcdir)/Profiles/InDel_Profiles/*
dist_profile2_DATA = $(srcdir)/Profiles/Base-Calling_Profiles/*
dist_profile3_DATA = $(srcdir)/Profiles/GC-depth_Profiles/*

bin_PROGRAMS = pirs
pirs_SOURCES =				\
	InputStream.cpp			\
	InputStream.h			\
	OutputStream.cpp		\
	OutputStream.h			\
	main.cpp			\
	pirs.h				\
	Random.h			\
	util.cpp			\
	util.h

pirs_LDADD = -lz

if HAVE_SSE2
pirs_SOURCES +=				\
	SFMT-src-1.4/SFMT.c		\
	SFMT-src-1.4/SFMT.h		\
	SFMT-src-1.4/SFMT-params.h	\
	SFMT-src-1.4/SFMT-params19937.h	\
	SFMT-src-1.4/SFMT-common.h
else
pirs_SOURCES +=				\
	mt19937-64.cpp			\
	mt19937-64.h
endif

if ENABLE_PIRS_DIPLOID
pirs_SOURCES += pirs_diploid.cpp	
endif

if ENABLE_THREADS
pirs_SOURCES += SharedQueue.h		\
		Lock.h
endif

if ENABLE_PIRS_SIMULATE
pirs_SOURCES +=				\
	BaseCallingProfile.cpp		\
	BaseCallingProfile.h		\
	GCBiasProfile.cpp		\
	GCBiasProfile.h			\
	IndelProfile.cpp		\
	IndelProfile.h			\
	MaskQvalsByEamss.cpp		\
	MaskQvalsByEamss.h		\
	Profile.cpp			\
	Profile.h			\
	Read.h				\
	SimulationParameters.h		\
	pirs_simulate.cpp
endif
