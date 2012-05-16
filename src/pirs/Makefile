ifeq ($(shell uname -s),Darwin)
CPPFLAGS += -march=core2 -Bstatic -L/opt/local/lib -I/opt/local/include
else
CPPFLAGS += -static
endif

#GCC_VERSION := `gcc -dumpversion`
GCC_VERSION := $(shell gcc -dumpversion)
GCC_MAJOR := $(shell echo $(GCC_VERSION) | cut -d'.' -f1)
#GCC_MINOR := $(shell echo $(GCC_VERSION) | cut -d'.' -f2)
#GCC_SUB := $(shell echo $(GCC_VERSION) | cut -d'.' -f3)
GCC4_VER := $(shell echo $(GCC_VERSION) | cut -d'.' -f2,3)

ifeq ($(GCC_MAJOR),4)
ifeq ($(shell expr 2.4 \<= $(GCC4_VER)),1)
CPPFLAGS += -march=native
endif
ifeq ($(shell expr 3.6 \<= $(GCC4_VER)),1)
CPPFLAGS += -std=c++0x
endif
endif


all: pirs

pirs: main.cpp simulate_diploid_genome.cpp simulate.cpp gzstream.cpp simulate_Illumina_reads.cpp load_file.cpp MaskQvalsByEamss.cpp
	g++ $(CPPFLAGS) -mtune=generic -O3 -o $@ $^ -lz

test: pirs
	cd test && sh pirs.sh
	ls -l test/

distclean:
	-rm test/Illumina_* test/ref22* test/*.out test/*.err test/EAMSS2_*

clean: distclean
	-rm pirs
