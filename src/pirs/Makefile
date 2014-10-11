CFLAGS = -pipe -mtune=generic -msse2 -fopenmp -DSFMT_MEXP=19937 -DHAVE_CONFIG_H -DPKGDATADIR='"/usr/local/share/pirs"' -I. -ISFMT-src-1.4

ifeq ($(shell uname -s),Darwin)
CFLAGS += -march=core2 -Bstatic -L/opt/local/lib -I/opt/local/include
else
CFLAGS += -static
endif

#GCC_VERSION := `gcc -dumpversion`
GCC_VERSION := $(shell gcc -dumpversion)
GCC_MAJOR := $(shell echo $(GCC_VERSION) | cut -d'.' -f1)
#GCC_MINOR := $(shell echo $(GCC_VERSION) | cut -d'.' -f2)
#GCC_SUB := $(shell echo $(GCC_VERSION) | cut -d'.' -f3)
GCC4_VER := $(shell echo $(GCC_VERSION) | cut -d'.' -f2,3)

CPPFLAGS += $(CFLAGS)

ifeq ($(GCC_MAJOR),4)
ifeq ($(shell expr 2.4 \<= $(GCC4_VER)),1)
CFLAGS += -march=native
CPPFLAGS += -march=native
endif
ifeq ($(shell expr 3.6 \<= $(GCC4_VER)),1)
#CPPFLAGS += -std=c++0x
endif
endif

LDLIB     = -lz -lpthread
LDFLAGS   = -Wl,-O1 -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both $(LDLIB)

SRCDIRS = . SFMT-src-1.4
SRCEXTS = .C .cc .cpp .CPP .c++ .cxx .c
OLevel= -O2
WARNFLAGS = -Wfloat-equal -Wall
SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))

# g++ -DHAVE_CONFIG_H -I. -DSFMT_MEXP=19937 -DPKGDATADIR='"/usr/local/share/pirs"' -Wall -Wextra -msse2  -fopenmp -g -O2 *.cpp SFMT-src-1.4/SFMT.c -lz

all: pirs
	@echo Souce Files:[$(SOURCES)]

pirs: $(SOURCES)
	g++ $(WARNFLAGS) $(CPPFLAGS) $(OLevel) -o $@ $^ $(LDFLAGS)

test: pirs
	cd test && sh pirs.sh
	ls -l test/

distclean:
	-rm test/Illumina_* test/ref22* test/*.out test/*.err test/EAMSS2_*

clean: distclean
	-rm pirs

.PHONY: distclean clean test
