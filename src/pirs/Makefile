ifeq ($(shell uname -s),Darwin)
CPPFLAGS += -march=core2 -Bstatic -L/opt/local/lib -I/opt/local/include
else
CPPFLAGS += -std=c++0x -march=native -static
endif

all: pirs

pirs: main.cpp simulate_diploid_genome.cpp simulate.cpp gzstream.cpp simulate_Illumina_reads.cpp load_file.cpp MaskQvalsByEamss.cpp
	g++ $(CPPFLAGS) -mtune=generic -O3 -o $@ $^ -lz

test: pirs
	cd test && sh pirs.sh
	ls -l test/

clean:
	-rm test/Illumina_* test/ref22* test/*.out test/*.err test/EAMSS2_*

allclean: clean
	-rm pirs
