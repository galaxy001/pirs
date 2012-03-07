all: pirs

pirs: main.cpp simulate_diploid_genome.cpp simulate.cpp gzstream.cpp simulate_Illumina_reads.cpp load_file.cpp
	g++ -std=c++0x -march=native -mtune=generic -O3 -static -o $@ $^ -lz -lgsl

clean:
	-rm pirs
