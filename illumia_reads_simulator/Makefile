all: simulate_diploid_genome simulate_Illumina_reads

simulate_diploid_genome: simulate_diploid_genome.cpp simulate.cpp gzstream.cpp
	g++ -std=c++0x -O3 -static -o $@ $^ -lz

simulate_Illumina_reads: simulate_Illumina_reads.cpp simulate.cpp load_file.cpp  gzstream.cpp
	g++ -std=c++0x -O3 -static -o $@ $^ -lz

clean:
	rm simulate_diploid_genome simulate_Illumina_reads
