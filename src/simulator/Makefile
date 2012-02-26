all: simulate_diploid_genome simulate_illumina_reads

simulate_diploid_genome: simulate_diploid_genome.cpp simulate.cpp gzstream.cpp
	g++ -std=c++0x -O3 -march=native -mtune=generic -static -o $@ $^ -lz

simulate_illumina_reads: simulate_illumina_reads.cpp simulate.cpp load_file.cpp  gzstream.cpp
	g++ -std=c++0x -O3 -march=native -mtune=generic -static -o $@ $^ -lz

clean:
	-rm simulate_diploid_genome simulate_illumina_reads
