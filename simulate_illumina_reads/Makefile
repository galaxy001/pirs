all: simulate_snp_indel_seq simulate_illumina_reads

simulate_snp_indel_seq: simulate_snp_indel_seq.cpp simulate.cpp gzstream.cpp
	g++ -std=c++0x -O3 -static -o $@ $^ -lz

simulate_illumina_reads: simulate_illumina_reads.cpp simulate.cpp gzstream.cpp
	g++ -std=c++0x -O3 -static -o $@ $^ -lz

clean:
	rm simulate_snp_indel_seq simulate_illumina_reads
