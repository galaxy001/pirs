#ifndef __SIMULATE_ILLUMINA_READS_H_
#define __SIMULATE_ILLUMINA_READS_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include "simulate.h"
#include "gzstream.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;

//set and check input and ouput file
void set_and_check_file(igzstream &infile, igzstream &infile2, ofstream &log, ofstream &insert_log, ofstream &error_log);

//preview error profile, get dimensions of error profile
void preview_error_profile (string exe_path, ofstream &log);

//load GC depth and translate to GC abundance for simulate GC bias
void load_GC_depth_profile (string exe_path, ofstream &log);

//load error profile and get the simulate matrix
void load_error_profile(string exe_path, ofstream &log);

//get genome sequence
void Get_genome(igzstream &inf,igzstream &inf2,ofstream &log1);

//contral the quantity of simulate reads
uint64_t contral_reads_quantity(string id_line,string id,string &sequ,string &sequ2,ofstream &log2,uint64_t read_genome);

//simulate fastq reads
uint64_t simulate_fq_reads(string &seq,uint64_t seqlen, uint64_t rd_pair, string id_seq, ofstream &log3, uint64_t reads_all);

//simulate fasta reads
uint64_t simulate_fa_reads(string &seq,uint64_t seqlen, uint64_t rd_pair, string id_seq, ofstream &log3, uint64_t reads_all);

#endif
