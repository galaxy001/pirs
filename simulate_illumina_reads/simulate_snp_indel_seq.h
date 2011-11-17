#ifndef __SIMULATE_SNP_INDEL_SEQ_H_
#define __SIMULATE_SNP_INDEL_SEQ_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include "simulate.h"
#include "gzstream.h"

using namespace std;

//get raw genome sequence.
void Get_raw_genome(igzstream &inf, ofstream &snp, ofstream &indel, ofstream &invertion);

//add snp and indel in raw seqence, and output result sequence.
void simulate_snp_indel_seq(string id_line,string id,string &sequ, ofstream &snp,ofstream &indel, ofstream &invertion);

#endif

