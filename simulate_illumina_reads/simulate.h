#ifndef __SIMULATE_H_
#define __SIMULATE_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;

extern char alphabet[128];
extern char alphabet2[128];
extern char Bases[5];
extern char c_bases[5];

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq);

//find random number location
int search_location(double *Arr, uint64_t ArrNum, double random_num);

//simulate normal distribution insertsize by Box-muller method
int simulate_insertsize(int mean, int sd);

//simulate GC bias
int simulate_GC_bias(string insert_str, double *GC_bias_abundance);

//Realization of snp 
char get_snp_match(char base, double snp_bias[]);

//Getting insert sequence when heterozygous indel rate is bigger than 0
string get_insertion(int num);

//get the reverse and complement sequence
string reversecomplementary (string read);

//Produce heterozygous SNPs
string Get_snp(string &seq,ofstream &snp,string id, double hetersnp_rate, double snp_bias[]);

//Getting invert sequence
string Get_invertion(string &seq, ofstream &invertion_file, string id, double SV_rate);

//Getting indel sequence
string Get_indel(string &seq,ofstream &indel,string id1,double heterindel_rate,double SV_rate);

#endif
