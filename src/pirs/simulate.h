#ifndef __SIMULATE_H_
#define __SIMULATE_H_

using namespace std;

extern char alphabet[128];
extern char alphabet2[128];
extern char Bases[5];
extern char c_bases[5];

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq);

//get base base on quality score
char get_base_by_Qscore(char ref_base, int qscore);

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

//simulate diploid seq
string simulate_diploid_seq(string id,string raw_seq, ofstream &snp_info, ofstream &indel_info, ofstream &inversion_info, double hetersnp_rate, double heterindel_rate, double SV_rate, double snp_bias[]);

//get the reverse and complement sequence
string reversecomplementary (string read);

//Get reads-indel error
void get_reads_indel(int read_len, map<int,string,less<int> > &indel1, map<int,string,less<int> > &indel2, int &r1_indel_len, int &r2_indel_len, int InDel_max_len, double** InDel_matrix, int* InDel_num);

//Get indel read
string ref2read(string seq_ref, map<int,string,less<int> > indel, bool* is_insertion_pos);

#endif
