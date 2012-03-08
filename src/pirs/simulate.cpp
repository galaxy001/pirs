#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "simulate.h"

//from ASCII of A(N) C G T to 0 1 2 3, auto dealing with upper or lower case.
//8bit char type, A=a=N=n=0, C=c=1, G=g=2, T=t=3, others as 4.
char alphabet[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};


//from ASCII of A C G T to 0 1 2 3, auto dealing with upper or lower case.
//8bit char type, A=a=0, C=c=1, G=g=2, T=t=3, others as 4.
char alphabet2[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//from 0 1 2 3 4 to A C G T N
char Bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//from 0 1 2 3 4 to T G C A N
char c_bases[5] ={
		'T', 'G', 'C', 'A', 'N'
};

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq)
{       
	int is_good = 1;
	for (int i = 0; i < seq.size(); i++)
	{   if (alphabet2[seq[i]] == 4)
			{   is_good = 0;
				break;
			}
	}
	return is_good;
}


//Rrealization of snp
char get_snp_match(char base, double snp_bias[]){
	char bases[4][3]={{'G','T','C'},  //A->G transition  A->C/T transvertion
			 {'C','G','A'},  	//T->C transition  T->G/A transvertion
			 {'T','A','G'},   //C->T transition  C->A/T transvertion
			 {'A','T','C'}};  //G->A transition  G->T/C transvertion
			 	
	int a;
	char n='N';
	switch (base)
	{
	case 'A': a=0;break;
	case 'T': a=1;break;
	case 'C': a=2;break;
	case 'G': a=3;break;
	default: return n;
	}
	
	double num = double(rand())/double(RAND_MAX);
	int i = 0;
	for(i=0; i<3; i++)
	{
		double p = snp_bias[i];
		if(num <= p){break;}
	}
	return bases[a][i];
	
}

//Produce heterozygous SNPs in multiploid
string Get_snp(string &seq,ofstream &snp,string id, double hetersnp_rate, double snp_bias[]){
	map<uint64_t,string> snp_lst;
	uint64_t seq_len=seq.size();
	uint64_t N = 1000000000000;//the max random value 
	uint64_t snp_num = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		if(seq[seq_index] == 'N'){continue;}
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= hetersnp_rate)
		{
			string ss = boost::lexical_cast <std::string>(seq[seq_index]) + "\t";
			//get snp base
			seq[seq_index]=get_snp_match(seq[seq_index], snp_bias);
			ss += boost::lexical_cast <std::string>(seq[seq_index]);
			//put into list file
			snp_lst[seq_index+1] = ss;
			snp_num++;
		}
	}

//	snp<<id<<" total SNP number: "<<snp_num<<endl;
	map<uint64_t, string>::const_iterator map_it = snp_lst.begin();
	while (map_it != snp_lst.end())
	{
		snp<<id<<"\t"<<map_it->first<<"\t"<<map_it->second<<endl;
		map_it++;
	}
	
	return seq;
}

string Get_invertion(string &seq, ofstream &invertion_file, string id, double SV_rate)
{
	uint64_t N = 1000000000000; //the max random value
	uint64_t seq_len=seq.size();
	double invertion_rate =  SV_rate/3;
//	cerr<< invertion_rate <<endl;
	int invertion_len[5] = {100, 200, 500, 1000, 2000};
	double array[5] = {0.70, 0.90, 0.97, 0.99, 1.0};
	map<uint64_t,int> invertion_lst;
	uint64_t invertion_num = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= invertion_rate)
		{
//			cerr << random_num<<endl;
			double random_num2 = (double)rand() / double(RAND_MAX);
//			cerr << random_num2<<endl;
			int j = 0;
			for(j=0; j<5; j++)
			{
				if(random_num2 <= array[j])
				{
					break;
				}
			}
			if(seq_index+invertion_len[j]>seq_len){continue;}
  		string sub_seq = seq.substr(seq_index,invertion_len[j]);
  		if(!check_seq(sub_seq)){continue;} //contain 'N' or other nonbases char
  		string rc_sub_seq = reversecomplementary(sub_seq);

			for(int i = 0; i < rc_sub_seq.size(); i++)
			{
				seq[seq_index+i] =  rc_sub_seq[i];
			}
			
  		invertion_lst[seq_index] = invertion_len[j];
  		invertion_num++;
		}
	}

//	invertion_file<<id<<" total invertion number: "<<invertion_num<<endl;
	map<uint64_t, int>::const_iterator map_it = invertion_lst.begin();
	while (map_it != invertion_lst.end())
	{
		invertion_file<<id<<"\t"<<map_it->first<<"\t"<<map_it->second<<endl;
		map_it++;
	}
	
	return seq;
}

//Getting insert sequence when heterozygous indel rate is bigger than 0
string get_insertion(int num){
	char base[]={'A','T','G','C'};
	string s;
	for (int a=0;a<num;a++)
	{
		int index=int(rand()%4);
		s.push_back(base[index]);
	}
	return s;
}

//Produce heterozygous indels in multiploid
string Get_indel(string &seq,ofstream &indel,string id1,double heterindel_rate,double SV_rate){
	uint64_t seq_len=seq.size();
	uint64_t N = 1000000000000; //the max random value 
	
	double small_insertion_rate = heterindel_rate/2;
	double samll_deletion_rate = heterindel_rate/2;
	//use the empirical distribution(1~6) from panda re-sequencing data
	int small_indel_len[6] = {1,2,3,4,5,6};
//	double array1[6] = {0.6482, 0.1717, 0.0720, 0.0729, 0.0218, 0.0134};
	double array1[6] = {0.6482, 0.8199, 0.8919, 0.9648, 0.9866, 1};
	
	double large_insertion_rate = SV_rate/3;
	double large_deletion_rate = SV_rate/3;
	//insertion,deletion and inversion share one third of total SV respectively
	//SV-length 100(70%),200(20%),500(7%),1000(2%),2000(1%)
	int large_indel_len[5] = {100, 200, 500, 1000, 2000};
//	double array2[5] = {0.70, 0.20, 0.07, 0.02, 0.01};
	double array2[5] = {0.70, 0.90, 0.97, 0.99, 1};
	
	map<uint64_t,string> indel_lst;
	//simulate small deletion
	uint64_t small_deletion_count = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= samll_deletion_rate)
		{
			double random_num2 = (double)rand() / double(RAND_MAX);
			int i = 0;
			for(i = 0; i < 6; i++)
			{
				if(random_num2 <= array1[i])
				{
					break;
				}
			}
			if(seq_index+small_indel_len[i]>seq_len){continue;}
			string sub_str = seq.substr(seq_index, small_indel_len[i]);
			if(!check_seq(sub_str)){continue;}
			small_deletion_count++;
			string ss = "-\t" + boost::lexical_cast <std::string>(small_indel_len[i]) + "\t";
			
			for (int k=0;k<small_indel_len[i];k++)
			{
				ss += seq[seq_index+k];
//					indel<<seq[seq_index+k];
				seq[seq_index+k]='D';
			}
			indel_lst[seq_index] = ss;
		}
	}
	
	//simulate large deletion
	uint64_t large_deletion_count = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= large_deletion_rate)
		{
			double random_num2 = (double)rand() / double(RAND_MAX);
			int i = 0;
			for(i = 0; i < 5; i++)
			{
				if(random_num2 <= array2[i])
				{
					break;
				}
			}
			if(seq_index+large_indel_len[i]>seq_len){continue;}
			string sub_str = seq.substr(seq_index, large_indel_len[i]);
			if(!check_seq(sub_str)){continue;}
			large_deletion_count++;
			string ss = "-\t" + boost::lexical_cast <std::string>(large_indel_len[i]) + "\t";
			
			for (int k=0;k<large_indel_len[i];k++)
			{
				ss += seq[seq_index+k];
//					indel<<seq[seq_index+k];
				seq[seq_index+k]='D';
			}
			indel_lst[seq_index] = ss;
		}
	}
	
	vector<short> insert(seq_len);
	string s;
	//simulate small insertion
	uint64_t small_insertion_count = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= small_insertion_rate)
		{
			double random_num2 = (double)rand() / double(RAND_MAX);
			int i = 0;
			for(i = 0; i < 6; i++)
			{
				if(random_num2 <= array1[i]){break;}	
			}
			insert[seq_index] = small_indel_len[i];
			small_insertion_count++;
		}
	}
	//simulate large insertion
	uint64_t large_insertion_count = 0;
	for(uint64_t seq_index = 0; seq_index < seq_len; seq_index++)
	{
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num <= large_insertion_rate)
		{
			double random_num2 = (double)rand() / double(RAND_MAX);
			int i = 0;
			for(i = 0; i < 5; i++)
			{
				if(random_num2 <= array2[i]){break;}	
			}
			insert[seq_index] = large_indel_len[i];
			large_insertion_count++;
		}
	}
	//modify the seq
	uint64_t total_insertion_count = 0;
	for (uint64_t i=0;i<seq_len;i++)
	{
		if (insert[i]>=1)
		{
			total_insertion_count++;
			string temp;
			temp=get_insertion(insert[i]);
			s+=temp;
			string ss = "+\t" + boost::lexical_cast <std::string>(int(insert[i])) + "\t" + temp;
			
			//判断该位点是否有deletion
			if(indel_lst[i][0] == 0)
			{
				indel_lst[i] = ss;
			}else{
				indel_lst[i] += "\n" + id1 + "\t" + boost::lexical_cast <std::string>(i) + "\t" + ss;
			}
//			indel<<id1<<"\t"<<i+1<<"\t+\t"<<int(insert[i])<<"\t"<<temp<<endl;
			if (seq[i]!='D')
			{
				s+=seq[i];
			}
		}else{
			if (seq[i]!='D')
			{
				s+=seq[i];
			}
		}
	}

//	indel <<"#small deletion times: "<<small_deletion_count << "; large deletion times: "<<large_deletion_count<<"; small insertion times: "<<small_insertion_count<<"; large insertion times: "<<large_insertion_count<<endl;
	map<uint64_t, string>::const_iterator map_it = indel_lst.begin();
	while (map_it != indel_lst.end())
	{
		indel<<id1<<"\t"<<map_it->first<<"\t"<<map_it->second<<endl;
		map_it++;
	}
	
	return s;
}


//Binary search the random number location
int search_location(double *Arr, uint64_t ArrNum, double random_num){

	uint64_t left = 0;
	uint64_t right = ArrNum;
	uint64_t middle = (left+right)/2;
	
	if(random_num < Arr[0]){return 0;}
	if(random_num > Arr[ArrNum-1]){ return ArrNum;}
	
	//return the first location of bigger than 0
	if(random_num == 0){
		for(uint64_t i = 0; i < ArrNum; i++)
		{
			if(Arr[i] > 0){return i;}
		}
	}
	
	while(!(random_num > Arr[middle-1] && random_num <= Arr[middle]))
	{
		if (left == middle || right == middle){
			cerr <<"left == middle || right == middle"<<endl;
			return middle;
		}
		if(random_num > Arr[middle]){
			left = middle;
		}else if(random_num < Arr[middle-1]){
			right = middle;
		}else if(random_num == Arr[middle-1]){  //return the first region that upper boundary equal to random_num
			for(uint64_t i = middle-1; i>0; i--)
			{
				if(Arr[i-1] != Arr[i]){return i;}
			}
		}
		middle = (left+right)/2;
	}
	
	return middle;
}

//get the reverse and complement sequence
string reversecomplementary (string read)
{	
	string rc_read;
	for (int64_t i=read.size()-1; i>=0; i--)
	{	
		rc_read.push_back(c_bases[alphabet[read[i]]]);
	}
	return rc_read;
}

//simulate normal distribution insertsize by Box-muller method
int simulate_insertsize(int mean, int sd)
{
	int insertsize = 0;
	double x1,x2,radius_pow2,y1 = 0.0;
	
	do
	{
		x1 = 2.0*double(rand())/double(RAND_MAX)-1.0;
		x2 = 2.0*double(rand())/double(RAND_MAX)-1.0;
		radius_pow2 = x1*x1+x2*x2;
	}while(radius_pow2>=1.0 || radius_pow2 == 0);
	
	radius_pow2 = sqrt(-2.0*log(radius_pow2)/radius_pow2);
	y1 = x1*radius_pow2;
	insertsize = int(mean+y1*sd);
	
	return insertsize;
}


//simulate GC bias
int simulate_GC_bias(string insert_str, double *GC_bias_abundance){
	int is_ignore = 0;
	int GC_num = 0;
	int insert_size = insert_str.size();
	//get GC rate
	for (int i=0; i<insert_size; i++)
	{	
		if(insert_str[i] == 'G' || insert_str[i] == 'C'){GC_num++;}
	}
	
	double GC_rate = double(GC_num)/double(insert_size);

	//get relative abundance
	double bias_abund = GC_bias_abundance[int(GC_rate*100)];
	
	//decide whether ignore this insert string.
	double num = double(rand())/double(RAND_MAX);
	if(num > bias_abund){is_ignore = 1;}
		
	return is_ignore;
}

//Set indel rate
void set_rate(int read_len, double total_rate, int max_num, vector <double>& rate)
{
  rate.resize(100);
  for(int i=1; i<=100; i++){
//    rate[i-1]= gsl_cdf_binomial_Q(i, total_rate, read_len);  //ART bug?
		rate[i-1]= gsl_ran_binomial_pdf(i, total_rate, read_len);
//    cerr<<"rate: "<<rate[i-1]<<endl;
    
  }
}

//simulate reads indel 
int simulate_reads_indel(vector<double> del_rate, vector<double> ins_rate, map<int,char,less<int> > &indel, int read_len)
{
    int ins_len=0, del_len=0;
    //deletion
    for(int i=(int)del_rate.size()-1; i>=0; i--){
    	double num = double(rand())/double(RAND_MAX);
        if(del_rate[i]>=num){
            del_len=i+1;
            for(int j=i; j>=0;){
            	double num1 = double(rand())/double(RAND_MAX);
                int pos=(int) floor((read_len-1)*num1); //invalid deletion positions: 0 or read_len-1
                if(pos==0) continue;
                if(indel.count(pos)==0){
                    indel[pos]='-';
                    j--;
                }
            }
            break;
        }
    }

    for(int i=ins_rate.size()-1; i>=0; i--){
    	double num = double(rand())/double(RAND_MAX);
        if(ins_rate[i]>=num){
            ins_len=i+1;
            for(int j=i; j>=0;){
            	double num1 = double(rand())/double(RAND_MAX);
                int pos=(int) floor(num1*read_len);
                if(indel.count(pos)==0){
                	double num2 = double(rand())/double(RAND_MAX);
                    short base=(short)ceil(num2*4);
                    switch(base){
                                case 1:
                                    indel[pos]='A';   break;
                                case 2:
                                    indel[pos]='C';   break;
                                case 3:
                                    indel[pos]='G';   break;
                                case 4:
                                    indel[pos]='T';  
                    }
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len-del_len);
}

//Get indel read
string ref2read(string seq_ref, map<int,char,less<int> > indel)
{
	string read;
    if(indel.size()==0){
        read=seq_ref;
        return read;
    }
    read.clear();
    int k=0;
    for(size_t i=0; i<seq_ref.size();){
        //cout<<i<<"\t"<<k<<endl;
        if(indel.count(k)==0){
            read.push_back(seq_ref[i]); i++; k++; 
        }
        else if(indel[k]=='-'){
            i++;k++;
        }
        else{
            read.push_back(indel[k]); k++;
        }
    }
    while(indel.count(k)>0){
        read.push_back(indel[k]);
        k++;
    }
    return read;
}
