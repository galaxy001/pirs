#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdint.h>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
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

//get snp base
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

//string get_diploid_seq(string raw_seq, string id, ofstream &indel_info, ofstream &snp_info, ofstream &inversion_info, double hetersnp_rate, double heterindel_rate, double SV_rate, double snp_bias[])
string simulate_diploid_seq(string id,string raw_seq, ofstream &snp_info, ofstream &indel_info, ofstream &inversion_info, 
	double hetersnp_rate, double heterindel_rate, double SV_rate, double snp_bias[])
{
	//convert lower case to upper case 
	boost::to_upper(raw_seq);
	
	cerr<< "Begin to simulate diploid sequence..."<<endl;
	
	//probability distribution
	//small insertion and deletion share half of total heterindel_rate
	double small_insertion_rate = heterindel_rate/2.0;
	double samll_deletion_rate = heterindel_rate/2.0;
	//insertion,deletion and inversion share one third of total SV respectively
	double sv_insertion_rate = SV_rate/3.0;
	double sv_deletion_rate = SV_rate/3.0;
	double sv_inversion_rate =  SV_rate/3.0;	
	
	//indel-length 1(64.82%), 2(17.17%), 3(7.2%), 4(7.29%), 5(2.18%), 6(1.34%)
	int indel_len[6] = {1,2,3,4,5,6};
	double array1[6] = {0.6482, 0.8199, 0.8919, 0.9648, 0.9866, 1};
	
	//SV-length 100(70%),200(20%),500(7%),1000(2%),2000(1%)
	int sv_len[5] = {100, 200, 500, 1000, 2000};
	double array2[5] = {0.70, 0.90, 0.97, 0.99, 1};
	
	//stat-counter
	uint64_t snp_num = 0, small_insert_time = 0, small_insert_num = 0, small_delet_time = 0, 
		small_delet_num = 0, sv_insert_time = 0, sv_insert_num = 0, sv_delet_time = 0, sv_delet_num = 0,
		sv_invert_time = 0, sv_invert_num = 0;
	
	string diploid_seq;
	
	uint64_t seq_len = raw_seq.size();
	
	//start to simulate diploid sequence
	//simulation order: SNP -> small-deletion -> small-insertion -> SV-deltion -> SV-insertion -> SV-inversion
	for(uint64_t raw_seq_index = 0, diploid_seq_index = 0; raw_seq_index < seq_len; diploid_seq_index++, raw_seq_index++)
	{
		//////////////////simulate snp//////////////////
		double random_num = (double)rand() / double(RAND_MAX);
		if(random_num < hetersnp_rate && raw_seq[raw_seq_index] != 'N')
		{
			//get snp base
			diploid_seq += get_snp_match(raw_seq[raw_seq_index], snp_bias);
			
			snp_num++;
			uint64_t pos = raw_seq_index + 1;
			snp_info << id << "\t" << pos << "\t" << raw_seq[raw_seq_index] << "\t" << diploid_seq[diploid_seq_index] <<endl;
			
		}else{
			diploid_seq += raw_seq[raw_seq_index];
		}
	
		//////////////////simulate deletion//////////////////
		random_num = (double)rand() / double(RAND_MAX);
		if(random_num < samll_deletion_rate)
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
			if(raw_seq_index + indel_len[i] <= seq_len - 1)
			{
  			string sub_str = raw_seq.substr(raw_seq_index + 1, indel_len[i]);
  			if(check_seq(sub_str))  //check whether contain 'N' or other nonbases char
  			{ 
    			string del_seq;
    			for (int k = 1; k <= indel_len[i]; k++)
    			{
    				del_seq += raw_seq[raw_seq_index+k];
    			}
    			indel_info << id << "\t" << raw_seq_index + 1 << "\t-\t" << indel_len[i] << "\t" << del_seq <<endl;
    			raw_seq_index += indel_len[i];
    			small_delet_time++;
    			small_delet_num += indel_len[i];
  			}
			}
		}else{
			//////////////////simulate insertion//////////////////
			random_num = (double)rand() / double(RAND_MAX);
  		if(random_num < small_insertion_rate) 
  		{
  			double random_num2 = (double)rand() / double(RAND_MAX);
  			int i = 0;
  			for(i = 0; i < 6; i++)
  			{
  				if(random_num2 <= array1[i]){break;}	
  			}
  			string insert_seq = get_insertion(indel_len[i]);
  			
  			indel_info << id << "\t" << raw_seq_index + 1 << "\t+\t" << indel_len[i] << "\t" << insert_seq <<endl;
  			diploid_seq += insert_seq;
  			diploid_seq_index += indel_len[i];

  			small_insert_time++;
  			small_insert_num += indel_len[i];
  		}else{
  			//////////////////simulate sv-deletion//////////////////
  			random_num = (double)rand() / double(RAND_MAX);
        if(random_num < sv_deletion_rate) 
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
          if(raw_seq_index + sv_len[i] <= seq_len - 1 )
          {
            string sub_str = raw_seq.substr(raw_seq_index + 1, sv_len[i]);
            if(check_seq(sub_str)) //check whether contain 'N' or other nonbases char
            {
        			string del_seq;
        			for (int k = 1; k <= sv_len[i]; k++)
        			{
        				del_seq += raw_seq[raw_seq_index+k];
        			}
        			indel_info << id << "\t" << raw_seq_index + 1 << "\t-\t" << sv_len[i] << "\t" << del_seq <<endl;
        			raw_seq_index += sv_len[i];
        			sv_delet_time++;
        			sv_delet_num += sv_len[i];
          	}
        	}
        }else{
        	//////////////////simulate sv-insertion//////////////////
  				random_num = (double)rand() / double(RAND_MAX);
      		if(random_num < sv_insertion_rate)
      		{
      			double random_num2 = (double)rand() / double(RAND_MAX);
      			int i = 0;
      			for(i = 0; i < 5; i++)
      			{
      				if(random_num2 <= array2[i]){break;}	
      			}
      			string insert_seq = get_insertion(sv_len[i]);
      			
      			indel_info << id << "\t" << raw_seq_index + 1 << "\t+\t" << sv_len[i] << "\t" << insert_seq <<endl;
      			diploid_seq += insert_seq;
      			diploid_seq_index += sv_len[i];
      			sv_insert_time++;
      			sv_insert_num += sv_len[i];
      		}else{
      			//////////////////simulate sv-inversion//////////////////
        		random_num = (double)rand() / double(RAND_MAX);
        		if(random_num < sv_inversion_rate)
        		{
        			double random_num2 = (double)rand() / double(RAND_MAX);
        			int j = 0;
        			for(j=0; j<5; j++)
        			{
        				if(random_num2 <= array2[j])
        				{
        					break;
        				}
        			}
        			if(raw_seq_index+sv_len[j] <= seq_len - 1)
        			{
            		string sub_seq = raw_seq.substr(raw_seq_index + 1, sv_len[j]);
            		if(check_seq(sub_seq)) //check whether contain 'N' or other nonbases char
            		{ 
            			inversion_info << id << "\t" << raw_seq_index + 1 << "\t" << sv_len[j] <<endl;
              		string rc_sub_seq = reversecomplementary(sub_seq);

            			diploid_seq += rc_sub_seq;
            			diploid_seq_index += sv_len[j];
            			raw_seq_index += sv_len[j];
        
              		sv_invert_time++;
              		sv_invert_num += sv_len[j];
            		}
          		}
        		}
      		}
        }
  		}
		}
	}
	
	cout<<"\n*************heterozygosis table for seq: "<<id<<"*************"<<endl
		<<"raw_seq_length\tnew_seq_length"<<endl
		<<uint64_t(raw_seq.size())<<"\t"<<uint64_t(diploid_seq.size())<<endl
		<<"heterSNP & heterINDEL:"<<endl
		<<"total_snp_num\tdeletion_time\ttotal_deletion_length\tinsertion_time\ttotal_insertion_length"<<endl
		<<snp_num<<"\t"<<small_delet_time<<"\t"<<small_delet_num<<"\t"<<small_insert_time<<"\t"<<small_insert_num<<endl
		<<"structural variation:"<<endl
		<<"deletion_time\ttotal_deletion_length\tinsertion_time\ttotal_insertion_length\tinversion_time\ttotal_inversion_length"<<endl
		<<sv_delet_time<<"\t"<<sv_delet_num<<"\t"<<sv_insert_time<<"\t"<<sv_insert_num<<"\t"<<sv_invert_time<<"\t"<<sv_invert_num<<endl;

	cerr<<"Have finished simulating seq: "<<id<<endl;
	
	return diploid_seq;
	
}

char get_base_by_Qscore(char ref_base, int qscore)
{
	double error_rate = pow(10.0, -qscore/10.0);
	double random_num = (double)rand() / double(RAND_MAX);
	char base[]={'A','T','G','C'};
	if(random_num < error_rate)
	{
		while(1)
		{
  		int index=int(rand()%4);
  		if(base[index] != ref_base)
  		{
  			return base[index];
  		}
  	}
	}else{
		return ref_base;
	}
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

void get_reads_indel(int read_len, map<int,string,less<int> > &indel1, map<int,string,less<int> > &indel2, int &r1_indel_len, int &r2_indel_len, int InDel_max_len, double** InDel_matrix, int *InDel_num)
{
	//read1 indel
	int del_num = 0;
	int ins_num = 0;
	for(int i = 0; i < read_len; i++)
	{
  	double num=double(rand())/double(RAND_MAX);
		int location = search_location(InDel_matrix[i], InDel_max_len*2+1, num);
		int indel_num = InDel_num[location]; 
		if(indel_num == 0) //no indel
		{
		}else if(indel_num < 0) //deletion
		{
			string tem;
			for(int l = indel_num; l < 0; l++)
			{
				tem+="-";del_num++;
			}
			indel1[i] = tem;
		}else{//insertion
			string tem;
			for(int k = 0; k < indel_num; k++)
			{
				ins_num++;
				double num2 = double(rand())/double(RAND_MAX);
        short base=(short)ceil(num2*4);
        switch(base){
          case 1:tem+="A";break;
          case 2:tem+="C";break;
          case 3:tem+="G";break;
          case 4:tem+="T";
        }
			}
			indel1[i] = tem;
			i+=indel_num;
		}
	}
	
	r1_indel_len = ins_num - del_num;
	
	//read2 indel
	ins_num = 0;
	del_num = 0;
	for(int i = 0; i < read_len; i++)
	{
  	double num=double(rand())/double(RAND_MAX);
		int location = search_location(InDel_matrix[i+read_len], InDel_max_len*2+1, num);
		int indel_num = InDel_num[location]; 
		if(indel_num == 0) //no indel
		{
		}else if(indel_num < 0) //deletion
		{
			string tem;
			for(int l = indel_num; l < 0; l++)
			{
				tem+="-";del_num++;
			}
			indel2[i] = tem;
		}else{//insertion
			string tem;
			for(int k = 0; k < indel_num; k++)
			{
				ins_num++;
				double num2 = double(rand())/double(RAND_MAX);
        short base=(short)ceil(num2*4);
        switch(base){
          case 1:tem+="A";break;
          case 2:tem+="C";break;
          case 3:tem+="G";break;
          case 4:tem+="T";
        }
			}
			indel2[i] = tem;
			i+=indel_num;
		}
	}
	r2_indel_len = 	ins_num - del_num;
}

//Get indel read
string ref2read(string seq_ref, map<int,string,less<int> > indel, bool* is_insertion_pos)
{
	string read;
  if(indel.size()==0){
    read=seq_ref;
    return read;
  }
  read.clear();
  int k=0;
  read.push_back(seq_ref[0]);
  for(int i=1; i<seq_ref.size();){
    //cout<<i<<"\t"<<k<<endl;
    if(indel.count(k)==0){
        read.push_back(seq_ref[i]); i++; k++; 
    }
    else if(indel[k][0]=='-'){
    	i+=indel[k].size();
    	k++;
    }
    else{
      read+=indel[k];
      for(int l =0 ; l < indel[k].size(); l++)
      {
      	is_insertion_pos[k+l+1] = 1;
      } 
      k+=indel[k].size();
    }
  }
  while(indel.count(k)>0){
    read+=indel[k];
    for(int l =0 ; l < indel[k].size(); l++)
    {
    	is_insertion_pos[k+l+1] = 1;
    } 
    k+=indel[k].size();
  }
  return read;
}
