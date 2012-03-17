#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <boost/algorithm/string.hpp>
#include "simulate.h"
#include "gzstream.h"

using namespace std;

string input;
double hetersnp_rate=0.001;
double heterindel_rate=0.0001;
double big_SV_rate=0.000001; // structural variation rate
double snp_transition_by_transvertion_rate = 2; //transition_number ：transvertion_rate_number = 2
int output_type = 1;
ofstream outfile;
ogzstream gz_outfile;
string output_prefix = "ref_sequence";
double snp_bias[3]; //label transition and transvertion ratio

//get raw genome sequence.
void Get_raw_genome(igzstream &inf, ofstream &snp, ofstream &indel, ofstream &invertion);

//add snp and indel in raw seqence, and output result sequence.
void simulate_snp_indel_seq(string id_line,string id,string &sequ, ofstream &snp,ofstream &indel, ofstream &invertion);


void SimDiploid_Usage(){
	cout<<"\nDescription:"<<endl;
	cout<<endl;
	cout<<"\tIt is a program for simulating heterozygosis in diploid, such as SNP, InDel and structural variation (insertion, deletion and inversion). ";
	cout<<"When simulating SNP, we can set the total SNP rate with option -s and set the value of transition divided by transvertion with option -a. ";
	cout<<"When simulating InDel: they go halves on the total rate respectively, 1~6bp bases as the InDel number ,and the ";
	cout<<"rate distribution of each type as below: 1bp-64.82%, 2bp-17.17%, 3bp-7.20%, 4bp-7.29%, 5bp-2.18%, 6bp-1.34%, which is the empirical distribution from panda re-sequencing data, we can ";
	cout<<"set the total InDel rate with option -d. When simulating SV(structural variation): insertion, deletion and inversion, each of them share 1/3 of total rate, here we use this rate ";
	cout<<"distribution for simplicity: 100bp-70%, 200bp-20%, 500bp-7%, 1000bp-2%, 2000bp-1%, we can set the total SV rate with option -v. ";
	cout<<"Input sequence must be set ,because there is no default value."<<endl;
	cout<<endl<<"Program: pirs diploid"<<endl;
	cout<<endl<<"Usage:\t./pirs diploid [options]"<<endl;
	cout<<"\t-i	<string>	input reference genome sequence *.fa/*.fa.gz"<<endl;
	cout<<"\t-s	<double>	the heterozygous SNP rate of the diploid genome,default:"<<hetersnp_rate<<endl;
	cout<<"\t-a	<double>	the value of transition divided by transvertion for heterSNP,default:"<<snp_transition_by_transvertion_rate<<endl;
	cout<<"\t-d	<double>	the InDel rate of the diploid genome,default:"<<heterindel_rate<<endl;
	cout<<"\t-v	<double>	the structural variation rate(large insertion,deletion,invertion) of the diploid genome,default:"<<big_SV_rate<<endl;
	cout<<"\t-c	<int>   	ouput file type, 0:text, 1:*.gz, default:"<< output_type <<endl;
	cout<<"\t-o	<string>	output file prefix default:"<<output_prefix<<endl;
	cout<<"\t-h	        	output help infomation"<<endl;
	cout<<endl<<"Example:"<<endl;
	cout<<"\t1. ./pirs diploid -i ref_sequence.fa -s 0.001 -d 0.0001 -v 0.000001 -o ref_sequence >SimDiploid.out 2>SimDiploid.err"<<endl;
	exit(-1);
}

void SimDiploid_Getopt(int argc,char *argv[]){
	int c;
	while ((c=getopt(argc,argv,"i:s:d:o:a:v:c:h"))!=-1)
	{
		switch(c){
			case 'i': input=optarg;break;
			case 's': hetersnp_rate=strtod(optarg,NULL);break;
			case 'd': heterindel_rate=strtod(optarg,NULL);break;
			case 'a': snp_transition_by_transvertion_rate=strtod(optarg,NULL);break;
			case 'v': big_SV_rate=strtod(optarg,NULL);break;
			case 'c': output_type=atoi(optarg);break;
			case 'o': output_prefix=optarg;break;
			case 'h': SimDiploid_Usage();break;
			default: SimDiploid_Usage();
		}
	}
}

int simulate_diploid_genome(int argc, char *argv[])
{
	time_t time_start, time_end;
	time_start = time(NULL);
	srand((unsigned)time(NULL));
	
	if (argc==1)
	{
		SimDiploid_Usage();
	}
	SimDiploid_Getopt(argc,argv);
	
	argv--;
	
	if(hetersnp_rate < 0 || hetersnp_rate > 1){cerr<<"Error: hetersnp_rate should be set in 0~1 "<<endl;exit(-1);}
	if(snp_transition_by_transvertion_rate < 0){cerr<<"Error: the rate of transition/transvertion should be set greater than 0 "<<endl;exit(-1);}
	if(heterindel_rate < 0 || heterindel_rate > 1){cerr<<"Error: heterindel_rate should be set in 0~1 "<<endl;exit(-1);}
	if(big_SV_rate < 0 || big_SV_rate > 1){cerr<<"Error: large SV rate should be set in 0~1 "<<endl;exit(-1);}
	if(hetersnp_rate == 0 && heterindel_rate == 0 && big_SV_rate == 0){cerr<<"No variation for reference, please input at least one of parameter greater than 0!"<<endl; exit(-1);}
	if(output_type != 0 && output_type != 1){cerr <<"Error: -c (output_type) should be set 0 or 1!"<<endl; exit(-1);}
	
	
	//get snp bias
	if(hetersnp_rate > 0){
		double snp_sum = snp_transition_by_transvertion_rate + 1.0;
		double transition_rate =  double(snp_transition_by_transvertion_rate)/snp_sum;
		double transvertion_rate = 1.0/snp_sum;
		snp_bias[0] = transition_rate;
		snp_bias[1] = transition_rate+transvertion_rate/2.0;
		snp_bias[2] = 1;
	}
	
	//set output file name
	igzstream infile;
	infile.open(input.c_str());
	if (!infile)
	{
		cerr<<"Error:unable to open input file:"<<input<<endl;
		exit(-1);
	}
	string snp_output = output_prefix+"_snp.lst";
	string indel_output = output_prefix+"_indel.lst";
	string invertion_output = output_prefix+"_invertion.lst";
	if(hetersnp_rate != 0){output_prefix = output_prefix+".snp";}
	if(heterindel_rate != 0){output_prefix = output_prefix+".indel";}
	if(big_SV_rate !=0){output_prefix = output_prefix+".invertion";}
		
	string output_ref_file;
	if(!output_type){
		output_ref_file = output_prefix + ".fa";
	}else{
		output_ref_file = output_prefix + ".fa.gz";
	}
	
	if(!output_type){
		outfile.open(output_ref_file.c_str());
  	if (!outfile)
  	{
  		cerr<<"Error:unable to open output file:"<<output_ref_file<<endl;
  		exit(-1);
  	}
	}else{
		gz_outfile.open(output_ref_file.c_str());
  	if (!gz_outfile)
  	{
  		cerr<<"Error:unable to open output file:"<<output_ref_file<<endl;
  		exit(-1);
  	}
	}

	ofstream SNP_File,Indel_File,Ivertion_File;

	if(hetersnp_rate>0){
  	SNP_File.open(snp_output.c_str());
  	if(!SNP_File)
  	{
  		cerr<<"Error:unable to open output file:"<<snp_output<<endl;
  		exit(-1);
  	}
	}
	
	if(heterindel_rate>0){
  	Indel_File.open(indel_output.c_str());
  	if(!Indel_File)
  	{
  		cerr<<"Error:unable to open output file:"<<indel_output<<endl;
  		exit(-1);
  	}
	}
	
	if(big_SV_rate>0){
  	Ivertion_File.open(invertion_output.c_str());
  	if(!Ivertion_File)
  	{
  		cerr<<"Error:unable to open output file:"<<invertion_output<<endl;
  		exit(-1);
  	}
	}
	
	Get_raw_genome(infile,SNP_File,Indel_File,Ivertion_File);
	
	infile.close();
	
	if(!output_type){
		outfile.close();
	}else{
		gz_outfile.close();
	}
	
	time_end = time(NULL);
	cerr<<"All done! Run time: "<<time_end-time_start<<"s."<<endl;
	
	return 0;
}

//get raw genome sequence.
void Get_raw_genome(igzstream &inf, ofstream &snp_file, ofstream &indel_file, ofstream &invertion_file )
{
	string line,id,id_line,seq;
	while (getline(inf,line,'\n'))
	{
		if (line[0]=='>')
		{
			if (seq!="")
			{	
				cerr<<"Have finished reading scaffold "<<id<<endl;
				//start to simulate one scaffold
				simulate_snp_indel_seq(id_line,id,seq,snp_file,indel_file,invertion_file);
				seq="";
			}
			id_line = line;
			line.erase(0,1);
//			id=line;
			int pos=line.find(" ");
			line=line.substr(0,pos);
			id=line;
		}else{
			seq+=line;
		}		
	}
	cerr<<"Have finished reading scaffold "<<id<<endl;
	//start to simulate one scaffold
	simulate_snp_indel_seq(id_line,id,seq,snp_file,indel_file, invertion_file);
}

//add snp and indel in raw seqence, and output result sequence.
void simulate_snp_indel_seq(string id_line,string id,string &sequence, ofstream &snp_file,ofstream &indel_file, ofstream &invertion_file)
{
	//convert lower case to upper case 
	boost::to_upper(sequence);

	if (hetersnp_rate>0 || heterindel_rate>0 || big_SV_rate>0) //heterozygous SNP and heterozygous indel exists in diploid
	{
		if (hetersnp_rate>0)
		{
			cerr<<"Begin to simulate snp"<<endl;
			sequence=Get_snp(sequence,snp_file,id,hetersnp_rate,snp_bias);
			cerr<<"Have finished simulating snp"<<endl;
		}
		if (big_SV_rate>0)
		{
			cerr<<"Begin to simulate invertion"<<endl;
			sequence=Get_invertion(sequence,invertion_file,id,big_SV_rate);
		}
		if (heterindel_rate>0 || big_SV_rate>0)
		{
			cerr<<"Begin to simulate indel"<<endl;
			sequence=Get_genome_indel(sequence,indel_file,id,heterindel_rate,big_SV_rate);
			cerr<<"Have finished simulating indel"<<endl;
		}

		if(!output_type){
			outfile<<id_line<<endl;
			uint64_t seq_len = sequence.size();
			
			//50bp base per line
			uint64_t counter = 0; 
			while(seq_len > 50)
			{
				outfile<<sequence.substr(counter*50,50)<<endl;
				seq_len -= 50;
				counter++;
			}
			outfile<<sequence.substr(counter*50)<<endl;
		}else{
			gz_outfile<<id_line<<endl;
			uint64_t seq_len = sequence.size();

			uint64_t counter = 0;
			while(seq_len > 50)
			{
				gz_outfile<<sequence.substr(counter*50,50)<<endl;
				seq_len -= 50;
				counter++;
			}
			gz_outfile<<sequence.substr(counter*50)<<endl;
		}

		cerr<<"Have finished simulating "<<id_line<<endl;
	}
}


