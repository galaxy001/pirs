#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdint.h>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"
#include "global.h"
#include "load_file.h"

//from ASCII of A C G T to 0 1 2 3, auto dealing with upper or lower case.
//8bit char type, A=a=0, C=c=1, G=g=2, T=t=3, others as 4.
char alphabet3[128] =
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

double qscore2error(int qscore)
{
	double error_rate = pow(10.0, -qscore / 10.0);
	return error_rate;
}

void Qscore_cal(double time, int *q2q, int quality_num)
{
	for(int i = 0; i < quality_num; i++)
	{
		int new_q = int(-10 * log10(qscore2error(i) * time) + 0.5);
		if(new_q <= 2)
		{
			q2q[i] = 2;
		}
		else if(new_q >= quality_num - 1)
		{
			q2q[i] = quality_num - 1;
		}
		else
		{
			q2q[i] = new_q;
		}
	}
}

//check and open the outfile file
void set_and_check_file(PARAMETER InputParameter, igzstream &infile, igzstream &infile2, ofstream &outfile1,	ofstream &outfile2, 
	ogzstream &gz_outfile1, ogzstream &gz_outfile2, ofstream &insert_log, ofstream &error_log, ogzstream &infor_outfile)
{
	//set and open output file
	infile.open(InputParameter.Input_ref1.c_str());
	if (!infile)
	{
		cerr<<"Error:unable to open input file:"<<InputParameter.Input_ref1<<endl;
		exit(-1);
	}
	
	if(InputParameter.Input_ref2!="")
	{
		infile2.open(InputParameter.Input_ref2.c_str());
  	if (!infile2)
  	{
  		cerr<<"Error:unable to open input file:"<<InputParameter.Input_ref2<<endl;
  		exit(-1);
  	}
	}
	
	string output_file1, output_file2, output_insert_distr, output_error_distr, output_infor;
	string infix_name = "_"+boost::lexical_cast <std::string>(InputParameter.Read_length)+"_"+boost::lexical_cast <std::string>(InputParameter.Insertsize_mean);
	//read file1
	if(InputParameter.Is_simulate_quality){
		output_file1 = InputParameter.Output_prefix+infix_name+"_1.fq";
		//read file2
		output_file2 = InputParameter.Output_prefix+infix_name+"_2.fq";
	}else{
		output_file1 = InputParameter.Output_prefix+infix_name+"_1.fa";
		//read file2
		output_file2 = InputParameter.Output_prefix+infix_name+"_2.fa";
	}
	//compress ouput
	if(InputParameter.Output_type == 1){
		output_file1+=".gz";
		output_file2+=".gz";
	}

	//insert size distribution file
	output_insert_distr = InputParameter.Output_prefix+infix_name+".insertsize.distr";
	//error rate distribution file
	output_error_distr = InputParameter.Output_prefix+infix_name+".error_rate.distr";
	//read information file
	output_infor = InputParameter.Output_prefix+infix_name+".read.info.gz";
	//open file
  if(!InputParameter.Output_type){
  	outfile1.open(output_file1.c_str());
  	outfile2.open(output_file2.c_str());
  }else{
  	gz_outfile1.open(output_file1.c_str());
  	gz_outfile2.open(output_file2.c_str());
  }

  insert_log.open(output_insert_distr.c_str());
  error_log.open(output_error_distr.c_str());
  infor_outfile.open(output_infor.c_str());
	
	//check file
	if(!InputParameter.Output_type)
	{
		if (!outfile1 || !outfile2)
		{
			cerr<<"Error:unable to open output file."<<endl;
			exit(1);
		}
	}else{
		if (!gz_outfile1 || !gz_outfile2)
		{
			cerr<<"Error:unable to open output file."<<endl;
			exit(1);
		}
	}

	if(!insert_log){cerr<<"Error:unalbe to create insertsize distribution file."<<endl;exit(1);}
	if(!error_log){cerr<<"Error:unalbe to create error rate distribution file."<<endl;exit(1);}
	if(!infor_outfile){cerr<<"Error:unalbe to create reads information file."<<endl;exit(1);}
}

//get the attribute of Base-calling  profile
void preview_BaseCalling_profile (PARAMETER InputParameter, string exe_path, int &ref_Base_num, 
	int &statistical_Cycle_num, int &seq_Base_num, int &quality_num, double &statistical_average_error_rate)
{
	string matrix_file;
	
	if(InputParameter.BaseCalling_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + BASE_CALLING_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.BaseCalling_profile;
	}
	
  igzstream infile;
  infile.open(matrix_file.c_str());
  if ( ! infile )
	{
		cerr << "fail to open input file" << matrix_file << endl;
	}
	uint64_t total_count_sum = 0;
	uint64_t total_error_sum = 0;
	uint64_t total_correct_sum = 0;

	string lineStr;
	bool isEnd = 0;
	while (getline( infile, lineStr, '\n' ))
	{ 
		if(isEnd){break;}
		if (lineStr[0] == '#')
		{ if (lineStr[1] == 'D' && lineStr[2] == 'i' && lineStr[3] == 'm')
			{	vector<string> lineVec;
				boost::split(lineVec,lineStr, boost::is_any_of(":,; \t\n"), boost::token_compress_on);
				ref_Base_num = atoi(lineVec[2].c_str());
				statistical_Cycle_num = atoi(lineVec[4].c_str());
				
				if(InputParameter.Read_length > statistical_Cycle_num/2){cerr<<"Error: according to the Base-calling profile, program can be simulate "<< statistical_Cycle_num/2 <<"bp read length at most, please set read length again!"<<endl;exit(-1);}
				seq_Base_num = atoi(lineVec[6].c_str());
				quality_num = atoi(lineVec[8].c_str());
			}
			else
			{
				continue;
			}
		}
		else
		{ 
			if(lineStr == "[DistMatrix]"){
				while(getline( infile, lineStr, '\n' ))
				{
					if(lineStr == "" || lineStr[0] == '#'){continue;}
					if(lineStr == "<<END"){isEnd = 1; break;}
					vector<string> lineVec;
    			boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
    			
    			char ref_base = lineVec[0][0];
    			int current_cycle = atoi(lineVec[1].c_str());
    				
    			if(current_cycle <= InputParameter.Read_length || (current_cycle >= statistical_Cycle_num/2 && current_cycle <= statistical_Cycle_num/2+InputParameter.Read_length)){
        		for(int i = 0; i < seq_Base_num; i++)
        		{
        			for(int j = 0; j < quality_num; j++)
        			{
        				string num = lineVec[i * quality_num + j + 2];
        				uint64_t freqnum = boost::lexical_cast<uint64_t>(num);
    						
        				total_count_sum+=freqnum;
        				if( i == alphabet3[ref_base])  // ACGT  
        				{
        					total_correct_sum+=freqnum;
        				}else{
        					total_error_sum+=freqnum;
        				}
        			}
        		}
      		}
				}
			}else{
				continue;
			}
		}
	}
	
	infile.close();
	
	//get BaseCalling profile average error rate
	statistical_average_error_rate = double(total_error_sum)/double(total_count_sum);
  
}

//transform quality score according error-rate setting by user.
void transform_quality_by_error_rate(PARAMETER InputParameter, string exe_path, int *Qval2Qval, int seq_Base_num, int quality_num, int statistical_Cycle_num)
{
	if(InputParameter.Error_rate <= pow(10.0, -(quality_num - 1)/10.0))
	{
		for(int i = 0; i < quality_num; i++)
		{
			Qval2Qval[i] = quality_num - 1; //max quality score = quality_num - 1;
		}
		return;
	}
	
	string matrix_file;
	uint64_t *Qscor_distr = new uint64_t[quality_num]; //number distribution of Quality score 
	double *Qscore_ratio_distr = new double[quality_num];
	for(int i = 0; i < quality_num; i++)
	{
		Qscor_distr[i] = 0;
		Qscore_ratio_distr[i] = 0;
	}
	
	if(InputParameter.BaseCalling_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + BASE_CALLING_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.BaseCalling_profile;
	}
	
  igzstream infile;
  infile.open(matrix_file.c_str());
  if ( ! infile )
	{
		cerr << "fail to open input file" << matrix_file << endl;
	}
	uint64_t total_count_sum = 0;
	uint64_t total_error_sum = 0;
	uint64_t total_correct_sum = 0;

	string lineStr;
	bool isEnd = 0;
	while (getline( infile, lineStr, '\n' ))
	{ 
		if(isEnd){break;}
		if(lineStr == "[DistMatrix]"){
			while(getline( infile, lineStr, '\n' ))
			{
				if(lineStr == "" || lineStr[0] == '#'){continue;}
				if(lineStr == "<<END"){isEnd = 1; break;}
				vector<string> lineVec;
  			boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
  			
  			char ref_base = lineVec[0][0];
  			int current_cycle = atoi(lineVec[1].c_str());
  				
  			if(current_cycle <= InputParameter.Read_length || (current_cycle >= statistical_Cycle_num/2 && current_cycle <= statistical_Cycle_num/2+InputParameter.Read_length)){
      		for(int i = 0; i < seq_Base_num; i++)
      		{
      			for(int j = 0; j < quality_num; j++)
      			{
      				string num = lineVec[i * quality_num + j + 2];
      				uint64_t freqnum = boost::lexical_cast<uint64_t>(num);
  						Qscor_distr[j] += freqnum;
  						
      				total_count_sum += freqnum;
      				if( i == alphabet3[ref_base])  // ACGT  
      				{
      					total_correct_sum+=freqnum;
      				}else{
      					total_error_sum+=freqnum;
      				}
      			}
      		}
    		}
			}
		}else{
			continue;
		}
		
	}
	infile.close();
	
	//ratio distribution of Quality score 
	for(int i = 0; i < quality_num; i++)
	{
		Qscore_ratio_distr[i] = double(Qscor_distr[i]) / double(total_count_sum);
		//cout<<i<<": "<<Qscore_ratio_distr[i]<<endl;
	}
	
	delete[] Qscor_distr;
	
	//transform to new quality
  double min = 0.00015;
  double max = 6000;
  double mid = 0;
  double n = 0;
  double error_tmp;

	while(1)
	{
		mid = (min + max) / 2;
		error_tmp = 0;
		Qscore_cal(mid, Qval2Qval, quality_num);
		
		for(int k = 0; k < quality_num; k++)
		{
			error_tmp += Qscore_ratio_distr[k] * qscore2error(Qval2Qval[k]);
		}
		double len = 0;
		if(error_tmp > InputParameter.Error_rate)
		{
			len = max - mid;
			max = mid;
		}
		else if(error_tmp < InputParameter.Error_rate)
		{
			len = mid - min;
			min = mid;
		}
		else
		{
			break;
		}
		
		if(fabs(InputParameter.Error_rate - error_tmp) < 0.0001 || len < 0.00001)
		{
			break;
		}
	}
	n = mid;
	
	cout<<"\nShift quality score for simulating substiution-error:"<<endl
		<<"\traw_qulaity_score\ttransferred_quality_score"<<endl;
	for(int i = 0; i < quality_num ; i++)
	{
		cout<<"\t"<<i<<"\t"<<Qval2Qval[i]<<endl;
	}
	
	delete[] Qscore_ratio_distr;
}

void preview_InDel_error_profile (PARAMETER InputParameter, string exe_path, int &Statistical_Cycle_num2, int &InDel_max_len, uint64_t &read1_count, uint64_t &read2_count)
{
	string matrix_file;
	
	if(InputParameter.InDel_error_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + INDEL_ERROR_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.InDel_error_profile;
	}
	
  igzstream infile;
  infile.open(matrix_file.c_str());
  if ( ! infile )
	{
		cerr << "fail to open input file" << matrix_file << endl;
	}
	
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{ 
		if(lineStr == ""){continue;}
		if(lineStr == "<<END"){break;}

		vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);				
		
		//"Read_Length = 100"
		if (lineVec[0] == "Read_Length")
		{
			Statistical_Cycle_num2 = atoi(lineVec[2].c_str()) * 2;
			if(InputParameter.Read_length > Statistical_Cycle_num2/2){cerr<<"Error: according to the InDel-error profile, program can be simulate "<< Statistical_Cycle_num2/2 <<"bp read length at most, please set read length again!"<<endl;exit(-1);}
		}
		
		//"Read_1_Count = 147091454"
		if (lineVec[0] == "Read_1_Count")
		{
			read1_count = boost::lexical_cast<uint64_t>(lineVec[2]);
		}
		
		//"Read_2_Count = 145785569"
		if (lineVec[0] == "Read_2_Count")
		{
			read2_count = boost::lexical_cast<uint64_t>(lineVec[2]);
		}
		
		//"MaxInDel_Length = 3"
		if(lineVec[0] == "MaxInDel_Length")
		{
			InDel_max_len = atoi(lineVec[2].c_str());
		}
	}
	
	infile.close();
}

//load InDel-error profile and get the InDel distribution table.
string load_InDel_error_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num2, int InDel_max_len, uint64_t read1_count,
 uint64_t read2_count, int* InDel_num, double** InDel_error_matrix)
{
	string matrix_file;
	
	if(InputParameter.InDel_error_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + INDEL_ERROR_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.InDel_error_profile;
	}

	igzstream MatrixFile;
	MatrixFile.open(matrix_file.c_str());
	if ( ! MatrixFile )
	{	cerr << "fail to open input file " << matrix_file <<", please make sure statistics file place in program directory!"<< endl;
		exit(-1);
	}
	
	string str_line;
	cerr <<"Start to construct InDel-error matrix..."<<endl
		<<"Loading file: "<<matrix_file<<endl;
	
	uint64_t read1_ins_num = 0;
	uint64_t read1_del_num = 0;
	uint64_t read2_ins_num = 0;
	uint64_t read2_del_num = 0;
	
	int simulate_cycle = 0;
	bool isEnd = 0;
	while(getline(MatrixFile, str_line, '\n'))
	{
		if(isEnd){break;}
		//[InDel]
    //Cycle   -3      -2      -1      1       2       3
		//1       0       6       2       4       0       0
    //2       0       2       13      32      5       0
    //3       0       5       224     31      3       0
    //4       0       1       278     34      7       0
		//.....
		if(str_line == "[InDel]")
		{
			
			while(getline(MatrixFile, str_line, '\n'))
			{
				if(str_line == "" || str_line[0] == '#' || str_line[0] == 'C'){continue;}
				if(str_line == "<<END"){isEnd = 1; break;}
    		vector<string> str_line_tokens;
    		boost::split(str_line_tokens,str_line, boost::is_any_of(" \t\n"), boost::token_compress_on);
    		if(str_line_tokens.size() != InDel_max_len*2 + 1)
    		{
    			cerr<<"Error: Length of InDel are not consistent, please check the InDel-error profile is in right format!"<<endl;
    			exit(-1);
    		}
    		int current_cycle = atoi(str_line_tokens[0].c_str());
    		
    		if(current_cycle <= InputParameter.Read_length || (current_cycle > statistical_Cycle_num2/2 && current_cycle <= statistical_Cycle_num2/2+InputParameter.Read_length))
    		{
    			simulate_cycle = current_cycle;
    			if(simulate_cycle >= statistical_Cycle_num2/2){
    				simulate_cycle = InputParameter.Read_length + (current_cycle - statistical_Cycle_num2/2);
    			}

					uint64_t indel_sum = 0;
					for(int i = 0; i < InDel_max_len*2;i++)
					{
						string num = str_line_tokens[i+1];
						uint64_t current_num = boost::lexical_cast<uint64_t>(num);
						
						int j = i;
						if(j >= InDel_max_len){j++;}  //insertion 
							
						if((simulate_cycle <= InputParameter.Read_length) && (InputParameter.Read_length - simulate_cycle < InDel_num[j])
							|| (simulate_cycle > InputParameter.Read_length) && (InputParameter.Read_length - (simulate_cycle - InputParameter.Read_length) < InDel_num[j]))
						{//modify read-end insertion base-num
							current_num = 0;
						}
					
						InDel_error_matrix[simulate_cycle-1][j] = current_num;
						
						
						if(simulate_cycle <= InputParameter.Read_length){  //read1 indel
  						if(i<InDel_max_len)
  						{
  							read1_del_num += abs(InDel_num[j]) * current_num;
  						}else{
  							read1_ins_num += InDel_num[j] * current_num;
  						}
						}else{  //read2 indel
  						if(i<InDel_max_len)
  						{
  							read2_del_num += abs(InDel_num[j]) * current_num;
  						}else{
  							read2_ins_num += InDel_num[j] * current_num;
  						}
						}
						
						
						indel_sum+=current_num;
					}
					
					//count 0-indel number
					if(simulate_cycle < InputParameter.Read_length){
						InDel_error_matrix[simulate_cycle-1][InDel_max_len] = read1_count - indel_sum;
					}else{
						InDel_error_matrix[simulate_cycle-1][InDel_max_len] = read2_count - indel_sum;
					}
    		}
			}
		}else{
			continue;
		}
	}
	
	double read1_ins_rate = double(read1_ins_num)/double(read1_count*InputParameter.Read_length);
	double read1_del_rate = double(read1_del_num)/double(read1_count*InputParameter.Read_length);
	double read2_ins_rate = double(read2_ins_num)/double(read2_count*InputParameter.Read_length);
	double read2_del_rate = double(read2_del_num)/double(read2_count*InputParameter.Read_length);
	
  cout << "\nIn InDel-error profile:\n"
		<< "\tprofile_Cycle_num: " << statistical_Cycle_num2 << endl
  	<< "\tread1 count: " << read1_count << endl
  	<< "\tread2 count: " << read2_count << endl
  	<< "\tlength of max InDel: " << InDel_max_len << endl
  	<< "\tinsertion-base rate of "<<InputParameter.Read_length << "-bp read1: " << read1_ins_rate << endl
  	<< "\tdeletion-base rate of "<<InputParameter.Read_length << "-bp read1: " << read1_del_rate << endl
  	<< "\tinsertion-base rate of "<<InputParameter.Read_length << "-bp read2: " << read2_ins_rate << endl
  	<< "\tdeletion-base rate of "<<InputParameter.Read_length << "-bp read2: " << read2_del_rate << endl;

	for(int j=0; j<statistical_Cycle_num2; j++)
	{
		//the cumulative rate
		double sum = 0;
		for(int k=0; k<InDel_max_len*2+1; k++)
		{
			sum += InDel_error_matrix[j][k];
		}
		if(sum == 0)
		{
			for(int k = 0; k < InDel_max_len*2+1; k++)
			{
				InDel_error_matrix[j][k] = 0;
			}
		}else{
			double accumulate_value = 0;
			for(int k = 0; k < InDel_max_len*2+1; k++)
			{
				accumulate_value += InDel_error_matrix[j][k];
				InDel_error_matrix[j][k] = accumulate_value/sum;
			}
		}
	}
	
  MatrixFile.close();
  
  cerr <<"Have finished constructing InDel-error simulation matrix"<<endl;
  
  return matrix_file;
}

//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile1(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, 
	int quality_num, double statistical_average_error_rate, double*** simulation_matrix)
{
	string matrix_file;
	if(InputParameter.BaseCalling_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + BASE_CALLING_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.BaseCalling_profile;
	}

	igzstream MatrixFile;
	MatrixFile.open(matrix_file.c_str());
	if ( ! MatrixFile )
	{	cerr << "fail to open input file " << matrix_file <<", make sure statistics file place in program directory!"<< endl;
		exit(-1);
	}
	
	string str_line;
	cerr <<"Start to construct simulation matrix..."<<endl
		<<"Loading file: "<<matrix_file<<endl;
	
	
	//user_error_rate/profile_error_rate
	double E_ratio = 0;
	if(InputParameter.Error_rate == -1){
		E_ratio = 1;
	}else{
		E_ratio = InputParameter.Error_rate/statistical_average_error_rate;
	}
	
	int simulate_cycle = 0;
	bool isEnd = 0;
	while(getline(MatrixFile, str_line, '\n'))
	{
		if(isEnd){break;}
		if(str_line == "[QTransMatrix]")
		{
			while(getline(MatrixFile, str_line, '\n'))
			{
				if(str_line == "" || str_line[0] == '#'){continue;}
				if(str_line == "<<END"){isEnd = 1; break;}
    		vector<string> str_line_tokens;
    		boost::split(str_line_tokens,str_line, boost::is_any_of(" \t\n"), boost::token_compress_on);
//    		char ref_base = str_line_tokens[0][0];
    		int current_cycle = atoi(str_line_tokens[0].c_str());
    		
    		if(current_cycle == 1){simulate_cycle = 0;}
    		if(current_cycle <= InputParameter.Read_length || (current_cycle > statistical_Cycle_num/2 && current_cycle <= statistical_Cycle_num/2+InputParameter.Read_length))
    		{
    			int pre_Q = atoi(str_line_tokens[1].c_str());
    			
    			simulate_cycle = current_cycle;
    			if(simulate_cycle >= statistical_Cycle_num/2){
    				simulate_cycle = InputParameter.Read_length + (current_cycle - statistical_Cycle_num/2);
    			}

    			for(int j = 0; j < quality_num; j++)
    			{
    				string num = str_line_tokens[j + 2];
    				uint64_t current_num = boost::lexical_cast<uint64_t>(num);

						simulation_matrix[simulate_cycle-1][pre_Q][j] += current_num;
					}
    		}
			}
		}else{
			continue;
		}
	}

	for(int j=0; j<InputParameter.Read_length * 2; j++)
	{
		for(int k=0; k<quality_num; k++)
		{
  		//the cumulative rate
  		double sum = 0;
  		for(int l = 0; l < quality_num; l++)
  		{
  			sum += simulation_matrix[j][k][l];
  		}
  		if(sum == 0){
  			for(int l = 0; l < quality_num; l++)
  			{
  				simulation_matrix[j][k][l] = 0;
  			}	
  		}else{
  			double accumulate_value = 0;
  			for(int l = 0; l < quality_num; l++)
  			{
  				accumulate_value += simulation_matrix[j][k][l];
  				simulation_matrix[j][k][l] = accumulate_value/sum;
  			}
  			
  		}
    }
  } 
  
  MatrixFile.close();
  
  cerr <<"Have finished constructing Base-calling simulation matrix1"<<endl;
  
  return matrix_file;
}

//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile2(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, 
	int quality_num, double statistical_average_error_rate, double**** simulation_matrix, double*** First_cycle_matrix)
{
	string matrix_file;
	if(InputParameter.BaseCalling_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + BASE_CALLING_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.BaseCalling_profile;
	}

	igzstream MatrixFile;
	MatrixFile.open(matrix_file.c_str());
	if ( ! MatrixFile )
	{	cerr << "fail to open input file " << matrix_file <<", make sure statistics file place in program directory!"<< endl;
		exit(-1);
	}
	
	string str_line;
	cerr <<"Start to construct simulation matrix..."<<endl
		<<"Loading file: "<<matrix_file<<endl;
	
	
	//user_error_rate/profile_error_rate
	double E_ratio = 0;
	if(InputParameter.Error_rate == -1){
		E_ratio = 1;
	}else{
		E_ratio = InputParameter.Error_rate/statistical_average_error_rate;
	}
	
	int simulate_cycle = 0;
	bool isEnd = 0;
	while(getline(MatrixFile, str_line, '\n'))
	{
		if(isEnd){break;}
		if(str_line == "[DistMatrix]")
		{
			while(getline(MatrixFile, str_line, '\n'))
			{
				if(str_line == "" || str_line[0] == '#'){continue;}
				if(str_line == "<<END"){isEnd = 1; break;}
    		vector<string> str_line_tokens;
    		boost::split(str_line_tokens,str_line, boost::is_any_of(" \t\n"), boost::token_compress_on);
    		char ref_base = str_line_tokens[0][0];
    		int current_cycle = atoi(str_line_tokens[1].c_str());
    		
    		if(current_cycle == 1){simulate_cycle = 0;}
    		if(current_cycle <= InputParameter.Read_length || (current_cycle > statistical_Cycle_num/2 && current_cycle <= statistical_Cycle_num/2+InputParameter.Read_length))
    		{
    			simulate_cycle++;

        	//convert profile error matrix to simulation error matrix
      		for(int i = 0; i < quality_num; i++)
      		{
      			uint64_t count_sum = 0;
      			
      			for(int j = 0; j < seq_Base_num; j++)
      			{
      				string num = str_line_tokens[j * quality_num + i + 2];
      				uint64_t current_num = boost::lexical_cast<uint64_t>(num);
      				count_sum+=current_num;
      				
      				if( j == alphabet3[ref_base])  // ACGT  
      				{
      					simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j] = current_num;

      					if(simulate_cycle == 1){
      						First_cycle_matrix[alphabet3[ref_base]][0][j * quality_num + i] =  current_num;
      					}
      					if(simulate_cycle == statistical_Cycle_num/2 + 1){
      						First_cycle_matrix[alphabet3[ref_base]][1][j * quality_num + i] =  current_num;
      					}
      				}else{
      					simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j] = current_num;

      					if(simulate_cycle == 1){
      						First_cycle_matrix[alphabet3[ref_base]][0][j * quality_num + i] =  current_num;
      					}
      					if(simulate_cycle == statistical_Cycle_num/2 + 1){
      						First_cycle_matrix[alphabet3[ref_base]][1][j * quality_num + i] =  current_num;
      					}
      				}
      			}
      		}
        
      		//the cumulative rate
      		if(simulate_cycle == 1){
        		double sum = 0;
        		for(int i = 0; i < seq_Base_num*quality_num; i++)
        		{
        				sum += First_cycle_matrix[alphabet3[ref_base]][0][i];
        		}
        		//cerr<<sum<<endl;
        		if(sum == 0){
        			for(int i = 0; i < seq_Base_num*quality_num; i++)
        			{
        				First_cycle_matrix[alphabet3[ref_base]][0][i] = 0;
        			}	
        		}else{
        			double accumulate_value = 0;
        			for(int i = 0; i < seq_Base_num*quality_num; i++)
        			{
        				accumulate_value += First_cycle_matrix[alphabet3[ref_base]][0][i];
        				First_cycle_matrix[alphabet3[ref_base]][0][i] = accumulate_value/sum;
        			}	
        		}
      		}else if(simulate_cycle == statistical_Cycle_num/2 + 1)
      		{
        		double sum = 0;
        		for(int i = 0; i < seq_Base_num*quality_num; i++)
        		{
        				sum += First_cycle_matrix[alphabet3[ref_base]][1][i];
        		}
        		//cerr<<sum<<endl;
        		if(sum == 0){
        			for(int i = 0; i < seq_Base_num*quality_num; i++)
        			{
        				First_cycle_matrix[alphabet3[ref_base]][1][i] = 0;
        			}	
        		}else{
        			double accumulate_value = 0;
        			for(int i = 0; i < seq_Base_num*quality_num; i++)
        			{
        				accumulate_value += First_cycle_matrix[alphabet3[ref_base]][1][i];
        				First_cycle_matrix[alphabet3[ref_base]][1][i] = accumulate_value/sum;
        			}	
        		}
      		}
      		
      		for(int i = 0; i < quality_num; i++)
      		{
      			double sum = 0;
        		for(int j = 0; j < seq_Base_num; j++)
        		{
        			sum += simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j];
        		}
        		if(sum == 0){
        			for(int j = 0; j < seq_Base_num; j++)
        			{
        				simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j] = 0;
        			}	
        		}else{
        			double accumulate_value = 0;
        			for(int j = 0; j < seq_Base_num; j++)
        			{
        				accumulate_value += simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j];
        				simulation_matrix[alphabet3[ref_base]][simulate_cycle-1][i][j] = accumulate_value/sum;
        			}	
        		}
      		}
    		}
			}
		}else{
			continue;
		}
	}
  
  MatrixFile.close();
  
  cerr <<"Have finished constructing Base-calling simulation matrix2"<<endl;
  
  return matrix_file;
}

//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile3(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, 
	int quality_num, double statistical_average_error_rate, double*** simulation_matrix1, double**** simulation_matrix2)
{
	string matrix_file;
	if(InputParameter.BaseCalling_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		string directory_path = exe_path.substr(0,index);
  		matrix_file = directory_path + BASE_CALLING_PROFILE;
  	}
	}else{
		matrix_file = InputParameter.BaseCalling_profile;
	}

	igzstream MatrixFile;
	MatrixFile.open(matrix_file.c_str());
	if ( ! MatrixFile )
	{	cerr << "fail to open input file " << matrix_file <<", make sure statistics file place in program directory!"<< endl;
		exit(-1);
	}
	
	string str_line;
	cerr <<"Start to construct simulation matrix..."<<endl
		<<"Loading file: "<<matrix_file<<endl;
	
	
	//user_error_rate/profile_error_rate
	double E_ratio = 0;
	if(InputParameter.Error_rate == -1){
		E_ratio = 1;
	}else{
		E_ratio = InputParameter.Error_rate/statistical_average_error_rate;
	}
	
	int simulate_cycle = 0;
	bool isEnd = 0;
	while(getline(MatrixFile, str_line, '\n'))
	{
		if(isEnd){break;}
		if(str_line == "[DistMatrix]")
		{
			while(getline(MatrixFile, str_line, '\n'))
			{
				if(str_line == "" || str_line[0] == '#'){continue;}
				if(str_line == "<<END"){isEnd = 1; break;}
    		vector<string> str_line_tokens;
    		boost::split(str_line_tokens,str_line, boost::is_any_of(" \t\n"), boost::token_compress_on);
    		char ref_base = str_line_tokens[0][0];
    		int current_cycle = atoi(str_line_tokens[1].c_str());
    		
    		if(current_cycle == 1){simulate_cycle = 0;}
    		if(current_cycle <= InputParameter.Read_length || (current_cycle > statistical_Cycle_num/2 && current_cycle <= statistical_Cycle_num/2+InputParameter.Read_length))
    		{
    			simulate_cycle++;
    			uint64_t count_sum = 0;
    			
      		//convert profile error matrix to simulation error matrix
      		for(int i = 0; i < seq_Base_num; i++)
      		{
      			for(int j = 0; j < quality_num; j++)
      			{
      				string num = str_line_tokens[i * quality_num + j + 2];
      				uint64_t current_num = boost::lexical_cast<uint64_t>(num);
      				count_sum+=current_num;
							simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][j] += current_num;
							simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][j][i] += current_num;
      			}
      		}
      		
      		//convert profile error matrix to simulation error matrix
      		for(int i = 0; i < quality_num; i++)
      		{
      			for(int j = 0; j < seq_Base_num; j++)
      			{
      				string num = str_line_tokens[j * quality_num + i + 2];
      				uint64_t current_num = boost::lexical_cast<uint64_t>(num);
      				count_sum+=current_num;
							simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][i] += current_num;
							simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][i][j] = current_num;
      			}
      		}
      		
      		//the cumulative rate
      		double sum = 0;
      		for(int i = 0; i < quality_num; i++)
      		{
      			sum += simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][i];
      		}
      		//cerr<<sum<<endl;
      		if(sum == 0){
      			for(int i = 0; i < quality_num; i++)
      			{
      				simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][i] = 0;
      			}	
      		}else{
      			double accumulate_value = 0;
      			for(int i = 0; i < quality_num; i++)
      			{
      				accumulate_value += simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][i];
      				simulation_matrix1[alphabet3[ref_base]][simulate_cycle-1][i] = accumulate_value/sum;
      			}	
      		}
      		
      		//the cumulative rate
      		
      		for(int i = 0; i < quality_num; i++)
      		{
      			double sum = 0;
        		for(int j = 0; j < seq_Base_num; j++)
        		{
        			sum += simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][i][j];
        		}
        		if(sum == 0){
        			for(int j = 0; j < seq_Base_num; j++)
        			{
        				simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][i][j] = 0;
        			}	
        		}else{
        			double accumulate_value = 0;
        			for(int j = 0; j < seq_Base_num; j++)
        			{
        				accumulate_value += simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][i][j];
        				simulation_matrix2[alphabet3[ref_base]][simulate_cycle-1][i][j] = accumulate_value/sum;
        			}	
        		}
      		}
    		}
			}
		}else{
			continue;
		}
	}
	
  MatrixFile.close();
  cerr <<"Have finished constructing Base-calling simulation matrix"<<endl;
  return matrix_file;
}

string load_GC_depth_profile (PARAMETER InputParameter, string exe_path, double* GC_bias_abundance)
{
	string depth_file;
	
	if(InputParameter.GC_depth_profile == ""){
  	int index = exe_path.find_last_of('/');
  	if(index == -1){
  		cerr<<"Error: program path wrong!"<<endl;
  	}
  	else{
  		//The default GC bias profile are determined based on the twice length of read
  		string GC_depth_profile_name;
			int window_size = InputParameter.Read_length * 2;
    	if(window_size < 125){
    		GC_depth_profile_name = GC_DEPTH100_PROFILE;
    	}else{
    		if(window_size < 175){
    			GC_depth_profile_name = GC_DEPTH150_PROFILE;
    		}else{
    			GC_depth_profile_name = GC_DEPTH200_PROFILE;
    		}
    	}
  		string directory_path = exe_path.substr(0,index);
  		depth_file = directory_path + GC_DEPTH_PROFILE_PATH+GC_depth_profile_name;
  	}
	}else{
		depth_file = InputParameter.GC_depth_profile;
	}
	
  igzstream infile;
  infile.open(depth_file.c_str());
  if ( ! infile )
	{
		cerr << "fail to open input file" << depth_file << endl;
	}else{
		cerr << "Loading the GC bias profile: "<< depth_file <<endl;
	}
	uint64_t total_count_sum = 0;
	uint64_t total_error_sum = 0;

	string lineStr;
	vector<double> GC_ratio_vec;
	vector<double> depth_vec;
	while (getline( infile, lineStr, '\n' ))
	{ if (lineStr[0] == '#' || lineStr == "" || lineStr[0] == ' ')
		{ 
			continue;
		}
		else
		{ 
			vector<string> lineVec;
			//#GC%    RefCnt  DepthCnt       Mean	SmoothedMean    Small   Q1      Mid     Q3      Big     Min     Max     Refcntcal
			//0.5     2285822 1371         0.334179  0.334179   0     0      0       0.2     0.49    0       23.48   0.0945074
			//............
			boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
			GC_ratio_vec.push_back(boost::lexical_cast<double>(lineVec[0]));  //GC%
			depth_vec.push_back(boost::lexical_cast<double>(lineVec[4]));   //SmoothedMean
		}
	}
	
	//find max depth
	double max_depth = 0;
	for(int i = 0; i < depth_vec.size(); i++)
	{
		if(depth_vec[i] > max_depth){
			max_depth = depth_vec[i];
		}
	}
	
	//convert to GC abundance
	for(int i = 1; i < GC_ratio_vec.size(); i++)
	{
		for(int j = int(GC_ratio_vec[i-1]); j < int(GC_ratio_vec[i]); j++)
		{
			GC_bias_abundance[j] = depth_vec[i-1]/max_depth;
  	}
  }
  for(int i = int(GC_ratio_vec[GC_ratio_vec.size()-1]); i <=100; i++)
  {
  	GC_bias_abundance[i] = depth_vec[depth_vec.size()-1]/max_depth;
  }
  
  cout<<"\nFor simulating GC bias:"<<endl
  	<<"\tGC%\tabundance"<<endl;
  for(int i = 0; i <=100; i++)
  {
  	cout<<"\t"<<i<<"\t"<<GC_bias_abundance[i]<<endl;
  }
  
  infile.close();
  
  cerr <<"Have finished constructing GC bias simulation matrix"<<endl;
  
  return depth_file;
}


