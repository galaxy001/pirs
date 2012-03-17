#ifndef __LOADFILE_H_
#define __LOADFILE_H_

#include "global.h"

using namespace std;

#define BASE_CALLING_PROFILE "/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz"
#define INDEL_ERROR_PROFILE "/Profiles/InDel_Profiles/phixv2.InDel.matrix"
#define GC_DEPTH_PROFILE_PATH "/Profiles/GC-depth_Profiles/"
#define GC_DEPTH100_PROFILE "humNew.gcdep_100.dat"
#define GC_DEPTH150_PROFILE "humNew.gcdep_150.dat"
#define GC_DEPTH200_PROFILE "humNew.gcdep_200.dat"

//check and open the outfile file
//void set_and_check_file(PARAMETER InputParameter, igzstream &infile, igzstream &infile2, ofstream &outfile1, ofstream &outfile2, ogzstream &gz_outfile1, ogzstream &gz_outfile2, ofstream &insert_log, ofstream &error_log);
//check and open the outfile file
void set_and_check_file(PARAMETER InputParameter, igzstream &infile, igzstream &infile2, ofstream &outfile1,	ofstream &outfile2, 
	ogzstream &gz_outfile1, ogzstream &gz_outfile2, ofstream &insert_log, ofstream &error_log, ogzstream &infor_outfile);
	
//get the attribute of error profile
void preview_BaseCalling_profile (PARAMETER InputParameter, string exe_path, int &ref_Base_num, int &statistical_Cycle_num, int &seq_Base_num, int &quality_num, double &statistical_average_error_rate);

//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, int quality_num, double statistical_average_error_rate, double*** simulation_matrix);

//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int ref_Base_num, int simulate_Cycle_num, int seq_Base_num, int quality_num, double statistical_average_error_rate, double*** simulation_matrix2);
	
//read in quality distribution file and get the quality distribution table.
string load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, 
	int quality_num, double statistical_average_error_rate, double**** simulation_matrix1, double*** First_cycle_matrix);
	
string load_GC_depth_profile (PARAMETER InputParameter, string exe_path, double* GC_bias_abundance);

//load InDel-error profile and get the InDel distribution table.
void preview_InDel_error_profile (PARAMETER InputParameter, string exe_path, int &Statistical_Cycle_num2, int &InDel_max_len, uint64_t &read1_count, uint64_t &read2_count);

//load InDel-error profile and get the InDel distribution table.
string load_InDel_error_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num2, int InDel_max_len, uint64_t read1_count,
 uint64_t read2_count, int* InDel_num, double** InDel_error_matrix);
 
#endif
