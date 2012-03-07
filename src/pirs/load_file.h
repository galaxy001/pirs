#ifndef __LOADFILE_H_
#define __LOADFILE_H_

#include "global.h"

using namespace std;

//check and open the outfile file
void set_and_check_file(PARAMETER InputParameter, igzstream &infile, igzstream &infile2, ofstream &outfile1, ofstream &outfile2, ogzstream &gz_outfile1, ogzstream &gz_outfile2, ofstream &insert_log, ofstream &error_log);

//get the attribute of error profile
void preview_BaseCalling_profile (PARAMETER InputParameter, string exe_path, int &ref_Base_num, int &statistical_Cycle_num, int &seq_Base_num, int &quality_num, double &statistical_average_error_rate);

//read in quality distribution file and get the quality distribution table.
void load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, int quality_num, double statistical_average_error_rate, double*** simulation_matrix);

//read in quality distribution file and get the quality distribution table.
void load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int ref_Base_num, int simulate_Cycle_num, int seq_Base_num, int quality_num, double statistical_average_error_rate, double*** simulation_matrix2);
	
//read in quality distribution file and get the quality distribution table.
void load_BaseCalling_profile(PARAMETER InputParameter, string exe_path, int statistical_Cycle_num, int seq_Base_num, 
	int quality_num, double statistical_average_error_rate, double**** simulation_matrix1, double*** First_cycle_matrix);
	
void load_GC_depth_profile (PARAMETER InputParameter, string exe_path, double* GC_bias_abundance);

#endif
