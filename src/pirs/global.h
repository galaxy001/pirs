#ifndef GLOBAL_H
#define GLOBAL_H

using namespace std;
using namespace boost; 


typedef struct{
  int Read_length;
  int Insertsize_mean;
  int Insertsize_sd;
	int Is_simulate_InDel;
  int Is_cyclization;
  int Is_simulate_GC_bias;
  int Is_simulate_quality;
  int Q_Mode;
  int Q_shift;
  int Output_type;
  double Coverage;
  double Error_rate;
	string Input_ref1;
  string Input_ref2;
  string BaseCalling_profile;
  string GC_depth_profile;
  string InDel_error_profile;
  string Output_prefix;
}PARAMETER;

#endif
