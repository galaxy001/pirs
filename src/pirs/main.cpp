#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <string.h>

using namespace std;

const char *PROGRAM = "pirs(Profile-based Illumina pair-end Reads Simulator)";
const char *AUTHOR  = "Jianying Yuan (BGI-Shenzhen)";
const char *VERSION = "v1.1.1";                      
const char *CONTACT = "yuanjianying@genomics.org.cn";

extern int simulate_Illumina_reads(int argc, char *argv[]);
extern int simulate_diploid_genome(int argc, char *argv[]);


static void display_usage();

int main(int argc, char *argv[])
{
	argc--;
	argv++;

	if(argc==0){
		display_usage();
		return 0;
	}
	
	if(strcmp("diploid",argv[0])==0)
		simulate_diploid_genome(argc,argv);
	else if(strcmp("simulate",argv[0])==0)
		simulate_Illumina_reads(argc,argv);
	else
		display_usage();
	
	return 0;
}

static void display_usage()
{
	cerr<<"\nProgram:\t"<< PROGRAM <<"\n";
  cerr<<"Version:\t"<<VERSION<<"\n";
  cerr<<"Author:\t\t"<<AUTHOR<<"\n";
  cerr<<"Contact: \t"<<CONTACT<<"\n";
	cerr<<"Compile Date:\t"<<__DATE__<<" time: "<< __TIME__<<"\n";

	cerr<<"\nUsage: pirs <command> [option]\n";
	cerr<<"    diploid     generate diploid genome.\n";
	cerr<<"    simulate    simulate Illumina reads.\n\n";
}
