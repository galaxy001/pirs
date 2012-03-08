#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <string.h>

using namespace std;

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
	cerr<<"\nProgram: pirs(Profile based Illumina pair-end Reads Simulator)\n";
	cerr<<"Version: 1.00\n";
	cerr<<"Author: Jianying Yuan (BGI-Shenzhen)\n";
	cerr<<"Contact: yuanjianying@genomics.org.cn";
	cerr<<"\nUsage: pirs <command> [option]\n";
	cerr<<"    diploid     generate diploid genome.\n";
	cerr<<"    simulate    simulate Illumina reads.\n\n";
}
