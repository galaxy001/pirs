#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <cstdlib>
#include "gzstream.h"
#include "self_util.h"
#include "stat_soap_coverage.h"

using namespace std;

void usage()
{
    cerr << "Description:" << endl
        << "    it is a program to stat about the gc_depth and base coverage ratio and base_depth with *.soap.coverage  result." << endl
        << "Program:" << endl
        << "Name:           gc_coverage_bias" << endl
        << "Compile Date:   2012-02-23" << endl
        << "Author:         Gan Jun, Yuan JY, Hu Xuesong" << endl
        << "Version:        1.02" << endl
        << "Contact:        ganjun@genomics.org.cn" << endl;

    cerr << "\n\nUsage:\tgc_coverage_bias [options] <*.soap.coverage>" << endl
        << "Option: -r <string>     a reference sequence file about FA format " << endl
        << "        -o <string>     the the prefix about output file" << endl
        << "        -w <string>     the window length[such as:100,200,300] " << endl
        << "        --gcdump        output the gc ratio in the window length " << endl
        << "        --depwindump    output the avg depth in the window length " <<    endl;

    cerr << "\n\nUsage:" << endl
        << "\t./gc_coverage_bias -r test.fa -o test -w 40,50 [--gcdump --depwindump] test.depth" << endl;
    exit(EXIT_FAILURE);
}

/*
map<string, string> DealReference(string str_ref_file_name);
map<string, vector<vector<uint64_t, uint64_t, uint64_t> >, uint64_t> DealSoapCoverage()
*/


int main(int argc, char* argv[])
{
    if(argc < 8)
    {
        usage();
    }

    string str_argv;
    string str_ref_file_name;
    string str_width;
    string str_output_prefix;
    bool b_gcdump = false;
    bool b_depwindump = false;
    vector<string> vec_soap_file_name;
    vector<string> vec_width;

    for(int i=1; i<argc; ++i)
    {
        str_argv = argv[i];

        if(str_argv.compare("-r") == 0)
        {
            ++i;
            str_ref_file_name = argv[i];
            continue;
        }

        if(str_argv.compare("-o") == 0)
        {
            ++i;
            str_output_prefix = argv[i];
            continue;
        }

        if(str_argv.compare("-w") == 0)
        {
            ++i;
            str_width = argv[i];
            continue;
        }

        if(str_argv.compare("--gcdump") == 0)
        {
            b_gcdump = true;
            continue;
        }

        if(str_argv.compare("--depwindump") == 0)
        {
            b_depwindump = true;
            continue;
        }

        vec_soap_file_name.push_back(argv[i]);
    }

    if(str_ref_file_name.empty())
    {
        usage();
    }else
    {
        igzstream in(str_ref_file_name.c_str());

        if(!in)
        {
            cerr << "can't open the reference file " << str_ref_file_name << ", please check!" << endl;
            exit(EXIT_FAILURE);
        }

        in.close();
    }
    
    if(str_output_prefix.empty())
    {
        usage();
    }

    if(str_width.empty())
    {
        usage();
    }
    else
    {
        vec_width = splitString(str_width, ",");
    }

    //cerr << vec_soap_file_name.size() << endl;
    for(vector<string>::iterator it = vec_soap_file_name.begin(); it != vec_soap_file_name.end(); ++it)
    {
        igzstream in((*it).c_str());
        
        if(!in)
        {
            cerr << "can't open the reference file " << *it << ", please check!" << endl;
            exit(EXIT_FAILURE);
        }

        in.close();
    }

    stat_soap_coverage test(str_ref_file_name, str_output_prefix, vec_soap_file_name, vec_width, b_gcdump, b_depwindump);
}

/*
map<string, string> DealReference(string str_ref_file_name)
{
	igzstream in(str_ref_file.name);
	string line;
	map<string, string> map_reference_base;
	string keyname;
	string sequence;

	while(getline(in, line))
	{
		TrimLeft(line);
		TrimRight(line);

		if(line[0] == '>')
		{
			 if(sequence.length() != 0)
			{
				 map_reference_base[keyname] = sequence;
			}
			int index;
			if(((index = line.find(" ")) == string::npos) || ((index = line.find("\t")) == string::npos) || ((index = line.find("\n")) == string::npos))
			{
				keyname = line.substr(1, index);
			}
			sequence.clear();
			continue;
		}

		sequence += line;
	}

	 if(sequence.length() != 0)
	{
		 map_reference_base[keyname] = sequence;
	}

	in.close();
}

*/
