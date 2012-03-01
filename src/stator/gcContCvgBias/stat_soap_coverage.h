#ifndef _STAT_SOAP_COVERAGE_H
#define _STAT_SOAP_COVERAGE_H
#include <string>
//#include <pair>
#include <map>
#include <vector>
#include <stdint.h>
using namespace std;

#define MIN_LOESS_COUNT 100

class stat_soap_coverage
{
    private:
        string str_ref_file_name;
        string ref_id_file_name;
        string str_output_prefix;
        vector<string> vec_soap_file_name;
        vector<string> vec_width;
        bool b_gcdump;
        bool b_depwindump;
        map<string, string> map_reference_base;
        map<string, vector<unsigned int> > map_soap_coverage;
        map<int, map<double, vector<double> > > map_width_soap_gc_depth;
        map<int, vector<double> > map_gc_keyname;
        map<string, vector<double> > map_stat_coverage;
        map<int, uint64_t> map_sumwincount;
        map<int, uint64_t> map_sumdepthcount;
        map<int, map<double, uint64_t> >map_wincount;
        map<double, map<string, uint64_t> > map_stat_depth;
        vector<string> vec_chr_keyname;
        map<int, uint64_t> winCountN;

    public:
        stat_soap_coverage(string str_ref_file_name, string ref_id_file_name, string str_output_prefix,
                vector<string> vec_soap_file_name, vector<string> vec_width,
                bool b_gcdump = false, bool b_depwindump = false);
        void DealReference();
        void DealSoapCoverage();
        void DealStat();
        void StatGC();
        void statDepth();
        void StatCoverage();
        void Run();
        ~stat_soap_coverage();
};
#endif
