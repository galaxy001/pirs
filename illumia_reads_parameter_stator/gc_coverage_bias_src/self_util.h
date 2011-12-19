#ifndef _SELF_UTIL_H_
#define _SELF_UTIL_H_

#include <string>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <cstring>
#include <ctype.h>
#include <algorithm>
using namespace std;

void TrimLeft(string& str);
void TrimRight(string& str);
vector<string> splitString(string str, string delims=" \t\n");
vector<int> splitStringToInt(string str, string delims=" \t\n");
int toInt(string str);
uint64_t toUint64(string str);
template <typename T>
string toStr(T t)
{
    stringstream ss;
    ss << t;
    string s;
    ss >> s;
    return s;
}

string ToUpper(string str);
string ToLower(string str);


#ifndef _BIOINFO_UTIL_
#define _BIOINFO_UTIL_

uint64_t seq2bit(string &kseq);

string bit2seq(uint64_t kbit, int kmerSize);
bool check_seq (string &seq);
void reverse_complement (string &in_str, string &out_str);

#endif

#endif
