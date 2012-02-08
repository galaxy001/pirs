#include <string>
#include <sstream>
#include <iostream>
#include "self_util.h"

using namespace std;

char alphabet[128] =
{
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

char bases[5] =
{
    'A', 'C', 'G', 'T', 'N'
};

char c_bases[5] =
{
    'T', 'G', 'C', 'A', 'N'
};

uint64_t seq2bit(string &kseq)
{
    uint64_t kbit=0;
    
    for(int i=0; i<kseq.size(); i++)
    {
        kbit=(kbit<<2)|alphabet[kseq[i]];
    }
    
    return kbit;
}

string bit2seq(uint64_t kbit, int kmerSize)
{
    string kseq;

    for(int i=0; i<kmerSize; i++)
    {
        kseq.push_back(bases[(kbit>>(kmerSize-1-i)*2)&0x3]);
    }

    return kseq;
}

bool check_seq(string &seq)
{
    bool is_good = true;
    for (int i = 0; i < seq.size(); i++)
    {
        if (alphabet[seq[i]] == 4)
        {
            is_good = false;
            break;
        }
    }
    
    return is_good;
}


void reverse_complement (string &in_str, string &out_str)
{
    for (int64_t i=in_str.size()-1; i>=0; i--)
    {
        out_str.push_back(c_bases[alphabet[in_str[i]]]);
    }
}

void TrimRight(string& str)
{
    int pos = str.length() - 1;
    if(pos >= 0)
    {
        char c = str[pos];
        while((c == '\t') || (c == '\n') || (c == ' '))
        {
            str = str.substr(0, pos);
            pos = str.length() - 1;
            if(pos < 0)
            {
                break;
            }
            c = str[pos];
        }
    }
}

void TrimLeft(string& str)
{
    if(str.length() > 0)
    {
        char c = str[0];
        while((c == '\t') || (c == '\n') || (c == ' '))
        {
            str = str.substr(1);
            if(str.length() == 0)
            {
                break;
            }
            c = str[0];
        }
    }
}

vector<string> splitString(string str, string delims)
{
    vector<string> vec;
    string temp;
    //const string delims(" \t\n");
    string::size_type begIdx, endIdx;
    begIdx = str.find_first_not_of(delims);
    while(begIdx != string::npos)
    {
        endIdx = str.find_first_of(delims, begIdx);
        if(endIdx == string::npos)
        {
            endIdx == str.length();
        }

        vec.push_back(str.substr(begIdx, endIdx-begIdx));
        begIdx = str.find_first_not_of(delims, endIdx);
    }
    return vec;
}

vector<int> splitStringToInt(string str, string delims)
{
    vector<int> vec;
    string temp;
    string::size_type begIdx, endIdx;
    begIdx = str.find_first_not_of(delims);
    while(begIdx != string::npos)
    {
        endIdx = str.find_first_of(delims, begIdx);
        if(endIdx == string::npos)
        {
            endIdx == str.length();
        }
        vec.push_back(toInt(str.substr(begIdx, endIdx-begIdx)));
        begIdx = str.find_first_not_of(delims, endIdx);
    }
    return vec;
}

int toInt(string str)
{
    int i;
    stringstream ss;
    //int prec = numeric_limits<int>::digits10;
    //ss.precision(prec);
    ss << str;
    ss >> i;
    return i;
}

uint64_t toUint64(string str)
{
    int i;
    stringstream ss;
    //int prec = numeric_limits<int>::digits10;
    //ss.precision(prec);
    ss << str;
    ss >> i;
    return i;
}

char toUpper(const char& ch) 
{ 
    if(islower(ch))
    {
        return ch & 0x5F; 
    }
}

char toLower(const char& ch) 
{ 
    if(isupper(ch))
    {
        return ch | 0x20; 
    }
} 
string ToUpper(string str)
{
    string temp = "";
    for(int i=0; i<str.length(); ++i)
    {
        temp += toUpper(str[i]);
    }
    return temp;
}

string ToLower(string str)
{
    string temp = "";
    //transform (str.begin(),str.end(), temp.begin(), toLower); 
    for(int i=0; i<str.length(); ++i)
    {
        temp += toUpper(str[i]);
    }
    return temp;
}
    
