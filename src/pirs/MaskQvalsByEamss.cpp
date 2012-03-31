/**
 ** Copyright (c) 2008-2010 Illumina, Inc.
 **
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 **
 ** This file is part of the Consensus Assessment of Sequence And VAriation
 ** (CASAVA) software package.
 **
 ** \file MaskQvalsByEamss.cpp
 **
 ** \brief Masks quality values of a read, using the EAMSS.
 **
 ** \author Come Raczy
 **/

#include <boost/assign.hpp>

#include "demultiplex/MaskQvalsByEamss.hh"

namespace casava
{
namespace demultiplex
{

const int MaskQvalsByEamss::lowScore = 1;
const int MaskQvalsByEamss::mediumScore = 0;
const int MaskQvalsByEamss::highScore = -2;
const char MaskQvalsByEamss::mediumThreshold = 'O';
const char MaskQvalsByEamss::highThreshold = '^';
const int MaskQvalsByEamss::minScore = 1;
const std::vector<std::string> MaskQvalsByEamss::motifList = std::vector<std::string>(0);

/**
 ** \brief Mask quality using the EAMSS algorithm
 **/
void MaskQvalsByEamss::operator()(std::string &qValues,
                                  const std::string &baseCalls) const
{
    std::pair<int, int> tmp = eamss(qValues);
    const int score = tmp.first;
    int position = tmp.second;
    const bool filterRead = (score >= minScore);
    if (!filterRead)
    {
        return;
    }
    // look for troublemaker motifs
    int extendedPosition = position;
    for (std::vector<std::string>::const_iterator motif = motifList.begin(); motifList.end() != motif; ++motif)
    {
        const int troublemakerPosition = findStr(baseCalls, *motif,
                                                 std::max(0, position - 15),
                                                 std::min(static_cast<int>(baseCalls.length()), position));
        if (-1 < troublemakerPosition)
        {
            extendedPosition = std::min(extendedPosition, troublemakerPosition);
        }
    }
    position = std::min(position, extendedPosition);
    // look for extended run of polyG, with at least 90% G bases
    unsigned int numG = 0;
    unsigned int numBases = 0;
    bool maskPolyG = false;
    int maskStart = position;
    for (int curPos = position; 0 <= curPos; --curPos)
    {
        if ('G' == baseCalls[curPos])
        {
            ++numG;
        }
        ++numBases;
        if (10 > numBases)
        {
            continue;
        }
        const float gFrac = static_cast<float>(numG)/static_cast<float>(numBases);
        if (gFrac >= 0.9f && 'G' == baseCalls[curPos]) // only start masking at G
        {
            maskPolyG = true;
            maskStart = curPos;
        }
        else if (gFrac < 0.9f)
        {
            break;
        }
    }
    if (maskPolyG)
    {
        position = maskStart;
    }
    for (int idx = position; static_cast<long int>(qValues.length()) > idx; ++idx)
    {
        qValues[idx] = 'B';
    }
}

/**
 ** \brief Compute the position with the highest EAMSS score
 **
 ** \return the pair (bestScore, bestPosition) or (highScore-1, -1) if qValues
 ** is empty
 **/
std::pair<int, int> MaskQvalsByEamss::eamss(const std::string &qValues) const
{
    int curScore = 0;
    // initialize the bestscore to something lower than the first value of curScore
    int bestScore = std::min(std::min(highScore, mediumScore), lowScore) - 1;
    int bestPosition = -1;
    for (int idx = qValues.length() - 1; 0 <= idx; --idx)
    {
        if (qValues[idx] >= highThreshold)
        {
            curScore += highScore;
        }
        else if (qValues[idx] >= mediumThreshold)
        {
            curScore += mediumScore;
        }
        else
        {
            curScore += lowScore;
        }
        if (curScore >= bestScore)
        {
            bestScore = curScore;
            bestPosition = idx;
        }
    }
    return std::pair<int, int>(bestScore, bestPosition);
}

int MaskQvalsByEamss::findStr(const std::string &targetString, const std::string &queryString, int start, int stop) const
{
    const std::string substring = targetString.substr(start, stop - start + 1);
    const size_t result = substring.find(queryString);
    return (substring.npos == result) ? -1 : (result+start);
}

} // namespace demultiplex
} // namespace casava
