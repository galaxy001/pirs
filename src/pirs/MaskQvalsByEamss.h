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

// Modified for pIRS: Pass qValues and baseCalls as C arrays (on the stack)

#ifndef CASAVA_DEMULTIPLEX_MASK_QVALS_BY_EAMSS_H
#define CASAVA_DEMULTIPLEX_MASK_QVALS_BY_EAMSS_H

#include <vector>

namespace casava
{
namespace demultiplex
{

enum QualityMaskMode {
	MODE_NONE,
	MODE_QUALITY,
	MODE_LOWERCASE,
};

const char *get_quality_mask_mode_name(enum QualityMaskMode mode);

class MaskQvalsByEamss
{
public:
	int operator()(std::vector<char> &qValues, std::vector<char> &baseCalls,
		       enum QualityMaskMode mode) const;
private:
	static const int lowScore;
	static const int mediumScore;
	static const int highScore;
	static const char mediumThreshold;
	static const char highThreshold;
	static const int minScore;
	static const std::vector<std::string> motifList;
	std::pair<int, int> eamss(const std::vector<char> &qValues) const;
	int findStr(const char targetString[], const std::string &queryString, int start, int stop) const;
};

} // namespace demultiplex
} // namespace casava

#endif // #ifndef CASAVA_DEMULTIPLEX_MASK_QVALS_BY_EAMSS_H
