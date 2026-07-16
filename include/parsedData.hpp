#pragma once
#include "utility.hpp"
#include <string>
#include <vector>
#include <algorithm>

class parsedData
{
public:
    // construction/destruction functions
    parsedData(const AlignLine_Data_t data);
    ~parsedData() = default;

    // get raw data functions
    std::string getQNAME() const { return data.QNAME; }
    std::string getRNAME() const { return data.RNAME; }
    std::string getRawCigar() const { return data.cigar; }
    int getRawFLAG() const { return std::stoi(data.FLAG); }
    std::string getReadSeq() const { return data.read_seq; }
    int getMAPQ() const { return data.MAPQ; }
    int getPosition() const { return data.POS; }
    std::string getQualityScores() const { return data.quality_scores; }
    const options_map_t& getOptionsMap() const { return data.options_map; }

    //tag functions
    bool isTagExists(const std::string &tag) const;
    std::string getTagValue(const std::string &tag) const;

    // CIGAR functions
    std::vector<std::pair<char, std::string>> getCIGARList();
    std::vector<uint> getInsertionList() { return getOpListInCIGAR('I'); }
    std::vector<uint> getSoftClipList() { return getOpListInCIGAR('S'); }
    std::vector<uint> getMatchList() { return getOpListInCIGAR('M'); }
    std::vector<uint> getDeletionList() { return getOpListInCIGAR('D'); }

    // MD functions
    std::string getMDTagValue() const { return getTagValue("MD"); }
    std::vector<std::string> getMDList() const;

    // status flag functions
    bool isStatusFlagSet(int flag) const;
    bool isReadPaired() const { return isStatusFlagSet(READ_PAIRED); }
    bool isProperPair() const { return isStatusFlagSet(PROPER_PAIR); }
    bool isReadUnmapped() const { return isStatusFlagSet(READ_UNMAPPED); }
    bool isMateUnmapped() const { return isStatusFlagSet(MATE_UNMAPPED); }
    bool isReadReverseStrand() const { return isStatusFlagSet(READ_REVERSE_STRAND); }
    bool isMateReverseStrand() const { return isStatusFlagSet(MATE_REVERSE_STRAND); }
    bool isFirstInPair() const { return isStatusFlagSet(FIRST_IN_PAIR); }
    bool isSecondInPair() const { return isStatusFlagSet(SECOND_IN_PAIR); }
    bool isNotPrimaryAlignment() const { return isStatusFlagSet(NOT_PRIMARY_ALIGNMENT); }
    bool isReadFailsQualityChecks() const { return isStatusFlagSet(READ_FAILS_QUALITY_CHECKS); }
    bool isDuplicateRead() const { return isStatusFlagSet(DUPLICATE_READ); }
    bool isSupplementaryAlignment() const { return isStatusFlagSet(SUPPLEMENTARY_ALIGNMENT); }
    bool isReversedSeq() const { return isStatusFlagSet(READ_REVERSE_STRAND); }

    // set functions
    void set_ReadSeq_reverseSeq();

private:
    AlignLine_Data_t data;
    std::vector<std::pair<char, std::string>> cigar_list;
    bool flag_isReadReversed;

    // functions
    std::vector<uint> getOpListInCIGAR(char Op);
};
