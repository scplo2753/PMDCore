#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "sam/enum_flags.hpp"
#include "sam/struct_record.hpp"

enum class parsedRecordError
{
    SUCCESS,
    CIGAR_EMPTY,
    BAD_CIGAR_FORMAT,
    INVALID_CIGAR_OPERATION,
    INVALID_CIGAR_STEP,
    MD_TAG_NOT_EXIST,
    MD_TAG_EMPTY
};

class parsedData;

using parsedRecordResult = std::variant<parsedData, parsedRecordError>;
using CIGARList_t = std::vector<std::pair<char, std::size_t>>;

class parsedData
{
   public:
    // entry function
    [[nodiscard]]
    static parsedRecordResult parseRawData(const recordLine_struct_t& record);

    // get raw data functions
    std::string getQNAME() const
    {
        return data.QNAME;
    }
    std::string getRNAME() const
    {
        return data.RNAME;
    }
    std::string getRawCigar() const
    {
        return data.cigar;
    }
    int getRawFLAG() const
    {
        return std::stoi(data.FLAG);
    }
    std::string getReadSeq() const
    {
        return data.read_seq;
    }
    int getMAPQ() const
    {
        return data.MAPQ;
    }
    int getPosition() const
    {
        return data.POS;
    }
    std::string getQualityScores() const
    {
        return data.quality_scores;
    }
    const options_map_t& getOptionsMap() const
    {
        return data.options_map;
    }

    // tag functions
    bool isTagExists(const std::string& tag) const;
    std::string getTagValue(const std::string& tag) const;

    // CIGAR functions
    CIGARList_t getCIGARList() const
    {
        return cigar_list;
    }
    std::vector<std::size_t> getInsertionList()
    {
        return getOpListInCIGAR('I');
    }
    std::vector<std::size_t> getSoftClipList()
    {
        return getOpListInCIGAR('S');
    }
    std::vector<std::size_t> getMatchList()
    {
        return getOpListInCIGAR('M');
    }
    std::vector<std::size_t> getDeletionList()
    {
        return getOpListInCIGAR('D');
    }

    // MD functions
    std::string getMDTagValue() const
    {
        return getTagValue("MD");
    }
    std::vector<std::string> getMDList() const;

    // status flag functions
    bool isStatusFlagSet(SamFlag flag) const;
    bool isReadPaired() const
    {
        return isStatusFlagSet(SamFlag::READ_PAIRED);
    }
    bool isProperPair() const
    {
        return isStatusFlagSet(SamFlag::PROPER_PAIR);
    }
    bool isReadUnmapped() const
    {
        return isStatusFlagSet(SamFlag::READ_UNMAPPED);
    }
    bool isMateUnmapped() const
    {
        return isStatusFlagSet(SamFlag::MATE_UNMAPPED);
    }
    bool isReadReverseStrand() const
    {
        return isStatusFlagSet(SamFlag::READ_REVERSE_STRAND);
    }
    bool isMateReverseStrand() const
    {
        return isStatusFlagSet(SamFlag::MATE_REVERSE_STRAND);
    }
    bool isFirstInPair() const
    {
        return isStatusFlagSet(SamFlag::FIRST_IN_PAIR);
    }
    bool isSecondInPair() const
    {
        return isStatusFlagSet(SamFlag::SECOND_IN_PAIR);
    }
    bool isNotPrimaryAlignment() const
    {
        return isStatusFlagSet(SamFlag::NOT_PRIMARY_ALIGNMENT);
    }
    bool isReadFailsQualityChecks() const
    {
        return isStatusFlagSet(SamFlag::READ_FAILS_QUALITY_CHECKS);
    }
    bool isDuplicateRead() const
    {
        return isStatusFlagSet(SamFlag::DUPLICATE_READ);
    }
    bool isSupplementaryAlignment() const
    {
        return isStatusFlagSet(SamFlag::SUPPLEMENTARY_ALIGNMENT);
    }
    bool isReversedSeq() const
    {
        return isStatusFlagSet(SamFlag::READ_REVERSE_STRAND);
    }

    // set functions
    void set_ReadSeq_reverseSeq();

   private:
    // construct function
    parsedData(recordLine_struct_t data, CIGARList_t cigarList)
        : data(std::move(data)), cigar_list(std::move(cigarList)), flag_isReadReversed(false) {};

    // functions
    std::vector<std::size_t> getOpListInCIGAR(char Op);
    [[nodiscard]] static parsedRecordError parseCIGAR(const recordLine_struct_t& raw_data, CIGARList_t& output);

    // factoy functions
    static bool factory_isTagNotExists(const recordLine_struct_t& data, const std::string& tag);
    static bool factory_isTagValueEmpty(const recordLine_struct_t& raw_data, const std::string& Tag);

    recordLine_struct_t data;
    CIGARList_t cigar_list;
    bool flag_isReadReversed;
};