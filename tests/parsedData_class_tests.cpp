#include <gtest/gtest.h>
#include <limits>
#include "parsedData_test_helpers.hpp"
#include "utilities/sequence_utils.hpp"
#include "seqProcedures.hpp"
#include "parsedData.hpp"
#include "arguments.hpp"

using std::string;
using std::string_view;

recordLine_struct_t standardLine{"read1", "0", "chr1", 100, 60,
    "3M1D4X3M", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "10"}}}};
recordLine_struct_t dataForFlagTest{"read2", "4095", "chr2", 200, 30,
    "5M5S", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "5"}}}};



TEST(ParsedDataClassTest, ConstructorTest) {
    auto standardRecord = requireParsedRecord(standardLine);
    EXPECT_EQ(standardRecord.getQNAME(), "read1");
    EXPECT_EQ(standardRecord.getRawFLAG(), 0);
    EXPECT_EQ(standardRecord.getRNAME(), "chr1");
    EXPECT_EQ(standardRecord.getPosition(), 100);
    EXPECT_EQ(standardRecord.getMAPQ(), 60);
    EXPECT_EQ(standardRecord.getRawCigar(), "3M1D4X3M");
    EXPECT_EQ(standardRecord.getReadSeq(), "ACGTACGTAC");
    EXPECT_EQ(standardRecord.getQualityScores(), "IIIIIIIIII");
    EXPECT_EQ(standardRecord.getOptionsMap().size(), 1);
    EXPECT_EQ(standardRecord.getOptionsMap().at("MD").value, "10");
}

TEST(ParsedDataClassTest,tag_functions_test){
    auto standardRecord = requireParsedRecord(standardLine);
    EXPECT_TRUE(standardRecord.isTagExists("MD"));
    EXPECT_EQ(standardRecord.getTagValue("MD"), "10");
    EXPECT_FALSE(standardRecord.isTagExists("NM"));
    EXPECT_EQ(standardRecord.getTagValue("NM"), "");
}

TEST(ParsedDataClassTest,test_of_isStatusFlagSet_with_false)
{
    auto standardRecord = requireParsedRecord(standardLine);
    EXPECT_FALSE(standardRecord.isStatusFlagSet(SamFlag::READ_UNMAPPED));
}

TEST(ParsedDataClassTest,test_of_isStatusFlagset_with_true)
{
    auto flagTestData = requireParsedRecord(dataForFlagTest);
    EXPECT_TRUE(flagTestData.isStatusFlagSet(SamFlag::READ_UNMAPPED));
}

TEST(ParsedDataClassTest, status_FLAG_functions_test_FALSE){
    auto standardRecord = requireParsedRecord(standardLine);
    EXPECT_FALSE(standardRecord.isReadPaired());
    EXPECT_FALSE(standardRecord.isProperPair());
    EXPECT_FALSE(standardRecord.isReadUnmapped());
    EXPECT_FALSE(standardRecord.isMateUnmapped());
    EXPECT_FALSE(standardRecord.isReadReverseStrand());
    EXPECT_FALSE(standardRecord.isMateReverseStrand());
    EXPECT_FALSE(standardRecord.isFirstInPair());
    EXPECT_FALSE(standardRecord.isSecondInPair());
    EXPECT_FALSE(standardRecord.isNotPrimaryAlignment());
    EXPECT_FALSE(standardRecord.isReadFailsQualityChecks());
    EXPECT_FALSE(standardRecord.isDuplicateRead());
    EXPECT_FALSE(standardRecord.isSupplementaryAlignment());
}

TEST(ParsedDataClassTest,status_FLAG_functions_test_TRUE){
    auto flagTestData = requireParsedRecord(dataForFlagTest);
    EXPECT_TRUE(flagTestData.isReadPaired());
    EXPECT_TRUE(flagTestData.isProperPair());
    EXPECT_TRUE(flagTestData.isReadUnmapped());
    EXPECT_TRUE(flagTestData.isMateUnmapped());
    EXPECT_TRUE(flagTestData.isReadReverseStrand());
    EXPECT_TRUE(flagTestData.isMateReverseStrand());
    EXPECT_TRUE(flagTestData.isFirstInPair());
    EXPECT_TRUE(flagTestData.isSecondInPair());
    EXPECT_TRUE(flagTestData.isNotPrimaryAlignment());
    EXPECT_TRUE(flagTestData.isReadFailsQualityChecks());
    EXPECT_TRUE(flagTestData.isDuplicateRead());
    EXPECT_TRUE(flagTestData.isSupplementaryAlignment());
}

TEST(ParsedDataClassTest,test_of_getMDList)
{
    auto standardRecord = requireParsedRecord(standardLine);
    EXPECT_EQ(standardRecord.getMDList().at(0),"10");
}

TEST(ParsedDataClassTest,test_of_getMDList_with_multiplePos)
{
    auto multiple_pos_in_MD_Field = requireParsedRecord(recordLine_struct_t{"MultiMD","0","chr2", 200, 30,
    "5M5S", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "15^A12BC"}}}});
    const std::vector<string> expected{"15", "^A", "12", "B", "C"};
    EXPECT_EQ(multiple_pos_in_MD_Field.getMDList().size(), 5);
    EXPECT_EQ(multiple_pos_in_MD_Field.getMDList(), expected);
}

TEST(PasedDataClassTest,test_of_getCIGARList_with_standard_record)
{
    auto standardRecord = requireParsedRecord(standardLine);
    CIGARList_t expected{{'M', 3}, {'D', 1}, {'X', 4}, {'M', 3}};
    EXPECT_EQ(standardRecord.getCIGARList(), expected);
}

TEST(ArgumentsFunctionsTest, BasicFilterReturnsTrueForCpGMismatch)
{
    FLAGS_basic = 2;
    FLAGS_CpG = true;
    FLAGS_requirebaseq = 0;

    EXPECT_TRUE(function_basicFilter("CGG", "TGG", 3, "AAA"));
}

TEST(ArgumentsFunctionsTest, BasicFilterReturnsFalseForNonMatchingBases)
{
    FLAGS_basic = 2;
    FLAGS_CpG = true;
    FLAGS_requirebaseq = 0;

    EXPECT_FALSE(function_basicFilter("CCG", "CCG", 3, "!!!"));
}

TEST(ArgumentsFunctionsTest, BasicFilterHandlesEmptyInputSafely)
{
    FLAGS_basic = 2;
    FLAGS_CpG = true;
    FLAGS_requirebaseq = 0;

    EXPECT_FALSE(function_basicFilter("", "", 0, ""));
    EXPECT_FALSE(function_basicTerminal("", "", ""));
}

TEST(ParsedDataClassTest,test_of_getOpListInCIGAR)
{
    auto standardRecord = requireParsedRecord(standardLine);
    std::vector<std::size_t> expected_softclip{0, 5};
    std::vector<std::size_t> expected_match{0, 1, 2, 8, 9, 10};
    std::vector<std::size_t> expected_deletion{3};

    EXPECT_EQ(standardRecord.getInsertionList(), std::vector<std::size_t>{});
    EXPECT_EQ(standardRecord.getMatchList(), expected_match);
    EXPECT_EQ(standardRecord.getDeletionList(), expected_deletion);
}

TEST(ParsedDataClassTest,set_ReadSeq_reverseSeq)
{
    std::vector<string> splited_record;
    recordLine_struct_t alignLine_data;
    string raw_record = "M_E200023363L1C022R03201788808	0	1	1363297	37	4M1D26M	*	0	0	ACGACCCGGGCTCCCCCGCCCGCAGCGCCG	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:2	XO:i:1	XG:i:1	MD:Z:4^C5C16C3	RG:Z:GXBS2-01";
    ASSERT_TRUE(validAndParse(raw_record,splited_record,alignLine_data));
    auto parsed_data = requireParsedRecord(alignLine_data);
    parsed_data.set_ReadSeq_reverseSeq();
    EXPECT_EQ(parsed_data.getReadSeq(), "CGGCGCTGCGGGCGGGGGAGCCCGGGTCGT");
}

TEST(ParsedDataFactoryTest, ReportsCigarErrorsWithoutThrowing)
{
    struct Case { std::string cigar; parsedRecordError expected; };
    const Case cases[] = {
        {"", parsedRecordError::CIGAR_EMPTY},
        {"*", parsedRecordError::CIGAR_EMPTY},
        {"M10M", parsedRecordError::BAD_CIGAR_FORMAT},
        {"10M5", parsedRecordError::BAD_CIGAR_FORMAT},
        {"2M?", parsedRecordError::BAD_CIGAR_FORMAT},
        {"10Q", parsedRecordError::INVALID_CIGAR_OPERATION},
        {"2M1Q", parsedRecordError::INVALID_CIGAR_OPERATION},
        {std::to_string(std::numeric_limits<std::size_t>::max()) + "0M",
         parsedRecordError::INVALID_CIGAR_STEP}
    };
    for (const auto& test : cases)
    {
        SCOPED_TRACE(test.cigar);
        auto raw = standardLine;
        raw.cigar = test.cigar;
        parsedRecordResult result = parsedRecordError::SUCCESS;
        ASSERT_NO_THROW(result = parsedData::parseRawData(raw));
        const auto* error = std::get_if<parsedRecordError>(&result);
        ASSERT_NE(error, nullptr);
        EXPECT_EQ(*error, test.expected);
        EXPECT_EQ(raw.cigar, test.cigar);
    }
}

TEST(ParsedDataFactoryTest, PreservesSourceAndOwnsParsedCigar)
{
    auto raw = standardLine;
    parsedRecordResult result = parsedRecordError::SUCCESS;
    ASSERT_NO_THROW(result = parsedData::parseRawData(raw));
    auto* parsed = std::get_if<parsedData>(&result);
    ASSERT_NE(parsed, nullptr);
    EXPECT_EQ(raw.cigar, standardLine.cigar);
    EXPECT_EQ(raw.read_seq, standardLine.read_seq);
    raw.cigar = "1M";
    raw.read_seq = "T";
    const CIGARList_t expected{{'M', 3}, {'D', 1}, {'X', 4}, {'M', 3}};
    EXPECT_EQ(parsed->getCIGARList(), expected);
    EXPECT_EQ(parsed->getReadSeq(), standardLine.read_seq);
}
