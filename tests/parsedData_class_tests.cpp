#include <gtest/gtest.h>
#include "utility.hpp"
#include "seqProcedures.hpp"
#include "parsedData.hpp"
#include "arguments.hpp"

using std::string;
using std::string_view;

AlignLine_Data_t standardLine{"read1", "0", "chr1", 100, 60, 
    "3M1D4X3M", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "10"}}}};
AlignLine_Data_t dataForFlagTest{"read2", "4095", "chr2", 200, 30, 
    "5M5S", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "5"}}}};

parsedData standardRecord(standardLine);
parsedData flagTestData(dataForFlagTest);

TEST(ParsedDataClassTest, ConstructorTest) {
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
    EXPECT_TRUE(standardRecord.isTagExists("MD"));
    EXPECT_EQ(standardRecord.getTagValue("MD"), "10");
    EXPECT_FALSE(standardRecord.isTagExists("NM"));
    EXPECT_EQ(standardRecord.getTagValue("NM"), "");
}

TEST(ParsedDataClassTest,test_of_isStatusFlagSet_with_false)
{
    EXPECT_FALSE(standardRecord.isStatusFlagSet(0x4));
}

TEST(ParsedDataClassTest,test_of_isStatusFlagset_with_true)
{
    EXPECT_TRUE(flagTestData.isStatusFlagSet(0x4));
}

TEST(ParsedDataClassTest, status_FLAG_functions_test_FALSE){
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
    EXPECT_EQ(standardRecord.getMDList().at(0),"10");
}

TEST(ParsedDataClassTest,test_of_getMDList_with_multiplePos)
{
    parsedData multiple_pos_in_MD_Field(AlignLine_Data_t{"MultiMD","0","chr2", 200, 30, 
    "5M5S", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "15^A12BC"}}}});
    const std::vector<string> expected{"15", "^A", "12", "B", "C"};
    EXPECT_EQ(multiple_pos_in_MD_Field.getMDList().size(), 5);
    EXPECT_EQ(multiple_pos_in_MD_Field.getMDList(), expected);
}

TEST(PasedDataClassTest,test_of_getCIGARList_with_standard_record)
{
    std::vector<std::pair<char,string>> expected{{'M', "3"},{'D', "1"},{'X', "4"}, {'M', "3"}};
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
    std::vector<uint> expected_softclip{0, 5};
    std::vector<uint> expected_match{0, 1, 2, 8, 9, 10};
    std::vector<uint> expected_deletion{3};

    EXPECT_EQ(standardRecord.getInsertionList(), std::vector<uint>{});
    EXPECT_EQ(standardRecord.getMatchList(), expected_match);
    EXPECT_EQ(standardRecord.getDeletionList(), expected_deletion);
}

TEST(ParsedDataClassTest,set_ReadSeq_reverseSeq)
{
    std::vector<string> splited_record;
    AlignLine_Data_t alignLine_data;
    string raw_record = "M_E200023363L1C022R03201788808	0	1	1363297	37	4M1D26M	*	0	0	ACGACCCGGGCTCCCCCGCCCGCAGCGCCG	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:2	XO:i:1	XG:i:1	MD:Z:4^C5C16C3	RG:Z:GXBS2-01";
    validAndParse(raw_record,splited_record,alignLine_data);
    parsedData parsed_data(alignLine_data);
    parsed_data.set_ReadSeq_reverseSeq();
    EXPECT_EQ(parsed_data.getReadSeq(), "CGGCGCTGCGGGCGGGGGAGCCCGGGTCGT");
}