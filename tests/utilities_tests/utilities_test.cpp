#include <gtest/gtest.h>
#include "utility.hpp"

using std::string;
using std::vector;

string sample_record = "M_E200023363L1C007R00600374775	0	1	1447869	37	21M1I10M	*	0	0	CCGGGCGGCGCGGGCGGGCGGCGCGGCACGGG	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:2	XO:i:1	XG:i:1	MD:Z:1G23G5	RG:Z:GXBS2-01";
string sample_read = "GCCACTCCACGCAGCATCGTGACCTGCTGGT";

TEST(utility,revcomp)
{

    string result = revcomp(sample_read);
    string expect = "ACCAGCAGGTCACGATGCTGCGTGGAGTGGC";
    EXPECT_EQ(result, expect);
}

TEST(utility,split)
{
    vector<string> result = split(sample_record);
    vector<string> expect = {"M_E200023363L1C007R00600374775", "0", "1", "1447869", "37", "21M1I10M", "*", "0", "0", "CCGGGCGGCGCGGGCGGGCGGCGCGGCACGGG", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "XT:A:U", "NM:i:3", "X0:i:1", "X1:i:0", "XM:i:2", "XO:i:1", "XG:i:1", "MD:Z:1G23G5", "RG:Z:GXBS2-01"};
    EXPECT_EQ(result, expect);
}

string sample_digitString = "12345678910";
string sample_negativeDigitString = "-12345678910";
string sample_numberWithPoint = "1.01234";
string sample_negativeNumberWithPoint = "-1.12345";
string sample_alphabetString = "abcdefg";
string sample_upperAlphabetString = "HIJKLMN";
string sample_numberAlphaMixed = "12345ABCDEFG";

TEST(utility,isStringDigit)
{
    EXPECT_TRUE(isStringDigit(sample_digitString));
    EXPECT_TRUE(isStringDigit(sample_negativeDigitString));
    EXPECT_TRUE(isStringDigit(sample_negativeNumberWithPoint));
    EXPECT_TRUE(isStringDigit(sample_numberWithPoint));
    EXPECT_FALSE(isStringDigit(sample_numberAlphaMixed));
    EXPECT_FALSE(isStringDigit(sample_alphabetString));
    EXPECT_FALSE(isStringDigit(sample_upperAlphabetString));
}

TEST(utility,isStringAlphabet)
{
    EXPECT_TRUE(isStringAlphabet(sample_alphabetString));
    EXPECT_TRUE(isStringAlphabet(sample_upperAlphabetString));
    EXPECT_FALSE(isStringAlphabet(sample_numberAlphaMixed));
    EXPECT_FALSE(isStringAlphabet(sample_digitString));
    EXPECT_FALSE(isStringAlphabet(sample_negativeDigitString));
    EXPECT_FALSE(isStringAlphabet(sample_negativeNumberWithPoint));
    EXPECT_FALSE(isStringAlphabet(sample_numberWithPoint));
}

string example_tab_in_start_end = "\t123456\t";
string example_have_multiple_tab_in_start_end = "\t\t123456\t\t";
string example_tab_in_start = "\t123456";
string example_tab_in_end = "123456\t";
string example_tab_between_sentence = "123\t456";

TEST(utility,strip)
{
    EXPECT_EQ(strip(example_tab_in_start_end, '\t'), "123456");
    EXPECT_EQ(strip(example_tab_in_start, '\t'), "123456");
    EXPECT_EQ(strip(example_tab_in_end, '\t'), "123456");
    EXPECT_EQ(strip(example_tab_between_sentence, '\t'),"123\t456");
    EXPECT_EQ(strip(example_have_multiple_tab_in_start_end, '\t'), "123456");
}

TEST(utility,lstrip)
{
    EXPECT_EQ(lstrip(example_tab_in_start_end, '\t'), "123456\t");
    EXPECT_EQ(lstrip(example_tab_in_start, '\t'), "123456");
    EXPECT_EQ(lstrip(example_tab_in_end, '\t'), "123456\t");
    EXPECT_EQ(lstrip(example_tab_between_sentence, '\t'),"123\t456");
    EXPECT_EQ(lstrip(example_have_multiple_tab_in_start_end, '\t'), "123456\t\t");
}

TEST(utility,rstrip)
{
    EXPECT_EQ(rstrip(example_tab_in_start_end, '\t'), "\t123456");
    EXPECT_EQ(rstrip(example_tab_in_start, '\t'), "\t123456");
    EXPECT_EQ(rstrip(example_tab_in_end, '\t'), "123456");
    EXPECT_EQ(rstrip(example_tab_between_sentence, '\t'),"123\t456");
    EXPECT_EQ(rstrip(example_have_multiple_tab_in_start_end, '\t'), "\t\t123456");
}