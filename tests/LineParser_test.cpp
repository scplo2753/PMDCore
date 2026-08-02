#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "parsedData.hpp"
#include "seqProcedures.hpp"

#ifndef PMDCORE_TEST_DATA_DIR
#define PMDCORE_TEST_DATA_DIR ""
#endif

using std::string;
using std::vector;
using splited_record_t = vector<string>;

string legalLine{"M_E200023363L1C041R03604840354	0	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};

TEST(valid_line, emptyLine)
{
    string line{
        ""};
    recordLine_struct_t data;
    splited_record_t splited_record;
    EXPECT_FALSE(validAndParse(line, splited_record, data));
}

TEST(valid_line,headerLine)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    EXPECT_FALSE(validAndParse("@HD VN:1.6  SO:coordinate", splited_record, data));
    EXPECT_FALSE(validAndParse("@SQ SN:4    LN:18585056", splited_record, data));
    EXPECT_FALSE(validAndParse("@PG ID:samtools PN:samtools VN:1.21 CL:samtools view -h temp.bam", splited_record, data));
}

TEST(valid_line,parseLowQualLine)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string lowQualityScore{"ancient_read_001	0	chr1	1000	1	1M	*	0	0	C	!"}; //test low quality score
    EXPECT_FALSE(validAndParse(lowQualityScore, splited_record, data));
    EXPECT_EQ(splitOneLine(lowQualityScore, splited_record, data), 0);
}

TEST(valid_line,legalLine)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    EXPECT_TRUE(validAndParse(legalLine, splited_record, data));
    EXPECT_EQ(splitOneLine(legalLine,splited_record,data),0);
}

TEST(valid_line,lackOfSessionLine)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string lackOfSession{"read1\t0\tchr1"};
    EXPECT_FALSE(validAndParse(lackOfSession, splited_record, data));
    EXPECT_EQ(splitOneLine(lackOfSession, splited_record, data), BAD_RECORD_HAVE_EMPTY_FIELD);
}

TEST(illegal_line_filter_test,illegalFlag)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string alpha_in_flag{"M_E200023363L1C041R03604840354	illegal_flag	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_flag{"M_E200023363L1C041R03604840354	-16	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_flag,splited_record,data),BAD_RECORD_ILLEGAL_FLAG);
    EXPECT_FALSE(validAndParse(alpha_in_flag,splited_record,data));
    EXPECT_EQ(splitOneLine(negative_in_flag,splited_record,data),BAD_RECORD_ILLEGAL_FLAG);
    EXPECT_FALSE(validAndParse(negative_in_flag,splited_record,data));
}

TEST(illegal_line_filter_test,illegalPos)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string alpha_in_pos{"M_E200023363L1C041R03604840354	0	16	alpha	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_pos{"M_E200023363L1C041R03604840354	0	16	-81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string pos_is_empty{"M_E200023363L1C041R03604840354	0	16		37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_pos, splited_record,data), BAD_RECORD_ILLEGAL_POS);
    EXPECT_FALSE(validAndParse(alpha_in_pos,splited_record, data));
    EXPECT_EQ(splitOneLine(negative_in_pos,splited_record, data), BAD_RECORD_ILLEGAL_POS);
    EXPECT_FALSE(validAndParse(negative_in_pos,splited_record, data));
    EXPECT_EQ(splitOneLine(pos_is_empty,splited_record, data), BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(pos_is_empty,splited_record, data));
}

TEST(illegal_line_filter_test,emptyRecord)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string emptyRecord{"\t\t\t\t\t\t\t\t\t\t"};
    EXPECT_EQ(splitOneLine(emptyRecord,splited_record,data),BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(emptyRecord,splited_record,data));
}

TEST(illegal_line_filter_test,illegalMapQ)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string alpha_in_MAPQ{"example2	32	1	705745	alpha	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_MAPQ{"example2	32	1	705745	-37	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string MAPQ_is_empty{"example2	32	1	705745			37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_MAPQ,splited_record,data),BAD_RECORD_ILLEGAL_MAPQ);
    EXPECT_FALSE(validAndParse(alpha_in_MAPQ,splited_record,data));
    EXPECT_EQ(splitOneLine(negative_in_MAPQ,splited_record,data),BAD_RECORD_ILLEGAL_MAPQ);
    EXPECT_FALSE(validAndParse(negative_in_MAPQ,splited_record,data));
    EXPECT_EQ(splitOneLine(MAPQ_is_empty,splited_record,data),BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(MAPQ_is_empty,splited_record,data));
}

TEST(illegal_line_filter_test,missingReadOrQual)
{
    recordLine_struct_t data;
    splited_record_t splited_record;
    string missingRead{"read1	0	chr1	1000	37	10M	*	0	0		IIIIIIIIII	MD:Z:10"};
    string missingQual{"read1	0	chr1	1000	37	10M	*	0	0	ACGTACGTAC		MD:Z:10"};

    EXPECT_EQ(splitOneLine(missingRead, splited_record, data), BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(missingRead, splited_record, data));

    EXPECT_EQ(splitOneLine(missingQual, splited_record, data), BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(missingQual, splited_record, data));
}

