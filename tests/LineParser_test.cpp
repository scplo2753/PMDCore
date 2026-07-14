#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include "parsedData.hpp"
#include "seqProcedures.hpp"

#ifndef PMDCORE_TEST_DATA_DIR
#define PMDCORE_TEST_DATA_DIR ""
#endif

using std::string;

string legalLine{"M_E200023363L1C041R03604840354	0	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};

TEST(valid_line, emptyLine)
{
    string line{
        ""};
    AlignLine_Data_t data;
    EXPECT_FALSE(validAndParse(line,data));
}

TEST(valid_line,headerLine)
{
    AlignLine_Data_t data;
    EXPECT_FALSE(validAndParse("@HD VN:1.6  SO:coordinate",data));
    EXPECT_FALSE(validAndParse("@SQ SN:4    LN:18585056",data));
    EXPECT_FALSE(validAndParse("@PG ID:samtools PN:samtools VN:1.21 CL:samtools view -h temp.bam", data));
}

TEST(valid_line,parseLowQualLine)
{
    AlignLine_Data_t data;
    string lowQualityScore{"ancient_read_001	0	chr1	1000	1	1M	*	0	0	C	!"}; //test low quality score
    EXPECT_FALSE(validAndParse(lowQualityScore,data));
    EXPECT_EQ(splitOneLine(lowQualityScore,data),0);
}

TEST(valid_line,legalLine)
{
    AlignLine_Data_t data;
    EXPECT_TRUE(validAndParse(legalLine,data));
    EXPECT_EQ(splitOneLine(legalLine,data),0);
}

TEST(valid_line,lackOfSessionLine)
{
    AlignLine_Data_t data;
    string lackOfSession{"read1\t0\tchr1"};
    EXPECT_FALSE(validAndParse(lackOfSession,data));
    EXPECT_EQ(splitOneLine(lackOfSession,data),BAD_RECORD_HAVE_EMPTY_FIELD);
}

TEST(illegal_line_filter_test,emptyLine)
{
    AlignLine_Data_t data;
    string empty{""};
    EXPECT_FALSE(validAndParse(empty,data));
}

TEST(illegal_line_filter_test,illegalFlag)
{
    AlignLine_Data_t data;
    string alpha_in_flag{"M_E200023363L1C041R03604840354	illegal_flag	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_flag{"M_E200023363L1C041R03604840354	-16	16	81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_flag,data),BAD_RECORD_ILLEGAL_FLAG);
    EXPECT_FALSE(validAndParse(alpha_in_flag,data));
    EXPECT_EQ(splitOneLine(negative_in_flag,data),BAD_RECORD_ILLEGAL_FLAG);
    EXPECT_FALSE(validAndParse(negative_in_flag,data));
}

TEST(illegal_line_filter_test,illegalPos)
{
    AlignLine_Data_t data;
    string alpha_in_pos{"M_E200023363L1C041R03604840354	0	16	alpha	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_pos{"M_E200023363L1C041R03604840354	0	16	-81001900	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string pos_is_empty{"M_E200023363L1C041R03604840354	0	16		37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_pos, data), BAD_RECORD_ILLEGAL_POS);
    EXPECT_FALSE(validAndParse(alpha_in_pos, data));
    EXPECT_EQ(splitOneLine(negative_in_pos, data), BAD_RECORD_ILLEGAL_POS);
    EXPECT_FALSE(validAndParse(negative_in_pos, data));
    EXPECT_EQ(splitOneLine(pos_is_empty, data), BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(pos_is_empty, data));
}

TEST(illegal_line_filter_test,emptyRecord)
{
    AlignLine_Data_t data;
    string emptyRecord{"\t\t\t\t\t\t\t\t\t\t"};
    EXPECT_EQ(splitOneLine(emptyRecord,data),BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(emptyRecord,data));
}

TEST(illegal_line_filter_test,illegalMapQ)
{
    AlignLine_Data_t data;
    string alpha_in_MAPQ{"example2	32	1	705745	alpha	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string negative_in_MAPQ{"example2	32	1	705745	-37	37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    string MAPQ_is_empty{"example2	32	1	705745			37	49M	*	0	0	CTGCTTAAAAATTAAATGCAGAAAAAAGAGTAAGCTAGAAGTGATGACA	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:49	RG:Z:GXBS2-01"};
    EXPECT_EQ(splitOneLine(alpha_in_MAPQ,data),BAD_RECORD_ILLEGAL_MAPQ);
    EXPECT_FALSE(validAndParse(alpha_in_MAPQ,data));
    EXPECT_EQ(splitOneLine(negative_in_MAPQ,data),BAD_RECORD_ILLEGAL_MAPQ);
    EXPECT_FALSE(validAndParse(negative_in_MAPQ,data));
    EXPECT_EQ(splitOneLine(MAPQ_is_empty,data),BAD_RECORD_HAVE_EMPTY_FIELD);
    EXPECT_FALSE(validAndParse(MAPQ_is_empty,data));
}

