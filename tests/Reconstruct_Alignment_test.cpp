#include <gtest/gtest.h>
#include "seqProcedures.hpp"
#include "parsedData.hpp"

TEST(ReconstructAlignmentTest, SimpleMatchWithoutInsertions)
{
    recordLine_struct_t raw_data{"read1", "0", "chr1", 1, 60, "10M", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "10"}}}};
    parsedData parsed(raw_data);
    alignnmentData_t result{};

    EXPECT_EQ(0, ReconstructAlignmentAndRefSeq(parsed, result));
    EXPECT_EQ("..........", result.alignment);
    EXPECT_EQ("ACGTACGTAC", result.ref_seq);
}

TEST(ReconstructAlignmentTest, MatchWithInsertionAddsGapToRefSeq)
{
    recordLine_struct_t raw_data{"read2", "0", "chr1", 1, 60, "5M1I4M", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "10"}}}};
    parsedData parsed(raw_data);
    alignnmentData_t result{};

    EXPECT_EQ(0, ReconstructAlignmentAndRefSeq(parsed, result));
    EXPECT_EQ("..........", result.alignment);
    EXPECT_EQ("ACGTA-GTAC", result.ref_seq);
}

TEST(ReconstructAlignmentTest, ReversedReadProducesReverseComplementRefSeq)
{
    recordLine_struct_t raw_data{"read3", "16", "chr1", 1, 60, "10M", "ACGTACGTAC", "IIIIIIIIII", {{"MD", {'Z', "10"}}}};
    parsedData parsed(raw_data);
    alignnmentData_t result{};

    EXPECT_EQ(0, ReconstructAlignmentAndRefSeq(parsed, result));
    EXPECT_EQ("..........", result.alignment);
    EXPECT_EQ("GTACGTACGT", result.ref_seq);
}

TEST(ReconstructAlignmentTest, MissingMDTagReturnsError)
{
    recordLine_struct_t raw_data{"read4", "0", "chr1", 1, 60, "10M", "ACGTACGTAC", "IIIIIIIIII", {}};
    parsedData parsed(raw_data);
    alignnmentData_t result{};

    EXPECT_EQ(-1, ReconstructAlignmentAndRefSeq(parsed, result));
    EXPECT_TRUE(result.alignment.empty());
    EXPECT_TRUE(result.ref_seq.empty());
}

