#include <gtest/gtest.h>
#include <stdexcept>
#include "calPMD.hpp"
#include "arguments.hpp"

class DeaminationFunction_test : public ::testing::Test
{
protected:
    void SetUp() override
    {
        FLAGS_CpG = false;
        FLAGS_platypus = false;
        FLAGS_deamination = true;
        FLAGS_range = 2;
        FLAGS_requirebaseq = 0;

        FLAGS_maskterminalbases = 0;
        FLAGS_maskterminaldeaminations = 0;
        IS_USED_maskterminalbases = false;
        IS_USED_maskterminaldeaminations = false;
    }
};

TEST_F(DeaminationFunction_test, Counts_CT_At_Five_Prime_Distance_Zero)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TAAA", "CAAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
}

TEST_F(DeaminationFunction_test, Counts_Forward_Bases_Only_Within_Configured_Range)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"CAGT", "CCCC"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 0);
    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('G')], 0);
    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('C')], 1);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('T')], 0);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('G')], 0);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('C')], 0);

    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 1);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('T')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('G')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('C')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('T')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('G')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('C')], 0);

    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(DeaminationFunction_test, Counts_Reverse_Bases_By_Three_Prime_Distance)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"CAGT", "GGGG"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('T')], 1);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('G')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('C')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('T')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('G')], 1);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('C')], 0);

    EXPECT_EQ(deam_statics.reverse.totals[0], 1);
    EXPECT_EQ(deam_statics.reverse.totals[1], 1);

    EXPECT_EQ(deam_statics.forward.totals[0],0);
    EXPECT_EQ(deam_statics.forward.totals[1],0);
}

TEST_F(DeaminationFunction_test, Skips_Reference_Bases_Other_Than_C_And_G)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"CAGT", "TTTT"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(DeaminationFunction_test, Skips_Ambiguous_Read_Base_At_Five_Prime_End)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"NNNN", "TTTT"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(DeaminationFunction_test, Counts_Forward_Base_At_Last_Position_Inside_Range)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TACG", "TCTT"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 1);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(DeaminationFunction_test,deam_test_07)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TACG", "TTCT"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(DeaminationFunction_test,deam_test_08)
{
    FLAGS_range = 30;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"AAAA", "CCCC"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals.size(), 30);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 1);
    EXPECT_EQ(deam_statics.forward.totals[2], 1);
    EXPECT_EQ(deam_statics.forward.totals[3], 1);
    EXPECT_EQ(deam_statics.forward.totals[4], 0);
    EXPECT_EQ(deam_statics.forward.totals[29], 0);
    
    EXPECT_EQ(deam_statics.reverse.totals.size(), 30);
}

TEST_F(DeaminationFunction_test, Counts_Both_Directions_When_Terminal_Ranges_Overlap)
{
    FLAGS_range = 10;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"ATATA", "CGCGC"};
    std::vector<double> modern_model_deam(5, 0.01);
    std::vector<double> ancient_model_deam(5, 0.01);
    std::string quals = "IIIII";
    std::string maskedseq = "AAAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals.size(), 10);
    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.totals[2], 1);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);
    EXPECT_EQ(deam_statics.forward.counts[4][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.totals[4], 1);
    EXPECT_EQ(deam_statics.reverse.totals.size(), 10);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.counts[1][reverse_index('T')],1);
    EXPECT_EQ(deam_statics.reverse.totals[1], 1);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.counts[3][reverse_index('T')],1);
    EXPECT_EQ(deam_statics.reverse.totals[3], 1);
    EXPECT_EQ(deam_statics.reverse.totals[4], 0);
}

TEST_F(DeaminationFunction_test, Accumulates_Counts_From_Multiple_Reads)
{
    FLAGS_range = 10;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"ATATA", "CGCGC"};
    std::vector<double> modern_model_deam(5, 0.01);
    std::vector<double> ancient_model_deam(5, 0.01);
    std::string quals = "IIIII";
    std::string maskedseq = "AAAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    real_data_t real_data2{"ATATA", "CGCGC"};
    calPMD pmd2(std::move(real_data2), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.totals.size(), 10);
    EXPECT_EQ(deam_statics.forward.totals[0], 2);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 2);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);
    EXPECT_EQ(deam_statics.forward.totals[4], 2);
    EXPECT_EQ(deam_statics.reverse.totals.size(), 10);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 2);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 2);
    EXPECT_EQ(deam_statics.reverse.totals[4], 0);
}