#include <gtest/gtest.h>
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