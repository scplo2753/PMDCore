#include <gtest/gtest.h>
#include "calPMD.hpp"
#include "arguments.hpp"

class Deamination_CpG_test:public ::testing::Test
{
    protected:
    void SetUp() override
    {
        FLAGS_CpG = true;
        FLAGS_noCpG = false;
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

TEST_F(Deamination_CpG_test,one_pair_of_CpG_in_the_end)
{
    FLAGS_range = 4;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TAAA", "AACG"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 0);
    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 1);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')], 1);
    EXPECT_EQ(deam_statics.reverse.totals[0], 1);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
}

TEST_F(Deamination_CpG_test, CpG_Uses_Reference_When_Read_Neighbor_Is_Not_G)
{
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TAAA", "CGAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);

    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
}

TEST_F(Deamination_CpG_test, CpG_Skips_Site_When_Only_Read_Neighbor_Is_G)
{
    FLAGS_range = 3;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TGAA", "CAAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')],0);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 0);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')],0);
    EXPECT_EQ(deam_statics.reverse.totals[0],0);
    EXPECT_EQ(deam_statics.reverse.totals[1],0);
    EXPECT_EQ(deam_statics.reverse.totals[2],0);
}

TEST_F(Deamination_CpG_test, CpG_Reverse_Skips_When_Only_Read_Previous_Base_Is_C)
{
    FLAGS_range = 4;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"CAAA", "AGAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(
        static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(
        static_cast<size_t>(FLAGS_range));

    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
}

TEST_F(Deamination_CpG_test, CpG_pair_at_begin_and_end)
{
    FLAGS_range = 6;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{"TGAATG", "CGCCCG"};
    std::vector<double> modern_model_deam(6, 0.01);
    std::vector<double> ancient_model_deam(6, 0.01);
    std::string quals = "IIIIII";
    std::string maskedseq = "AAAAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 1);
    EXPECT_EQ(deam_statics.forward.counts[4][forward_index('T')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 0);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);
    EXPECT_EQ(deam_statics.forward.totals[4], 1);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('G')],1);
    EXPECT_EQ(deam_statics.reverse.counts[4][reverse_index('G')],1);
    EXPECT_EQ(deam_statics.reverse.totals[0], 1);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
    EXPECT_EQ(deam_statics.reverse.totals[4], 1);
}