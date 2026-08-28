#include <gtest/gtest.h>
#include "calPMD.hpp"
#include "arguments.hpp"

class Deamination_noCpG_test:public ::testing::Test
{
    protected:
    void SetUp() override
    {
        FLAGS_CpG = false;
        FLAGS_noCpG = true;
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

TEST_F(Deamination_noCpG_test,one_pair_of_CpG_in_the_end)
{
    FLAGS_range = 4;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{
        .real_read = "TAAA",
        .real_ref_seq = "AACG"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 0);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 0);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('A')], 0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
}

TEST_F(Deamination_noCpG_test, noCpG_reference_only_have_C)
{
    FLAGS_range = 4;
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{
        .real_read = "TACC",
        .real_ref_seq = "CCCC"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 1);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('C')], 1);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('C')], 0);
    EXPECT_EQ(deam_statics.forward.totals[0], 1);
    EXPECT_EQ(deam_statics.forward.totals[1], 1);
    EXPECT_EQ(deam_statics.forward.totals[2], 1);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);

    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
}

TEST_F(Deamination_noCpG_test, all_C_and_G_sites_are_CpG_or_terminal)
{
    FLAGS_range = 4;
    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{
        .real_read = "TACC",
        .real_ref_seq = "GCGC"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[0][forward_index('T')], 0);
    EXPECT_EQ(deam_statics.forward.counts[1][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('C')], 0);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('C')], 0);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 0);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);

    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
}

TEST_F(Deamination_noCpG_test,CpG_pair_at_begin_and_end)
{
    FLAGS_range = 6;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{
        .real_read = "TGAATG",
        .real_ref_seq = "CGCCCG"};
    std::vector<double> modern_model_deam(6, 0.01);
    std::vector<double> ancient_model_deam(6, 0.01);
    std::string quals = "IIIIII";
    std::string maskedseq = "AAAAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('A')], 1);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 1);
    EXPECT_EQ(deam_statics.forward.totals[3], 1);
    EXPECT_EQ(deam_statics.forward.totals[4], 0);
    EXPECT_EQ(deam_statics.forward.totals[5], 0);

    EXPECT_EQ(deam_statics.reverse.counts[0][reverse_index('G')],0);
    EXPECT_EQ(deam_statics.reverse.counts[4][reverse_index('G')],0);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 0);
    EXPECT_EQ(deam_statics.reverse.totals[3], 0);
    EXPECT_EQ(deam_statics.reverse.totals[4], 0);
    EXPECT_EQ(deam_statics.reverse.totals[5], 0);
}

TEST_F(Deamination_noCpG_test,CpG_pair_at_begin_and_end_fill_G)
{
    FLAGS_range = 6;

    platypus_statics_dicts_t platypus_statics_dict{};
    real_data_t real_data{
        .real_read = "TGAATG",
        .real_ref_seq = "CGGGCG"};
    std::vector<double> modern_model_deam(6, 0.01);
    std::vector<double> ancient_model_deam(6, 0.01);
    std::string quals = "IIIIII";
    std::string maskedseq = "AAAAAA";

    platypus_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    deamination_statics_t deam_statics(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, platypus_statics_dict, denom, deam_statics);

    EXPECT_EQ(deam_statics.forward.counts[2][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.counts[3][forward_index('A')], 0);
    EXPECT_EQ(deam_statics.forward.totals[0], 0);
    EXPECT_EQ(deam_statics.forward.totals[1], 0);
    EXPECT_EQ(deam_statics.forward.totals[2], 0);
    EXPECT_EQ(deam_statics.forward.totals[3], 0);
    EXPECT_EQ(deam_statics.forward.totals[4], 0);
    EXPECT_EQ(deam_statics.forward.totals[5], 0);

    EXPECT_EQ(deam_statics.reverse.counts[2][reverse_index('A')],1);
    EXPECT_EQ(deam_statics.reverse.counts[3][reverse_index('A')],1);
    EXPECT_EQ(deam_statics.reverse.totals[0], 0);
    EXPECT_EQ(deam_statics.reverse.totals[1], 0);
    EXPECT_EQ(deam_statics.reverse.totals[2], 1);
    EXPECT_EQ(deam_statics.reverse.totals[3], 1);
    EXPECT_EQ(deam_statics.reverse.totals[4], 0);
    EXPECT_EQ(deam_statics.reverse.totals[5], 0);
}