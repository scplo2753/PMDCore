#include <gtest/gtest.h>
#include "calPMD.hpp"
#include "utility.hpp"
#include "arguments.hpp"

class calPMDTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        FLAGS_CpG = false;
        FLAGS_ss = false;
        FLAGS_platypus = false;
        FLAGS_requirebaseq = 0;
        FLAGS_threshold = -20000.0;
        FLAGS_upperthreshold = 1000000.0;
        FLAGS_maskterminaldeaminations = 0;
        IS_USED_maskterminaldeaminations = false;
        FLAGS_polymorphism_ancient = 0.001;
        FLAGS_polymorphism_contamination = 0.001;
    }
};

TEST_F(calPMDTest, threshold_filter_returns_true_for_PerfectMatch)
{
    statics_dicts_t statics_dict{};
    real_data_t real_data{"AAAA", "AAAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    statics_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, statics_dict, denom);

    EXPECT_TRUE(pmd.threshold_filter());
    EXPECT_TRUE(statics_dict.mismatch_dict.empty());
    EXPECT_TRUE(statics_dict.mismatch_dict_rev.empty());
    EXPECT_TRUE(statics_dict.mismatch_dict_CpG.empty());
    EXPECT_TRUE(statics_dict.mismatch_dict_CpG_rev.empty());
}

TEST_F(calPMDTest, platypus_increments_mismatch_dictionaries_for_CT_Mismatch)
{
    FLAGS_platypus = true;

    statics_dicts_t statics_dict{};
    real_data_t real_data{"TTTT", "CTTT"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "TTTT";

    statics_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, statics_dict, denom);

    EXPECT_EQ(statics_dict.mismatch_dict["CT0"], 1);
    EXPECT_EQ(statics_dict.mismatch_dict_rev["CT3"], 1);
    EXPECT_TRUE(pmd.threshold_filter());
}

TEST_F(calPMDTest, ThresholdFilterRespectsConfiguredBounds)
{
    FLAGS_threshold = 1.0;
    FLAGS_upperthreshold = 2.0;

    statics_dicts_t statics_dict{};
    real_data_t real_data{"AAAA", "AAAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq = "AAAA";

    statics_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, statics_dict, denom);

    EXPECT_FALSE(pmd.threshold_filter());
}

TEST_F(calPMDTest, EmptyMaskedSequenceIsAllowedWhenMaskingIsDisabled)
{
    statics_dicts_t statics_dict{};
    real_data_t real_data{"AAAA", "AAAA"};
    std::vector<double> modern_model_deam(4, 0.01);
    std::vector<double> ancient_model_deam(4, 0.01);
    std::string quals = "IIII";
    std::string maskedseq;

    statics_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    calPMD pmd(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, statics_dict, denom);

    EXPECT_TRUE(pmd.get_maskedSeq().empty());
}

TEST_F(calPMDTest, CpGCheckHandlesReferenceShorterThanRead)
{
    FLAGS_CpG = true;

    statics_dicts_t statics_dict{};
    real_data_t real_data{"TT", "C"};
    std::vector<double> modern_model_deam(2, 0.01);
    std::vector<double> ancient_model_deam(2, 0.01);
    std::string quals = "II";
    std::string maskedseq;

    statics_denominator_table_t denom(static_cast<size_t>(FLAGS_range));
    EXPECT_NO_THROW(calPMD(std::move(real_data), modern_model_deam, ancient_model_deam, quals, maskedseq, statics_dict, denom));
}
