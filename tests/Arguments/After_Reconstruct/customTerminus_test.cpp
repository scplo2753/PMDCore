/**
 * @file customTerminus_test.cpp
 *
 * | Test case | Positions | Reference | Read | Expected status |
 * |---|---|---|---|---|
 * | Non-negative position matches | `0` | `CAA` | `TAA` | `MATCHED` |
 * | Negative position matches | `-1` | `AAG` | `AAA` | `MATCHED` |
 * | Valid positions do not match | `0,-1` | `CAG` | `CAG` | `NOT_MATCHED` |
 * | Match precedes out-of-range position | `0,100` | `CAA` | `TAA` | `INVALID_POSITION` |
 * | Out-of-range position precedes match | `100,0` | `CAA` | `TAA` | `INVALID_POSITION` |
 * | Non-match precedes out-of-range position | `0,100` | `CAA` | `CAA` | `INVALID_POSITION` |
 * | All positions are out of range | `100,-100` | `CAA` | `TAA` | `INVALID_POSITION` |
 * | Maximum non-negative position matches | `2` | `AAC` | `AAT` | `MATCHED` |
 * | Non-negative position equals length | `3` | `AAC` | `AAT` | `INVALID_POSITION` |
 * | Minimum negative position matches | `-3` | `GAA` | `AAA` | `MATCHED` |
 * | Negative position exceeds length | `-4` | `GAA` | `AAA` | `INVALID_POSITION` |
 */
#include "arguments.hpp"
#include <algorithm>
#include <array>
#include <gtest/gtest.h>
#include <span>
#include <string>
#include <string_view>
#include <vector>

class custom_terminus_test : public ::testing::Test {
protected:
  void set_custom_terminus(std::string value,
                           bool single_stranded = false) {
    char program[] = "PMDCore";
    char option[] = "--customterminus";

    arguments_.clear();
    arguments_.push_back(program);
    arguments_.push_back(option);
    arguments_.push_back(value);

    if (single_stranded) {
      arguments_.push_back("--ss");
    }

    argv_.clear();
    for (auto &argument : arguments_) {
      argv_.push_back(argument.data());
    }

    initCMDParse(static_cast<int>(argv_.size()), argv_.data());
  }

  void SetUp() override {}

  void TearDown() override {
    FLAGS_customterminus.clear();
    IS_USED_customterminus = false;
    FLAGS_ss = false;
    IS_USED_ss = false;
  }

  std::string non_negative_ct_match_line =
      "ct_non_negative_ct_match_01	0	chr1	1	60	20M	"
      "*	0	0	TACGTGATCCGATGCTAGCA	"
      "IIIIIIIIIIIIIIIIIIII	MD:Z:0C19";
  std::string last_non_negative_ct_match_line =
      "ct_last_non_negative_ct_match_01	0	chr1	1	60	20M	"
      "*	0	0	TACGTGATCCGATGCTAGCT	"
      "IIIIIIIIIIIIIIIIIIII	MD:Z:19C0";
  std::string negative_ga_match_line =
      "ct_negative_ga_match_01	16	chr1	1	60	20M	"
      "*	0	0	TACGTGATCCGATGCTAGCA	"
      "IIIIIIIIIIIIIIIIIIII	MD:Z:0C19";
  std::string minimum_negative_ga_match_line =
      "ct_minimum_negative_ga_match_01	0	chr1	1	60	20M	"
      "*	0	0	AACGTGATCCGATGCTAGCT	"
      "IIIIIIIIIIIIIIIIIIII	MD:Z:0G19";
  std::string no_mismatch_line =
      "ct_no_mismatch_01	0	chr1	1	60	20M	"
      "*	0	0	TACGTGATCCGATGCTAGCA	"
      "IIIIIIIIIIIIIIIIIIII	MD:Z:20";
  std::vector<std::string> arguments_;
  std::vector<char *> argv_;
};

TEST_F(custom_terminus_test, parses_option_after_custom_terminus) {
  set_custom_terminus("0,1,-1,-2", true);
  const std::span<const int> terminus_positions = get_customTerminusPositions();

  EXPECT_TRUE(IS_USED_customterminus);
  EXPECT_TRUE(FLAGS_ss);

  constexpr std::array<int, 4> expected{0, 1, -1, -2};
  ASSERT_EQ(terminus_positions.size(), expected.size());
  EXPECT_TRUE(std::ranges::equal(terminus_positions, expected));
}

TEST_F(custom_terminus_test, parses_custom_terminus_positions) {
  set_custom_terminus("0,1,-1,-2");
  const std::span<const int> terminus_positions = get_customTerminusPositions();

  EXPECT_TRUE(IS_USED_customterminus);
  EXPECT_FALSE(FLAGS_customterminus.empty());
  EXPECT_FALSE(terminus_positions.empty());

  constexpr std::array<int, 4> expected{0, 1, -1, -2};

  ASSERT_EQ(terminus_positions.size(), expected.size());
  EXPECT_TRUE(std::ranges::equal(terminus_positions, expected));
}

TEST_F(custom_terminus_test, matches_non_negative_position) {
  set_custom_terminus("0");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "CACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  testing::internal::CaptureStdout();

  CustomTerminusStatus result =
      function_customterminus(real_read, reference_sequence, quality,
                              non_negative_ct_match_line);

  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result, CustomTerminusStatus::MATCHED);
  EXPECT_EQ(output, non_negative_ct_match_line + "\n");
}

TEST_F(custom_terminus_test, matches_negative_position_without_ss) {
  set_custom_terminus("-1");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "TACGTGATCCGATGCTAGCG";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  testing::internal::CaptureStdout();

  const CustomTerminusStatus result =
      function_customterminus(real_read, reference_sequence, quality,
                              negative_ga_match_line);

  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result, CustomTerminusStatus::MATCHED);
  EXPECT_EQ(output, negative_ga_match_line + "\n");
}

TEST_F(custom_terminus_test, returns_not_matched_for_valid_positions) {
  set_custom_terminus("0,-1");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "TACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  testing::internal::CaptureStdout();

  std::string origin = no_mismatch_line;

  CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality, no_mismatch_line);

  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_TRUE(output.empty());
  EXPECT_EQ(result, CustomTerminusStatus::NOT_MATCHED);
  EXPECT_EQ(no_mismatch_line, origin + "\tLS:Z:0");
}

TEST_F(custom_terminus_test, returns_invalid_when_out_of_range_follows_match) {
  set_custom_terminus("0,1,-100");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "CACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result =
      function_customterminus(real_read, reference_sequence, quality,
                              non_negative_ct_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}

TEST_F(custom_terminus_test, returns_invalid_when_out_of_range_precedes_match) {
  set_custom_terminus("-100,0");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "CACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result =
      function_customterminus(real_read, reference_sequence, quality,
                              non_negative_ct_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}

TEST_F(custom_terminus_test,
       returns_invalid_when_non_match_precedes_out_of_range) {
  set_custom_terminus("1,100");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "CACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality, non_negative_ct_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}

TEST_F(custom_terminus_test, returns_invalid_when_all_positions_out_of_range) {
  set_custom_terminus("99,100");

  std::string_view real_read = "TACGTGATCCGATGCTAGCA";
  std::string_view reference_sequence = "CACGTGATCCGATGCTAGCA";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality, non_negative_ct_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}

TEST_F(custom_terminus_test, matches_maximum_non_negative_position) {
  set_custom_terminus("19");

  std::string_view real_read = "TACGTGATCCGATGCTAGCT";
  std::string_view reference_sequence = "TACGTGATCCGATGCTAGCC";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  testing::internal::CaptureStdout();

  const CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality,
      last_non_negative_ct_match_line);

  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result, CustomTerminusStatus::MATCHED);
  EXPECT_EQ(output, last_non_negative_ct_match_line + "\n");
}

TEST_F(custom_terminus_test,
       returns_invalid_for_non_negative_position_at_length) {
  set_custom_terminus("20");

  std::string_view real_read = "TACGTGATCCGATGCTAGCT";
  std::string_view reference_sequence = "TACGTGATCCGATGCTAGCC";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality,
      last_non_negative_ct_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}

TEST_F(custom_terminus_test, matches_minimum_negative_position) {
  set_custom_terminus("-20");

  std::string_view real_read = "AACGTGATCCGATGCTAGCT";
  std::string_view reference_sequence = "GACGTGATCCGATGCTAGCT";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  testing::internal::CaptureStdout();

  const CustomTerminusStatus result = function_customterminus(
      real_read, reference_sequence, quality,
      minimum_negative_ga_match_line);

  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result, CustomTerminusStatus::MATCHED);
  EXPECT_EQ(output, minimum_negative_ga_match_line + "\n");
}

TEST_F(custom_terminus_test,
       returns_invalid_for_negative_position_beyond_length) {
  set_custom_terminus("-21");

  std::string_view real_read = "AACGTGATCCGATGCTAGCT";
  std::string_view reference_sequence = "GACGTGATCCGATGCTAGCT";
  std::string_view quality = "IIIIIIIIIIIIIIIIIIII";

  CustomTerminusStatus result =
      function_customterminus(real_read, reference_sequence, quality,
                              minimum_negative_ga_match_line);
  EXPECT_EQ(result, CustomTerminusStatus::INVALID_POSITION);
}
