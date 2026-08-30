#include <gtest/gtest.h>
#include "utilities/string_utils.hpp"
#include <vector>

using std::vector;

TEST(StringToInt,normal)
{
    int zero = 99999;
    int negative = 9999;
    int postive = -1;
    int minusZero = -10;
    vector<StringToIntStatus> status(4);
    status[0] = StringToInt("0", zero);
    status[1] = StringToInt("-1", negative);
    status[2] = StringToInt("1", postive);
    status[3] = StringToInt("-0", minusZero);
    EXPECT_EQ(status[0], StringToIntStatus::SUCCESS);
    EXPECT_EQ(status[1], StringToIntStatus::SUCCESS);
    EXPECT_EQ(status[2], StringToIntStatus::SUCCESS);
    EXPECT_EQ(status[3], StringToIntStatus::SUCCESS);
    EXPECT_EQ(zero, 0);
    EXPECT_EQ(negative, -1);
    EXPECT_EQ(postive, 1);
    EXPECT_EQ(minusZero, 0);
}

TEST(StringToInt,try_to_parse_illegal_number_type)
{
    int floatNumber = 0;
    int floatNumber2 = 0;
    int scienceNumber = 0;
    vector<StringToIntStatus> status(3);
    status[0] = StringToInt("1.234", floatNumber);
    status[1] = StringToInt(".5", floatNumber2);
    status[2] = StringToInt("1e-5", scienceNumber);
    EXPECT_EQ(status[0], StringToIntStatus::TRAILING_CHARACTER);
    EXPECT_EQ(status[1], StringToIntStatus::INVALID_CHARACTER);
    EXPECT_EQ(status[2], StringToIntStatus::TRAILING_CHARACTER);
}

TEST(StringToInt,not_numeric_string)
{
    int alpha = INT_MAX;
    int blank = INT_MAX;
    StringToIntStatus alpha_status = StringToInt("a", alpha);
    StringToIntStatus blank_status = StringToInt(" ", blank);
    EXPECT_EQ(alpha_status, StringToIntStatus::INVALID_CHARACTER);
    EXPECT_EQ(blank_status, StringToIntStatus::INVALID_CHARACTER);
    EXPECT_EQ(alpha, INT_MAX);
    EXPECT_EQ(blank, INT_MAX);
}

TEST(StringToInt, empty_string)
{
    int result = INT_MAX;
    StringToIntStatus status = StringToInt("", result);
    EXPECT_EQ(status, StringToIntStatus::EMPTY_INPUT);
    EXPECT_EQ(result, INT_MAX);
}

TEST(StringToInt, out_of_range)
{
    int intmax = INT_MIN;
    int intmin = INT_MAX;
    int plusmax = 123;
    int minusmin = 456;
    const std::string above_max = std::to_string(static_cast<long long>(INT_MAX) + 1LL);
    const std::string below_min = std::to_string(static_cast<long long>(INT_MIN) - 1LL);
    StringToIntStatus status[4];

    status[0] = StringToInt(std::to_string(INT_MAX), intmax);
    status[1] = StringToInt(std::to_string(INT_MIN), intmin);
    status[2] = StringToInt(above_max, plusmax);
    status[3] = StringToInt(below_min, minusmin);

    EXPECT_EQ(status[0], StringToIntStatus::SUCCESS);
    EXPECT_EQ(status[1], StringToIntStatus::SUCCESS);
    EXPECT_EQ(intmax, INT_MAX);
    EXPECT_EQ(intmin, INT_MIN);
    EXPECT_EQ(status[2], StringToIntStatus::OUT_OF_RANGE);
    EXPECT_EQ(status[3], StringToIntStatus::OUT_OF_RANGE);
    EXPECT_EQ(plusmax, 123);
    EXPECT_EQ(minusmin, 456);
}