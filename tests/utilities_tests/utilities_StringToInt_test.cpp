#include <gtest/gtest.h>
#include <climits>
#include "utilities/string_utils.hpp"
#include <vector>

using std::vector;

TEST(ParseInteger,normal)
{
    int zero = 99999;
    int negative = 9999;
    int postive = -1;
    int minusZero = -10;
    vector<ParseIntegerStatus> status(4);
    status[0] = ParseInteger("0", zero);
    status[1] = ParseInteger("-1", negative);
    status[2] = ParseInteger("1", postive);
    status[3] = ParseInteger("-0", minusZero);
    EXPECT_EQ(status[0], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(status[1], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(status[2], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(status[3], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(zero, 0);
    EXPECT_EQ(negative, -1);
    EXPECT_EQ(postive, 1);
    EXPECT_EQ(minusZero, 0);
}

TEST(ParseInteger,try_to_parse_illegal_number_type)
{
    int floatNumber = 0;
    int floatNumber2 = 0;
    int scienceNumber = 0;
    vector<ParseIntegerStatus> status(3);
    status[0] = ParseInteger("1.234", floatNumber);
    status[1] = ParseInteger(".5", floatNumber2);
    status[2] = ParseInteger("1e-5", scienceNumber);
    EXPECT_EQ(status[0], ParseIntegerStatus::TRAILING_CHARACTER);
    EXPECT_EQ(status[1], ParseIntegerStatus::INVALID_CHARACTER);
    EXPECT_EQ(status[2], ParseIntegerStatus::TRAILING_CHARACTER);
}

TEST(ParseInteger,not_numeric_string)
{
    int alpha = INT_MAX;
    int blank = INT_MAX;
    ParseIntegerStatus alpha_status = ParseInteger("a", alpha);
    ParseIntegerStatus blank_status = ParseInteger(" ", blank);
    EXPECT_EQ(alpha_status, ParseIntegerStatus::INVALID_CHARACTER);
    EXPECT_EQ(blank_status, ParseIntegerStatus::INVALID_CHARACTER);
    EXPECT_EQ(alpha, INT_MAX);
    EXPECT_EQ(blank, INT_MAX);
}

TEST(ParseInteger, empty_string)
{
    int result = INT_MAX;
    ParseIntegerStatus status = ParseInteger("", result);
    EXPECT_EQ(status, ParseIntegerStatus::EMPTY_INPUT);
    EXPECT_EQ(result, INT_MAX);
}

TEST(ParseInteger, out_of_range)
{
    int intmax = INT_MIN;
    int intmin = INT_MAX;
    int plusmax = 123;
    int minusmin = 456;
    const std::string above_max = std::to_string(static_cast<long long>(INT_MAX) + 1LL);
    const std::string below_min = std::to_string(static_cast<long long>(INT_MIN) - 1LL);
    ParseIntegerStatus status[4];

    status[0] = ParseInteger(std::to_string(INT_MAX), intmax);
    status[1] = ParseInteger(std::to_string(INT_MIN), intmin);
    status[2] = ParseInteger(above_max, plusmax);
    status[3] = ParseInteger(below_min, minusmin);

    EXPECT_EQ(status[0], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(status[1], ParseIntegerStatus::SUCCESS);
    EXPECT_EQ(intmax, INT_MAX);
    EXPECT_EQ(intmin, INT_MIN);
    EXPECT_EQ(status[2], ParseIntegerStatus::OUT_OF_RANGE);
    EXPECT_EQ(status[3], ParseIntegerStatus::OUT_OF_RANGE);
    EXPECT_EQ(plusmax, 123);
    EXPECT_EQ(minusmin, 456);
}
