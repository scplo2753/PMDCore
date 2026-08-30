#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "utilities/string_utils.hpp"

using std::string;
using std::vector;

TEST(transformWhileSplit,normal)
{
    const string customterminus("0,1,2,-1,-2");
    vector<int> transformed{99,99,99,99};
    const auto origin = transformed;

    vector<int> expected{0, 1, 2, -1, -2};

    auto [index, offset, status] = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(status, StringToIntStatus::SUCCESS);

    EXPECT_EQ(index, expected.size());
    EXPECT_EQ(offset, customterminus.size());
    EXPECT_EQ(transformed, expected);

    EXPECT_NE(transformed, origin);
}

TEST(transformWhileSplit,with_space)
{
    const string customterminus("0, ,2, ,-2");
    vector<int> transformed{};

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::INVALID_CHARACTER);
    EXPECT_EQ(returns.field_index, 1);
    EXPECT_EQ(returns.offset, 2);

    EXPECT_TRUE(transformed.empty());
}

TEST(transformWhileSplit,single_number)
{
    const string customterminus("20");
    vector<int> transformed;
    vector<int> expected{20};

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.field_index, expected.size());
    EXPECT_EQ(returns.offset, customterminus.size());
    EXPECT_EQ(transformed, expected);
}

TEST(transformWhileSplit,empty_input)
{
    const string customterminus("");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::EMPTY_INPUT);

    EXPECT_EQ(returns.field_index, 0U);
    EXPECT_EQ(returns.offset, 0U);

    EXPECT_EQ(transformed, origin);
}

TEST(transformWhileSplit,empty_between_two_number)
{
    const string customterminus("1,,-2");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::EMPTY_INPUT);

    EXPECT_EQ(returns.field_index, 1U);
    EXPECT_EQ(returns.offset, 2U);

    EXPECT_EQ(transformed, origin);
}

TEST(transformWhileSplit,first_element_is_empty)
{
    const string customterminus(",-2");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::EMPTY_INPUT);

    EXPECT_EQ(returns.field_index, 0U);
    EXPECT_EQ(returns.offset, 0U);

    EXPECT_EQ(transformed, origin);
}

TEST(transformWhileSplit,last_element_is_empty)
{
    const string customterminus("-2,");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::EMPTY_INPUT);

    EXPECT_EQ(returns.field_index, 1U);
    EXPECT_EQ(returns.offset, 3U);

    EXPECT_EQ(transformed, origin);
}

TEST(transformWhileSplit,trailiing_character_in_middle)
{
    const string customterminus("-2,1.23,-1");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::TRAILING_CHARACTER);

    EXPECT_EQ(returns.field_index, 1U);
    EXPECT_EQ(returns.offset, 3U);

    EXPECT_EQ(transformed, origin);
}

TEST(transformWhileSplit,out_of_range)
{
    const string customterminus("0,999999999999999999999,-1");
    vector<int> transformed{20};
    const auto origin = transformed;

    auto returns = transformWhileSplit(customterminus, ',', transformed, StringToInt, StringToIntStatus::SUCCESS);

    EXPECT_EQ(returns.status, StringToIntStatus::OUT_OF_RANGE);

    EXPECT_EQ(returns.field_index, 1U);
    EXPECT_EQ(returns.offset, 2U);

    EXPECT_EQ(transformed, origin);
}

enum class floatMultipleStatus
{
    SUCCESS = 0,
    EMPTY_INPUT,
};

[[nodiscard]]
floatMultipleStatus multiple_f(std::string_view input,float &output) noexcept
{
    if(input.empty())
    {
        return floatMultipleStatus::EMPTY_INPUT;
    }
    float result{};
    auto [ptr, ec] = std::from_chars(input.data(), input.data() + input.size(), result);

    output = result * 2;
    return floatMultipleStatus::SUCCESS;
}

TEST(transformWhileSplit,custom_function)
{
    const string customterminus("5,0.3,0.1");
    vector<float> transformed{};

    auto returns = transformWhileSplit(customterminus, ',', transformed, multiple_f, floatMultipleStatus::SUCCESS);

    EXPECT_EQ(returns.status, floatMultipleStatus::SUCCESS);

    EXPECT_EQ(returns.field_index, 3);
    EXPECT_EQ(returns.offset, customterminus.size());
    EXPECT_EQ(transformed, (vector<float>{10, 0.6, 0.2}));
}