#include <gtest/gtest.h>
#include "statics.hpp"

TEST(StaticsTest, split_key_splits_correctly) {
    std::string key = "CT5";
    std::string name;
    int index;

    split_key(key, name, index);

    EXPECT_EQ(name, "CT");
    EXPECT_EQ(index, 5);
}