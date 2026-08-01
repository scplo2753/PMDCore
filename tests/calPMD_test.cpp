#include <gtest/gtest.h>
#include "calPMD.hpp"

class calPMDTest : public ::testing::Test {
protected:
    calPMD *calPMD_ptr;
    void SetUp() override{
        calPMD_ptr=new calPMD()
    }
}