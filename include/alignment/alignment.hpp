#pragma once
#include <string>

/******
 * @brief struct of reconstructed reference and alignment data
 * @details first var is reference, the second one is read's alignment data
*/
struct alignnmentData_t
{
    std::string ref_seq; ///<reconstructed reference sequence
    std::string alignment; ///<read's alignment data
};

struct real_data_t
{
    std::string real_read;
    std::string real_ref_seq;
};
