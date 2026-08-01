#pragma once
#include "utility.hpp"
#include "arguments.hpp"

/**
 * @brief wrapper of get_ancient_model_deam, using FLAGS_PMDpparam and FLAGS_PMDconstant as pval and constant param
 * 
 * @return std::vector<double> ancient model deam
 */
static std::vector<double> ancientModelDeam_wrapper()
{
    return get_ancient_model_deam(FLAGS_PMDpparam, FLAGS_PMDconstant);
}