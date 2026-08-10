#pragma once
#include <cstddef>
#include <vector>
#include "utilities/probability_utils.hpp"

/**
 * @brief geometric distribution
 *
 * @param[in] pval success probability, in PMDS it's PMD param
 * @param[in] kval pos, pos is 1 are not defined
 * @param[in] constant PMD constant
 * @return double
 */
inline double geometric(double pval, double kval, double constant)
{
    return (std::pow((1.0 - pval),(kval - 1.0))) * pval + constant;
}

inline std::vector<double> get_modern_model_deam(std::size_t length, double probability)
{
    return std::vector<double>(length, probability);
}

/**
 * @brief Get the vector of probabilities of ancient model deam object
 *
 * @param[in] pval FLAGS_PMDpparam
 * @param[in] constant FLAGS_PMDconstant
 * @return std::vector<double>
 */
inline std::vector<double> get_ancient_model_deam(double pval, double constant)
{
    std::vector<double> ancient_model_deam(999);
    for (int pos = 1; pos < 1000; ++pos)
    {
        //ancient_model_deam.at(pos - 1) = geometric(FLAGS_PMDpparam, pos, FLAGS_PMDconstant);
        ancient_model_deam.at(pos - 1) = geometric(pval, pos, constant);
    }
    return ancient_model_deam;
}
