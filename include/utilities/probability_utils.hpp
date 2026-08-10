#pragma once
#include <cmath>

/**
 * @brief error probabilities P to Phred quality score Q
 * \f$-10*\log_{10}\f$
 * @param[in] P error probabilities
 * @return double output phred qulity score
 */
inline double prob2phred(double P)
{
    return -10.0 * std::log10(P);
}

/**
 * @brief Phred quality score Q to error probabilities P
 *
 * @param[in] Q Phred quality score
 * @return double output error probabilities
 */
inline double phred2prob(double Q)
{
    return std::pow(10.0, -Q / 10.0);
}