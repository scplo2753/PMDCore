#pragma once
#include <cmath>

/**
 * @brief error probabilities P to Phred quality score Q 
 * \f$-10*\log_{10}
 * @param[in] P error probabilities 
 * @return double output phred qulity score
 */
static double prob2phred(double P)
{
    return -10.0 * std::log10(P);
}

static double phred2prob(double Q)
{
    return std::pow(10.0, -Q / 10.0);
}

static double geometric(double pval,double kval,double constant)
{
    return (std::pow((1.0 - pval),(kval - 1.0))) * pval + constant;
}

static double py_round(double value, int ndigits = 0)
{
    if (std::isnan(value) || std::isinf(value) || value == 0.0)
    {
        return value;
    }

    long double factor = std::pow(10.0L, static_cast<long double>(ndigits));
    long double shifted = static_cast<long double>(value) * factor;
    long double rounded = std::nearbyint(shifted);
    return static_cast<double>(rounded / factor);
}