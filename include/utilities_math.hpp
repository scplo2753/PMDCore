#pragma once
#include <cmath>

/**
 * @brief error probabilities P to Phred quality score Q 
 * \f$-10*\log_{10}\f$
 * @param[in] P error probabilities 
 * @return double output phred qulity score
 */
static double prob2phred(double P)
{
    return -10.0 * std::log10(P);
}

/**
 * @brief Phred quality score Q to error probabilities P
 * 
 * @param[in] Q Phred quality score
 * @return double output error probabilities
 */
static double phred2prob(double Q)
{
    return std::pow(10.0, -Q / 10.0);
}

/**
 * @brief geometric distribution
 * 
 * @param[in] pval success probability, in PMDS it's PMD param
 * @param[in] kval pos, pos is 1 are not defined
 * @param[in] constant PMD constant
 * @return double 
 */
static double geometric(const double pval,const double kval,const double constant)
{
    return (std::pow((1.0 - pval),(kval - 1.0))) * pval + constant;
}

static long double py_round(const double value,const int ndigits = 0)
{
    if (std::isnan(value) || std::isinf(value) || value == 0.0)
    {
        return value;
    }

    long double factor = std::pow(10.0L, static_cast<long double>(ndigits));
    long double shifted = static_cast<long double>(value) * factor;
    long double rounded = std::nearbyint(shifted);
    return (rounded / factor);
}

/**
 * @brief Get the vector of probabilities of ancient model deam object
 * 
 * @param[in] pval FLAGS_PMDpparam
 * @param[in] constant FLAGS_PMDconstant
 * @return std::vector<double> 
 */
static std::vector<double> get_ancient_model_deam(const double pval,const double constant)
{
    std::vector<double> ancient_model_deam(999);
    for (int pos = 1; pos < 1000; ++pos)
    {
        //ancient_model_deam.at(pos - 1) = geometric(FLAGS_PMDpparam, pos, FLAGS_PMDconstant);
        ancient_model_deam.at(pos - 1) = geometric(pval, pos, constant);
    }
    return ancient_model_deam;
}