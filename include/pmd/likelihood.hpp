#pragma once
#include <cstddef>
#include <string_view>
#include <vector>
#include "utilities/probability_utils.hpp"


/**
 * @brief Computes the probability of a matching observation
 * @details implements Eq. (1)
 * P. Skoglund, B.H. Northoff, M.V. Shunkov, A.P. Derevianko, S. Pääbo, J. Krause, & M. Jakobsson,
 *  Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal,
 * Proc. Natl. Acad. Sci. U.S.A. 111 (6) 2229-2234,
 * https://doi.org/10.1073/pnas.1318934111 (2014).
 *
 * @param fpos
 * @param fmodel
 * @param fqual
 * @param polymorphism
 * @return double
 */
inline double L_match(std::size_t fpos, const std::vector<double> &fmodel, std::string_view fqual, double polymorphism)
{
    double P_damage = fmodel[fpos];
    double P_error = phred2prob(fqual[fpos] - 33) / 3.0;
    const double P_poly = polymorphism;
    double P_match = (1.0 - P_damage) * (1.0 - P_error) * (1.0 - P_poly) + (P_damage * P_error * (1.0 - P_poly)) + (P_error * P_poly * (1.0 - P_damage));
    return P_match;
}

inline double L_mismatch(std::size_t fpos, const std::vector<double> &fmodel, std::string_view fqual, double polymorphism)
{
    return 1.0f - L_match(fpos, fmodel, fqual, polymorphism);
}

inline double L_match_ss(std::size_t fpos, std::size_t zpos, const std::vector<double> &fmodel, std::string_view fqual, double polymorphism)
{
    double P_damage = fmodel[fpos];
    double P_damage2 = fmodel[zpos];
    double P_error = phred2prob(fqual[fpos] - 33) / 3.0;
    const double P_poly = polymorphism;
    double P_match = (1.0 - P_damage) * (1.0 - P_error) * (1.0 - P_poly) * (1.0 - P_damage2) + (P_damage * P_error * (1.0 - P_poly) * (1.0 - P_damage2)) + (P_damage2 * P_error * (1.0 - P_poly) * (1.0 - P_damage)) + (P_error * P_poly * (1.0 - P_damage2) * (1.0 - P_damage));
    return P_match;
}

inline double L_mismatch_ss(std::size_t fpos, std::size_t zpos, const std::vector<double> &fmodel, std::string_view fqual, double polymorphism)
{
    return 1.0f - L_match_ss(fpos, zpos, fmodel, fqual, polymorphism);
}
