#pragma once
#include <cmath>

inline long double py_round(const double value,const int ndigits = 0)
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