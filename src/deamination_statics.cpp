#include "deamination_statics.hpp"
#include <array>
#include <iostream>
#include <string_view>
#include "utilities/math_utils.hpp"

constexpr std::array<std::string_view, 8> pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG"};

void print_deamination_statistics(const deamination_statics_t &statistics)
{
    const std::size_t range = statistics.forward.counts.size();
    if (statistics.forward.totals.size() != range ||
        statistics.reverse.counts.size() != range ||
        statistics.reverse.totals.size() != range)
    {
        throw std::invalid_argument(
            "deamination statistics have inconsistent ranges");
    }

    std::cout << "z\t";
    for (const auto &pair : pairs)
        std::cout << pair << '\t';
    std::cout << '\n';

    for (std::size_t row = 0; row < range; ++row)
    {
        std::cout << row << '\t';
        const std::uint64_t forward_total = statistics.forward.totals[row];
        const std::uint64_t reverse_total = statistics.reverse.totals[row];

        for (const auto count : statistics.forward.counts[row])
        {
            if (forward_total == 0 || count == 0)
                std::cout << "0.00000\t";
            else
                std::cout << py_round(
                                 static_cast<double>(count) /
                                     static_cast<double>(forward_total),
                                 5)
                          << '\t';
        }

        for (const auto count : statistics.reverse.counts[row])
        {
            if (reverse_total == 0 || count == 0)
                std::cout << "0.00000\t";
            else
                std::cout << py_round(
                                 static_cast<double>(count) /
                                     static_cast<double>(reverse_total),
                                 5)
                          << '\t';
        }
        std::cout << '\n';
    }
}
