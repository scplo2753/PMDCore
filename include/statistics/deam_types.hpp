#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

constexpr std::size_t forward_index(char base)
{
    switch (base)
    {
    case 'T':
        return 0;
    case 'A':
        return 1;
    case 'G':
        return 2;
    case 'C':
        return 3;
    default:
        return 4;
    }
}

constexpr std::size_t reverse_index(char base)
{
    switch (base)
    {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    default:
        return 4;
    }
}

struct deamination_direction_counts_t
{
    std::vector<std::array<std::uint64_t, 4>> counts;
    std::vector<std::uint64_t> totals;

    explicit deamination_direction_counts_t(std::size_t range)
        : counts(range, std::array<std::uint64_t, 4>{0, 0, 0, 0}),
          totals(range, std::uint64_t{0})
    {}

    deamination_direction_counts_t &operator+=(const deamination_direction_counts_t &right)
    {
        if (counts.size() != right.counts.size() ||
            totals.size() != right.totals.size())
        {
            throw std::invalid_argument("cannot merge deamination counts with different ranges");
        }

        for (std::size_t row = 0; row < counts.size(); ++row)
        {
            for (std::size_t col = 0; col < counts[row].size(); ++col)
            {
                counts[row][col] += right.counts[row][col];
            }
            totals[row] += right.totals[row];
        }

        return *this;
    }
};

struct deamination_statics_t
{
    deamination_direction_counts_t forward;
    deamination_direction_counts_t reverse;

    explicit deamination_statics_t(std::size_t range) : forward(range), reverse(range)
    {}

    deamination_statics_t &operator+=(const deamination_statics_t &right)
    {
        forward += right.forward;
        reverse += right.reverse;
        return *this;
    }
};
