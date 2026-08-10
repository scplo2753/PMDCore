#pragma once
#include <string>
#include <string_view>

inline constexpr char complement(char nucleotide)
{
    switch (nucleotide)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    case 'a':
        return 't';
    case 't':
        return 'a';
    case 'c':
        return 'g';
    case 'g':
        return 'c';
    case 'R':
        return 'Y';
    case 'Y':
        return 'R';
    case 'N':
        return 'N';
    case '-':
        return '-';
    default:
        return 'N';
    }
}

inline std::string revcomp(std::string_view input)
{
    std::string output;
    output.reserve(input.size());
    for (auto it = input.rbegin(); it != input.rend(); ++it)
    {
        output += complement(*it);
    }
    return output;
}
