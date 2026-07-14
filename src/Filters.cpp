#include "Filters.hpp"
#include <string_view>

inline bool lengthFilter(const AlignLine_Data_t &data)
{
    if (FLAGS_readlength > 0 && data.read_seq.size() != FLAGS_readlength)
    {
        return false;
    }
    if (FLAGS_maxlength > 0 && data.read_seq.size() > FLAGS_maxlength)
    {
        return false;
    }
    if (FLAGS_minlength > 0 && data.read_seq.size() < FLAGS_minlength)
    {
        return false;
    }
    return true;
}

bool cigarFilter(const AlignLine_Data_t &data)
{
    if (FLAGS_noinsertion && data.cigar.find('I') != std::string::npos)
    {
        return false;
    }
    if (FLAGS_nodeletion && data.cigar.find('D') != std::string::npos)
    {
        return false;
    }
    if (FLAGS_onlyinsertions && data.cigar.find('I') == std::string::npos)
    {
        return false;
    }
    if (FLAGS_onlydeletions && data.cigar.find('D') == std::string::npos)
    {
        if (FLAGS_noinsertion && data.cigar.find('I') != std::string::npos)
        {
            return false;
        }
        if (FLAGS_nodeletion && data.cigar.find('D') != std::string::npos)
        {
            return false;
        }
        if (FLAGS_onlyinsertions && data.cigar.find('I') == std::string::npos)
        {
            return false;
        }
        if (FLAGS_onlydeletions && data.cigar.find('D') == std::string::npos)
        {
            return false;
        }
        if (FLAGS_noindels && (data.cigar.find('I') != std::string::npos || data.cigar.find('D') != std::string::npos))
        {
            return false;
        }
        if (FLAGS_noclips && (data.cigar.find('S') != std::string::npos ||
                              data.cigar.find('H') != std::string::npos ||
                              data.cigar.find('N') != std::string::npos ||
                              data.cigar.find('P') != std::string::npos))
        {
            return false;
        }
        if (FLAGS_onlyclips && data.cigar.find('S') == std::string::npos)
        {
            return false;
        }
        if (data.cigar.find('H') != std::string::npos ||
            data.cigar.find('N') != std::string::npos ||
            data.cigar.find('P') != std::string::npos)
        {
            std::cout << "Not support M, I, S, D, this alignment will be skipped." << std::endl;
            return false;
        }
        return false;
    }
    if (FLAGS_noindels && (data.cigar.find('I') != std::string::npos || data.cigar.find('D') != std::string::npos))
    {
        return false;
    }
    if (FLAGS_noclips && (data.cigar.find('S') != std::string::npos ||
                          data.cigar.find('H') != std::string::npos ||
                          data.cigar.find('N') != std::string::npos ||
                          data.cigar.find('P') != std::string::npos))
    {
        return false;
    }
    if (FLAGS_onlyclips && data.cigar.find('S') == std::string::npos)
    {
        return false;
    }
    if (data.cigar.find('H') != std::string::npos ||
        data.cigar.find('N') != std::string::npos ||
        data.cigar.find('P') != std::string::npos)
    {
        std::cout << "Not support M, I, S, D, this alignment will be skipped." << std::endl;
        return false;
    }
    return true;
}

bool ArgsFilter(const AlignLine_Data_t &data)
{
    if (!lengthFilter(data))
    {
        return false;
    }
    if (!FLAGS_chromosome.empty() && data.RNAME != FLAGS_chromosome)
    {
        return false;
    }

    if (!cigarFilter(data))
    {
        return false;
    }

    if (FLAGS_requiremapq > 0 && data.MAPQ < FLAGS_requiremapq) // maybe need to validate the MAPQ value
    {
        return false;
    }
    return true;
}

bool badRefSeq_Vailder(std::string_view real_ref_seq, std::string_view line)
{
    if (std::find(real_ref_seq.begin(), real_ref_seq.end(), 'A') == real_ref_seq.end())
        if (std::find(real_ref_seq.begin(), real_ref_seq.end(), 'T') == real_ref_seq.end())
            if (std::find(real_ref_seq.begin(), real_ref_seq.end(), 'C') == real_ref_seq.end())
                if (std::find(real_ref_seq.begin(), real_ref_seq.end(), 'G') == real_ref_seq.end())
                {
                    std::cout << "bad reference sequence reconstruction: " << real_ref_seq << std::endl;
                    std::cout << "SAM line: " << line << std::endl;
                    return false;
                }
    return true;
}