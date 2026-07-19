#include "Filters.hpp"
#include <string_view>

/**
 * @brief This function checks if the length of the read sequence falls within the specified minimum and maximum length thresholds. It returns true if the read length is valid, and false otherwise.
 * @param[in] data The AlignLine_Data_t structure containing the read sequence and other alignment information.
 * @returns true if the read length is valid, false otherwise
 */
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

/**
 * @brief This function applies a series of filters to the alignment data based on the specified command-line flags. It checks the read length, chromosome name, CIGAR string, and mapping quality against the provided thresholds. If any of the filters fail, it returns false; otherwise, it returns true.
 * @param[in] data The AlignLine_Data_t structure containing the read sequence and other alignment information.
 * @returns true if the alignment data passes all filters, false otherwise
 */
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

/**
 * @brief This function checks if the reference sequence is valid by ensuring it contains at least one of the nucleotide bases (A, T, C, G). If the reference sequence is invalid, it prints an error message along with the corresponding SAM line and returns false; otherwise, it returns true.
 * @param[in] real_ref_seq The reference sequence to be validated.
 * @param[in] line The original SAM line associated with the reference sequence.
 * @returns true if the reference sequence is valid, false otherwise
 */
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