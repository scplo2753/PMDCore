#include "Filters.hpp"
#include <string_view>

/**
 * @brief This function checks if the length of the read sequence falls within the specified minimum and maximum length thresholds. It returns true if the read length is valid, and false otherwise.
 * @param[in] data The recordLine_struct_t structure containing the read sequence and other alignment information.
 * @returns true if the read length is valid, false otherwise
 */
inline bool lengthFilter(const recordLine_struct_t &data)
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

bool cigarFilter(const recordLine_struct_t &data)
{
    // Single pass over the CIGAR string to collect presence flags.
    bool hasI = false, hasD = false, hasS = false;
    bool hasH = false, hasN = false, hasP = false;
    const std::string &c = data.cigar;
    for (char ch : c)
    {
        switch (ch)
        {
        case 'I': hasI = true; break;
        case 'D': hasD = true; break;
        case 'S': hasS = true; break;
        case 'H': hasH = true; break;
        case 'N': hasN = true; break;
        case 'P': hasP = true; break;
        default: break;
        }
        // early exit if all flags are set
        if (hasI && hasD && hasS && hasH && hasN && hasP)
            break;
    }

    bool hasHNP = hasH || hasN || hasP;

    if (FLAGS_noinsertion && hasI) return false;
    if (FLAGS_nodeletion && hasD) return false;
    if (FLAGS_onlyinsertions && !hasI) return false;
    if (FLAGS_onlydeletions && !hasD) return false;
    if (FLAGS_noindels && (hasI || hasD)) return false;
    if (FLAGS_noclips && (hasS || hasHNP)) return false;
    if (FLAGS_onlyclips && !hasS) return false;

    if (hasHNP)
    {
        std::cout << "Not support M, I, S, D, this alignment will be skipped." << std::endl;
        return false;
    }

    return true;
}

/**
 * @brief This function applies a series of filters to the alignment data based on the specified command-line flags. It checks the read length, chromosome name, CIGAR string, and mapping quality against the provided thresholds. If any of the filters fail, it returns false; otherwise, it returns true.
 * @param[in] data The recordLine_struct_t structure containing the read sequence and other alignment information.
 * @returns true if the alignment data passes all filters, false otherwise
 */
bool ArgsFilter(const recordLine_struct_t &data)
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