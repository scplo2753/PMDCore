#include "arguments.hpp"
#include <algorithm>
#include <ranges>

#define DEFINE_bool(name, default_val, help_text) \
    bool FLAGS_##name = default_val;

#define DEFINE_int32(name, default_val, help_text) \
    int FLAGS_##name = default_val; \
    bool IS_USED_##name=false;

#define DEFINE_string(name, default_val, help_text) \
    std::string FLAGS_##name = default_val; \
    bool IS_USED_##name=false;

#define DEFINE_double(name, default_val, help_text) \
    double FLAGS_##name = default_val;  \
    bool IS_USED_##name=false;

#define DEFINE_uint32(name, default_val, help_text) DEFINE_int32(name, default_val, help_text)

#include "args.list"

#undef DEFINE_uint32
#undef DEFINE_bool
#undef DEFINE_int32
#undef DEFINE_string
#undef DEFINE_double

void initCMDParse(int argc, char *argv[])
{
    argparse::ArgumentParser parser("PMDCore");

    #define DEFINE_bool(name, default_val, help_text)           \
        parser.add_argument("--" #name).help(help_text).flag(); \
        parser.add_argument("--no" #name).default_value(default_val).implicit_value(true);

    #define DEFINE_int32(name, default_val, help_text) \
        parser.add_argument("--" #name).default_value(default_val).help(help_text).scan<'i', int>(); \
        IS_USED_##name=parser.is_used("FLAGS_" #name);

    #define DEFINE_string(name, default_val, help_text) \
        parser.add_argument("--" #name).default_value(default_val).help(help_text); \
        IS_USED_##name=parser.is_used("FLAGS_" #name);

    #define DEFINE_double(name, default_val, help_text) \
        parser.add_argument("--" #name).default_value(default_val).help(help_text).scan<'g', double>(); \
        IS_USED_##name=parser.is_used("FLAGS_" #name);

    #define DEFINE_uint32(name, default_val, help_text) DEFINE_int32(name, default_val, help_text)

    #include "args.list"

    #undef DEFINE_uint32
    #undef DEFINE_bool
    #undef DEFINE_int32
    #undef DEFINE_string
    #undef DEFINE_double

        parser.parse_args(argc, argv);

    #define DEFINE_bool(name, default_val, help_text) \
        FLAGS_##name = default_val;                   \
        if (parser.get<bool>("--" #name))             \
        {                                             \
            FLAGS_##name = true;                      \
        }                                             \
        if (parser.get<bool>("--no" #name))           \
        {                                             \
            FLAGS_##name = false;                     \
        }

    #define DEFINE_int32(name, default_val, help_text) \
        FLAGS_##name = parser.get<int>("--" #name);

    #define DEFINE_uint32(name, default_val, help_text) DEFINE_int32(name, default_val, help_text)

    #define DEFINE_double(name, default_val, help_text) \
        FLAGS_##name = parser.get<double>("--" #name);

    #define DEFINE_string(name, default_val, help_text) \
        FLAGS_##name = parser.get<std::string>("--" #name);

    #include "args.list"

    #undef DEFINE_bool
    #undef DEFINE_int32
    #undef DEFINE_uint32
    #undef DEFINE_double
#undef DEFINE_string
}

void function_basecomposition(parsedData &data, alignnmentData_t &real_alignmentData)
{
    uint backoffset = 10;
    int endpos, startpos;
    int position = data.getPosition();
    int len_real_read = sizeof(real_alignmentData.ref_seq);
    if (data.isReversedSeq())
    {
        endpos = position;
        startpos = position + len_real_read;
    }
    else
    {
        endpos = position + len_real_read;
        startpos = position;
    }

    // statics of 5' end
}

/**
 * @brief Active if FLAGS_basic is set, corresponding to the first n bases of read, check if there is C>T mismatch with reference sequence and have sufficient base quality
 * 
 * @param[in] real_ref_seq real reference sequence
 * @param[in] real_read real read
 * @param[in] read_len length of read
 * @param[in] qual_seq  quality scores
 * @return true if first n bases have C>T mismatch with base and have sufficient base quality
 * @return false 
 */
bool function_basicFilter(std::string_view real_ref_seq, std::string_view real_read, uint read_len, std::string_view qual_seq)
{
    if (!(FLAGS_basic > 0))
    {
        return false;
    }

    const int limit = std::min<int>({static_cast<int>(real_ref_seq.size()), static_cast<int>(real_read.size()), FLAGS_basic});
    for (uint pos = 0; pos < limit; ++pos)
    {
        if (pos >= read_len)
        {
            return false;
        }

        const char read_base = real_read[pos];
        const char ref_base = real_ref_seq[pos];
        const char qual_char = pos < qual_seq.size() ? qual_seq[pos] : '\0';

        if (read_base == 'N' || ref_base == 'N' || read_base == '-' || ref_base == '-')
        {
            return false;
        }

        const bool quality_ok = qual_char != '\0' && (qual_char - 33) > FLAGS_requirebaseq;
        const bool is_ct_mismatch = ref_base == 'C' && read_base == 'T';

        if (FLAGS_CpG)
        {
            if (is_ct_mismatch && quality_ok && pos + 1 < real_ref_seq.size() && real_ref_seq[pos + 1] == 'G')
            {
                return true;
            }
        }
        else if (is_ct_mismatch && quality_ok)
        {
            return true;
        }
    }

    return false;
}

/**
 * @brief Active if FLAGS_terminal is set, prints the SAM line if a C>T mismatch with sufficient base quality is observed in the first and last base
 * 
 * @param[in] real_read real read
 * @param[in] real_ref_seq read reference sequence
 * @param[in] qual_seq quality scores
 * @return true if a C>T mismatch with sufficient base quality is observed in the first and last base
 * @return false 
 */
bool function_basicTerminal(std::string_view real_read, std::string_view real_ref_seq, std::string_view qual_seq)
{
    if (real_read.empty() || real_ref_seq.empty() || qual_seq.empty())
    {
        return false;
    }

    const int threshold = FLAGS_requirebaseq + 33;

    if (real_ref_seq[0] == 'C' && real_read[0] == 'T' && qual_seq[0] > threshold)
    {
        return true;
    }

    const char last_read = real_read.back();
    const char last_ref = real_ref_seq.back();
    const char last_qual = qual_seq.back();
    if (last_ref == 'G' && last_read == 'A' && last_qual > threshold)
    {
        return true;
    }

    return false;
}

void function_in_thread_pool_maskterminaldeam_or_maskterminalbases(const std::string &maskedseq, std::string &real_read, bool is_reverse, const std::vector<std::string> &splited_line, std::string &origin_line)
{
    std::string readp = maskedseq;
    if (is_reverse)
    {
        readp = revcomp(readp);
    }

    real_read = maskedseq;

    if (splited_line.empty())
    {
        return;
    }

    std::string rebuilt;
    rebuilt.reserve(maskedseq.size() + splited_line.size() * 8 + 32);

    const size_t first_count = std::min<size_t>(9, splited_line.size());
    for (size_t i = 0; i < first_count; ++i)
    {
        if (i > 0)
        {
            rebuilt.push_back('\t');
        }
        rebuilt += splited_line[i];
    }

    if (!first_count && splited_line.size() > 0)
    {
        rebuilt += splited_line[0];
    }

    rebuilt.push_back('\t');
    rebuilt += readp;

    if (splited_line.size() > 10)
    {
        rebuilt.push_back('\t');
        for (size_t i = 10; i < splited_line.size(); ++i)
        {
            if (i > 10)
            {
                rebuilt.push_back('\t');
            }
            rebuilt += splited_line[i];
        }
    }

    origin_line = std::move(rebuilt);
}