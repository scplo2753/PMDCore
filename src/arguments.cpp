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

// if return true, display currespond line
bool function_basicFilter(std::string_view real_ref_seq, std::string_view real_read, std::string_view qual_seq)
{
    if (!(FLAGS_basic > 0))
    {
        return false;
    }
    int limit = std::min<int>({(int)real_ref_seq.size(), (int)real_read.size(), FLAGS_basic});
    for (int pos = 0; pos < limit; ++pos)
    {
        if (real_read[pos] == 'N' ||
            real_ref_seq[pos] == 'N' ||
            real_read[pos] == '-' ||
            real_ref_seq[pos] == '-')
            return false;
        if (real_ref_seq[pos] == 'C' && real_read[pos] == 'T' && qual_seq[pos] - 33 > FLAGS_requirebaseq)
        {
            if (FLAGS_CpG)
            {
                if (pos + 1 < real_ref_seq.size() && real_ref_seq[pos + 1] == 'G')
                    return true;
            }
            else
            {
                return true;
            }
        }
    }
    return false;
}

bool function_basicTerminalFilter(std::string_view real_read, std::string_view real_ref_seq, std::string_view qual_seq)
{
    const int threshold = FLAGS_requirebaseq + 33;
    char a = real_read[0];
    char b = real_ref_seq[0];
    if (a == 'T' && b == 'C' && (qual_seq[0] > threshold)) [[unlikely]]
    {
        return true;
    }
    a = real_read.back();
    b = real_ref_seq.back();
    if (a == 'A' && b == 'G' && (*qual_seq.rbegin() > threshold)) [[unlikely]]
    {
        return true;
    }
    return false;
}