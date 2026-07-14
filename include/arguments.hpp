#pragma once
#include "argparse/argparse.hpp"
#include "utility.hpp"
#include "parsedData.hpp"

// X-Macro
#define DEFINE_bool(name, default_val, help_text) \
    extern bool FLAGS_##name;

#define DEFINE_int32(name, default_val, help_text) \
    extern int FLAGS_##name;        \
    extern bool IS_USED_##name;

#define DEFINE_string(name, default_val, help_text) \
    extern std::string FLAGS_##name; \
    extern bool IS_USED_##name;

#define DEFINE_double(name, default_val, help_text) \
    extern double FLAGS_##name;     \
    extern bool IS_USED_##name;
#define DEFINE_uint32(name, default_val, help_text) DEFINE_int32(name, default_val, help_text)

#include "args.list"

#undef DEFINE_uint32
#undef DEFINE_bool
#undef DEFINE_int32
#undef DEFINE_string
#undef DEFINE_double

void initCMDParse(int argc, char *argv[]);

// need for implement
void function_basecomposition(parsedData &data, alignnmentData_t &real_alignmentData);

bool function_basicFilter(std::string_view real_ref_seq, std::string_view real_read, std::string_view qual_seq);

bool function_basicTerminalFilter(std::string_view real_read, std::string_view real_ref_seq, std::string_view qual_seq);