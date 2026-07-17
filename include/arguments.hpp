/**
 * @file arguments.hpp
 * @author 2753
 * @brief this file is used for define cmdline arguments and implement some basic filter functions
 * @version 0.1
 * @date 2026-07-17
 * 
 * 
 */
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

bool function_basicFilter(std::string_view real_ref_seq, std::string_view real_read, uint read_len, std::string_view qual_seq);

bool function_basicTerminal(std::string_view real_read, std::string_view real_ref_seq, std::string_view qual_seq);

void function_in_thread_pool_maskterminaldeam_or_maskterminalbases(const std::string &maskedseq, std::string &real_read, bool is_reverse, const std::vector<std::string> &splited_line, std::string &origin_line);