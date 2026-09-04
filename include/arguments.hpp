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
#include <string>
#include <string_view>
#include <vector>
#include <span>
#include "alignment/alignment.hpp"
#include "parsedData.hpp"

// X-Macro
#define DEFINE_bool(name, default_val, help_text) \
    extern bool FLAGS_##name; \
    extern bool IS_USED_##name;

#define DEFINE_int32(name, default_val, help_text) \
    extern int FLAGS_##name;        \
    extern bool IS_USED_##name;

#define DEFINE_string(name, default_val, help_text) \
    extern std::string FLAGS_##name; \
    extern bool IS_USED_##name;

#define DEFINE_POSITION_LIST(name, default_val, help_text)                     \
    extern std::string FLAGS_##name;                                             \
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
#undef DEFINE_POSITION_LIST
#undef DEFINE_double

enum class CustomTerminusStatus {
  MATCHED,
  NOT_MATCHED,
  INVALID_POSITION
};

void initCMDParse(int argc, char *argv[]);

// need for implement
void function_basecomposition(parsedData &data, alignnmentData_t &real_alignmentData);

bool function_basicFilter(std::string_view real_ref_seq, std::string_view real_read, std::size_t read_len, std::string_view qual_seq);

bool function_basicTerminal(std::string_view real_read, std::string_view real_ref_seq, std::string_view qual_seq);

void function_in_thread_pool_maskterminaldeam_or_maskterminalbases(const std::string &maskedseq, std::string &real_read, bool is_reverse, const std::vector<std::string> &splited_line, std::string &origin_line);

/**
 * @brief Checks custom terminal positions for deamination mismatches.
 *
 * The positions specified by @c --customterminus are examined in their
 * original order. Non-negative positions are indexed from the beginning of
 * the aligned reference sequence. Negative positions are resolved relative
 * to the end of the aligned reference sequence.
 *
 * At non-negative positions, the function detects C-to-T mismatches. At
 * negative positions, it detects C-to-T mismatches when @c --ss is enabled
 * and G-to-A mismatches otherwise. A mismatch is accepted only when the
 * corresponding Phred+33 quality score is greater than or equal to
 * @c --requirebaseq.
 *
 * If at least one qualifying mismatch is found, the SAM record is written to
 * standard output. If no qualifying mismatch is found, the trailing newline
 * is removed from @p line and the @c LS:Z:0 tag is appended.
 *
 * @param[in] real_read
 *     Aligned read sequence containing the observed bases.
 * @param[in] real_ref_seq
 *     Aligned reference sequence containing the reference bases. Its length
 *     is used to resolve negative custom positions.
 * @param[in] quals
 *     Phred+33 quality characters corresponding to aligned sequence
 *     positions.
 * @param[in,out] line
 *     SAM record associated with the sequences. The record is not modified
 *     when a match is found, but is written to standard output. When no match
 *     is found, its trailing newline is removed and @c LS:Z:0 is appended.
 *     It remains unchanged if a requested position is invalid.
 *
 * @retval CustomTerminusStatus::MATCHED
 *     At least one requested position contains a qualifying deamination
 *     mismatch. The SAM record has been written to standard output.
 * @retval CustomTerminusStatus::NOT_MATCHED
 *     All requested positions are valid, but none contains a qualifying
 *     mismatch. The @c LS:Z:0 tag has been appended to @p line.
 * @retval CustomTerminusStatus::INVALID_POSITION
 *     At least one resolved position falls outside @p real_ref_seq,
 *     @p real_read, or @p quals. The SAM record is neither modified nor
 *     written.
 *
 * @pre The custom positions have been parsed by inputParams_validator().
 * @pre Quality characters use Phred+33 encoding.
 *
 * @note A negative position is resolved using the length of
 *       @p real_ref_seq, and the resulting absolute index is then used for
 *       all three input sequences.
 *
 * @see CustomTerminusStatus
 * @see get_customTerminusPositions()
 */
[[nodiscard]]
CustomTerminusStatus function_customterminus(std::string_view real_read,
                             std::string_view real_ref_seq,
                             std::string_view quals, std::string &line);

[[nodiscard]]
std::span<const int> get_customTerminusPositions() noexcept;

namespace ParamChecks{
inline bool isUsing_customterminus() { return IS_USED_customterminus; }

inline bool isUsing_maskTerminalBases() { return IS_USED_maskterminalbases; }

inline bool isUsing_maskTerminalDeaminations() {
  return IS_USED_maskterminaldeaminations;
    }
}