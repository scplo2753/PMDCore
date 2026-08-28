#pragma once
#include "alignment/alignment.hpp"
#include "statistics/statistics_types.hpp"
#include "statistics/deam_types.hpp"
#include <string>
#include <string_view>
#include <map>
#include <vector>

using std::map;
using std::string;
using std::string_view;

struct L_MD_t
{
    double L_D = 1.0;
    double L_M = 1.0;
    double L_D_max = 1.0;
    double L_M_max = 1.0;
    double L_D_min = 1.0;
    double L_M_min = 1.0;
};

/**
 * @class calPMD
 * @brief This class is responsible for calculating the PMD (Post Mortem Damage) score based on the provided real read, reference sequence, quality scores, and other parameters. It utilizes the ancient and modern deamination models to compute the degradation score and update mismatch dictionaries accordingly.
 * 
 */
class calPMD
{
public:
    /**
     * @brief Constructs a calPMD object and initializes its members based on the provided parameters.
     *        It calculates the PMD score by iterating through the real read and reference sequence, applying the ancient and modern deamination models, and updating mismatch dictionaries accordingly.
     * @param[in] real_data A real_data_t object containing the real read and reference sequence
     * @param[in] modern_model_deam A vector of doubles representing the modern deamination model
     * @param[in] ancient_model_deam A vector of doubles representing the ancient deamination model
     * @param[in] quals A string_view representing the quality scores of the read
     * @param[in] maskedseq_input A string representing the masked sequence
     * @param[in] platypus_statics_dict A platypus_statics_dicts_t object for managing mismatch dictionaries
     * @param[in,out] platypus_denominator_table Platypus nucleotide totals
     * @param[in,out] deam_statics_table Deamination counts for the current worker
     * @note The constructor assumes that the input sequences and quality scores are valid and properly formatted.
     */
    calPMD(real_data_t &&real_data, 
        const std::vector<double> &modern_model_deam, const std::vector<double> &ancient_model_deam, 
        std::string_view quals, 
        const std::string &maskedseq, 
        platypus_statics_dicts_t &platypus_statics_dict,
        platypus_denominator_table_t &platypus_denominator_table,
        deamination_statics_t &deam_statics_table);

    ~calPMD() = default;

    calPMD(const calPMD &) = delete;
    calPMD& operator=(const calPMD&) = delete;
    calPMD(calPMD&&) = delete;
    calPMD& operator=(calPMD&&) = delete;

    /**
     * @brief This function checks if the calculated likelihood ratio (LR) falls within the specified threshold range. 
     *        It returns true if the LR is greater than or equal to the lower threshold and less than the upper threshold, indicating that the PMD score meets the filtering criteria.
     * @returns true if the LR is within the threshold range, false otherwise
     */
    bool threshold_filter();
    const std::string& get_maskedSeq() const noexcept { return maskedseq; }

private:
    std::string real_read;
    std::string real_ref_seq;
    std::string quals;
    std::string temp_quals;
    std::string maskedseq;
    const std::vector<double> &ancient_model_deam;
    const std::vector<double> &modern_model_deam;
    size_t start_pos;
    size_t backStart_pos;
    size_t real_read_length;

    match_dict_t& mismatch_dict_CpG;
    match_dict_t& mismatch_dict;
    match_dict_t& mismatch_dict_CpG_rev;
    match_dict_t& mismatch_dict_rev;
    platypus_statics_dicts_t& platypus_statics_dict;  // to read and update mismatch dictionaries

    platypus_denominator_table_t &platypus_denominator_table; // to read and update nucleotide totals
    deamination_statics_t &deam_statics_table;

    L_MD_t L_MD;

    double LR;

    void calPMD_loop();

    /**
     * @brief Records the 5-prime part of the refactored Platypus statistics.
     *
     * This function implements the forward half of the reference platypus()
     * operation and is called directly by calPMD_loop().
     *
     * @param start_distance Zero-based distance from the 5-prime end.
     * @param real_ref_seq_pos Reference base at the current aligned position.
     * @param real_read_pos Observed read base at the current aligned position.
     */
    void platypus_forward(size_t start_distance, const char &real_ref_seq_pos, const char &real_read_pos);

    /**
     * @brief Records the 3-prime part of the refactored Platypus statistics.
     *
     * This function implements the reverse half of the reference platypus()
     * operation and is called directly by calPMD_loop().
     *
     * @param start_distance Zero-based distance from the 5-prime end, used for
     *        checking the preceding reference base.
     * @param backStart_distance Zero-based distance from the 3-prime end.
     * @param real_ref_seq_pos Reference base at the current aligned position.
     * @param real_read_pos Observed read base at the current aligned position.
     */
    void platypus_backward(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos);

    /**
     * @brief Records Platypus statistics for both the 5-prime and 3-prime ends
     *        using the original PMDtools processing flow.
     *
     * This function intentionally preserves the structure and operation order of
     * the original implementation. It serves as a reference implementation for
     * validating the results produced by the refactored platypus_forward() and
     * platypus_backward() functions.
     *
     * The function updates both forward and reverse mismatch dictionaries and
     * nucleotide denominator tables. CpG observations are placed in their
     * corresponding CpG-specific tables.
     *
     * @param start_distance Zero-based distance from the 5-prime end.
     * @param backStart_distance Zero-based distance from the 3-prime end.
     * @param real_ref_seq_pos Reference base at the current aligned position.
     * @param real_read_pos Observed read base at the current aligned position.
     *
     * @note This is a reference implementation rather than the implementation used
     *       by calPMD_loop(). Preserve its original control flow when modifying the
     *       refactored Platypus counting code.
     */
    void platypus(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos);

    /**
     * @brief Records a reference-C or reference-G observation in the
     *        deamination statistics table.
     *
     * Reference-C observations are recorded in the forward table using their
     * zero-based distance from the 5-prime end. Reference-G observations are
     * recorded in the reverse table using their zero-based distance from the
     * 3-prime end.
     *
     * Context filtering depends on the active mode:
     * - FLAGS_CpG records only C positions followed by G and G positions
     *   preceded by C.
     * - FLAGS_noCpG excludes C positions followed by G and G positions
     *   preceded by C.
     * - When neither mode is enabled, all eligible reference-C and reference-G
     *   positions are recorded.
     *
     * Positions outside the corresponding statistics table and observations with
     * unsupported read bases are skipped without stopping iteration.
     *
     * A terminal reference C cannot be classified by its following-base context.
     * When FLAGS_CpG or FLAGS_noCpG is enabled, encountering such a position
     * causes the caller to stop processing the current read.
     *
     * @param start_distance Zero-based distance from the 5-prime end and index of
     *        the current aligned position.
     * @param backStart_distance Zero-based distance from the 3-prime end.
     * @param real_ref_seq_pos Reference base at the current aligned position.
     * @param real_read_pos Observed read base at the current aligned position.
     *
     * @retval true Continue processing subsequent aligned positions.
     * @retval false Stop processing the current read because the CpG context of a
     *         terminal reference C cannot be determined.
     */
    bool deamination(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos);

    /**
     * @brief Updates the ancient and modern likelihoods for one aligned position.
     *
     * For a reference C, the function evaluates C-to-T mismatches and C-to-C
     * matches using the distance from the 5-prime end. For a reference G in
     * double-stranded mode, it evaluates G-to-A mismatches and G-to-G matches
     * using the distance from the 3-prime end.
     *
     * FLAGS_CpG restricts evaluation to CpG-context positions, while FLAGS_noCpG
     * excludes CpG-context positions. When FLAGS_ss is enabled, reference-G
     * processing is disabled and the single-stranded likelihood contribution is
     * additionally applied to eligible reference-C observations.
     *
     * Eligible deamination mismatches may also update maskedseq when terminal
     * deamination masking is enabled.
     *
     * @param start_distance Zero-based distance from the 5-prime end and index of
     *        the current aligned position.
     * @param backStart_distance Zero-based distance from the 3-prime end.
     * @param real_ref_seq_pos Reference base at the current aligned position.
     * @param real_read_pos Observed read base at the current aligned position.
     * @param qualsRev Reversed quality-score string used for reference-G
     *        likelihood calculations.
     *
     * @retval 0 The position was handled normally. This can also mean that the
     *         reference/read base combination required no likelihood update.
     * @retval -1 Skip the current position and continue processing the read.
     * @retval -2 Stop processing the current read because the following-base
     *         context required for a terminal reference C is unavailable.
     */
    int computeDegradationScore(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos, std::string &qualsRev);

    /**
     * @brief Masks an eligible terminal deamination position in maskedseq.
     *
     * The function has no effect unless terminal-deamination masking was
     * explicitly enabled through IS_USED_maskterminaldeaminations.
     *
     * The current aligned position is masked with 'N' when it is within
     * FLAGS_maskterminaldeaminations of the 5-prime end. It is also masked when
     * it is within that distance of the 3-prime end and either single-stranded
     * mode is enabled or the observation is being processed in a reverse
     * deamination context.
     *
     * @param start_distance Zero-based distance from the 5-prime end and index
     *        of the position to mask in maskedseq.
     * @param backStart_distance Zero-based distance from the 3-prime end.
     * @param is_reverse_context true when called for a reverse-context
     *        deamination observation, such as a reference-G/read-A mismatch;
     *        false for a forward-context observation.
     *
     * @pre start_distance must be a valid maskedseq index whenever the position
     *      satisfies the masking conditions.
     */
    void function_maskterminaldeam_init_maskedseq(size_t start_distance, size_t backStart_distance, bool is_reverse_context);

    /**
     * @brief Selects the nucleotide-total vector associated with a base.
     *
     * @param base Nucleotide whose vector should be selected.
     * @param nucleotide_total_table Table containing the A, C, G, and T total
     *        vectors.
     *
     * @return Pointer to the vector corresponding to base.
     * @retval nullptr base is not one of 'A', 'C', 'G', or 'T'.
     */
    std::vector<double>* choose_nucleo_total_table_vector(const char &base,statics_nucleo_total_table_t &nucleotide_total_table);
};
