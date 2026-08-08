#pragma once
#include "utility.hpp"
#include "statics_types.hpp" 
#include <string>
#include <string_view>
#include <map>

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
     * @param[in] statics_dict A statics_dicts_t object for managing mismatch dictionaries
     * @param[in] statics_nucleo_total_table_t object for managing nucleotide totals
     * @note The constructor assumes that the input sequences and quality scores are valid and properly formatted.
     */
    calPMD(real_data_t &&real_data, 
        const std::vector<double> &modern_model_deam, const std::vector<double> &ancient_model_deam, 
        std::string_view quals, 
        const std::string &maskedseq, 
        statics_dicts_t &statics_dict,
        statics_denominator_table_t &nucleo_total_table);

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
    statics_dicts_t& statics_dict;  // to read and update mismatch dictionaries

    statics_denominator_table_t& statics_denominator_table; // to read and update nucleotide totals

    L_MD_t L_MD;

    double LR;

    void calPMD_loop();
    void platypus_forward(size_t start_distance, const char &real_ref_seq_pos, const char &real_read_pos);
    void platypus_backward(size_t start_distance,size_t backstart_distance, const char &real_ref_seq_pos, const char &real_read_pos);
    void platypus(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos);
    int computeDegradationScore(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos, std::string &qualsRev);
    void function_maskterminaldeam_init_maskedseq(size_t start_distance, size_t backstart_distance, bool is_reverse_context);
    std::vector<double>* choose_nucleo_total_table_vector(const char &base,statics_nucleo_total_table_t &denominator_table);
};
