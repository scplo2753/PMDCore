#pragma once
#include "utility.hpp"
#include <string>
#include <string_view>
#include <map>

using std::map;
using std::string;
using std::string_view;

struct L_MD_t
{
    float L_D = 1.0;
    float L_M = 1.0;
    float L_D_max = 1.0;
    float L_M_max = 1.0;
    float L_D_min = 1.0;
    float L_M_min = 1.0;
};

/**
 * @class calPMD
 * @brief This class is responsible for calculating the PMD (Post Mortem Damage) score based on the provided real read, reference sequence, quality scores, and other parameters. It utilizes the ancient and modern deamination models to compute the degradation score and update mismatch dictionaries accordingly.
 * 
 */
class calPMD
{
public:
    calPMD(real_data_t &&real_data, const std::vector<double> &modern_model_deam, const std::vector<double> &ancient_model_deam, std::string_view quals, const std::string &maskedseq, statics_dicts_t &statics_dict);
    ~calPMD() = default;

    bool threshold_filter();

private:
    std::string real_read;
    std::string real_ref_seq;
    std::string quals;
    std::string temp_quals;
    std::string maskedseq;
    const std::vector<double> &ancient_model_deam;
    const std::vector<double> &modern_model_deam;
    int start_pos;
    int backStart_pos;
    int real_read_length;

    match_dict_t& mismatch_dict_CpG;
    match_dict_t& mismatch_dict;
    match_dict_t& mismatch_dict_CpG_rev;
    match_dict_t& mismatch_dict_rev;
    statics_dicts_t& statics_dict;  // 为了存取mutex

    L_MD_t L_MD;

    double LR;

    void calPMD_loop();
    void platypus(const int &start_distance, const int &backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos,int addition);
    int computeDegradationScore(int start_distance, int backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos, std::string &qualsRev);
    void function_maskterminaldeam_init_maskedseq(int start_distance, int backstart_distance, bool is_reverse_context);
};