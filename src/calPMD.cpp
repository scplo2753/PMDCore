#include "calPMD.hpp"
#include "arguments.hpp"
#include <algorithm>

/**
 * @brief Constructs a calPMD object and initializes its members based on the provided parameters. It calculates the PMD score by iterating through the real read and reference sequence, applying the ancient and modern deamination models, and updating mismatch dictionaries accordingly.
 * @details The constructor takes in a real_data_t object containing the real read and reference sequence, vectors for modern and ancient deamination models, quality scores, a masked sequence, and a statics_dicts_t object for managing mismatch dictionaries. It initializes the member variables and calls the calPMD_loop() function to perform the PMD calculation.
 * @param[in] real_data A real_data_t object containing the real read and reference sequence
 * @param[in] modern_model_deam A vector of doubles representing the modern deamination model
 * @param[in] ancient_model_deam A vector of doubles representing the ancient deamination model
 * @param[in] quals A string_view representing the quality scores of the read
 * @param[in] maskedseq_input A string representing the masked sequence
 * @param[in] statics_dict A statics_dicts_t object for managing mismatch dictionaries
 * @note The constructor assumes that the input sequences and quality scores are valid and properly formatted.
 */
calPMD::calPMD(real_data_t &&real_data, const std::vector<double> &modern_model_deam, const std::vector<double> &ancient_model_deam, std::string_view quals, const std::string &maskedseq_input, statics_dicts_t &statics_dict) : real_read(std::move(real_data.real_read)),
                                                                                                                                                                                                                      real_ref_seq(std::move(real_data.real_ref_seq)),
                                                                                                                                                                                                                      quals(quals),
                                                                                                                                                                                                                      temp_quals(quals),
                                                                                                                                                                                                                      ancient_model_deam(ancient_model_deam),
                                                                                                                                                                                                                      modern_model_deam(modern_model_deam),
                                                                                                                                                                                                                      mismatch_dict(statics_dict.mismatch_dict),
                                                                                                                                                                                                                      mismatch_dict_CpG(statics_dict.mismatch_dict_CpG),
                                                                                                                                                                                                                      mismatch_dict_rev(statics_dict.mismatch_dict_rev),
                                                                                                                                                                                                                      mismatch_dict_CpG_rev(statics_dict.mismatch_dict_CpG_rev),
                                                                                                                                                                                                                      maskedseq(maskedseq_input),
                                                                                                                                                                                                                      statics_dict(statics_dict)
{
    real_read_length = std::min({real_read.length(), real_ref_seq.length(), quals.length()});
    start_pos = 0;
    backStart_pos = real_read.length() > 0 ? real_read.length() - 1 : 0;

    if (FLAGS_CpG)
    {
        if (real_ref_seq.find_first_of("CG") == string_view::npos)
        {
            L_MD.L_D = 1.0;
            L_MD.L_M = 1.0;
        }
    }
    // else if (FLAGS_UDGhalf)
    //{
    //     /* code */
    // }
    else if (real_ref_seq.find('C') == std::string::npos && real_ref_seq.find('G') == std::string::npos)
    {
        L_MD.L_D = 1.0;
        L_MD.L_M = 1.0;
    }
    //---end----
    calPMD_loop();
}

/**
 * @brief This function iterates through the real read and reference sequence, calculating the PMD score based on the provided ancient and modern deamination models. It updates mismatch dictionaries accordingly and computes the degradation score for each position in the read.
 * @details The function loops through each position in the real read, checking for valid bases and
 * applying the ancient and modern deamination models to compute the degradation score. It also updates mismatch dictionaries based on the observed mismatches and their positions. The function handles both forward and reverse sequences, taking into account the specified flags for CpG context and other options.
 * @note The function assumes that the input sequences and quality scores are valid and properly formatted.
 */
void calPMD::calPMD_loop()
{
    int addition = 0;
    // a=real_read_pos
    // b=real_ref_seq_pos
    std::string qualsRev(quals.data(), quals.size());
    std::reverse(qualsRev.begin(), qualsRev.end());

    for (int site = 0; site < real_read_length; ++site)
    {
        if (real_read[site] == 'N' || real_ref_seq[site] == 'N' || real_ref_seq[site] == '-')
            continue;
        int start_distance = site - start_pos;
        int backStart_distance = backStart_pos - site;

        // if FLAGS_adjustbaseq_all
        // if FLAGS_adjustss:

        if (quals[site] - 33 < FLAGS_requirebaseq)
        {
            ///@todo implement if FLAGS_adjustbaseq
            continue;
        }

        if(FLAGS_platypus)
        {
            platypus(start_distance, backStart_distance, real_ref_seq[site], real_read[site], addition);
        }
        ///@todo implement options.deamination

        int result = computeDegradationScore(start_distance, backStart_distance, real_ref_seq[site], real_read[site], qualsRev);
        if (result == -1)
            continue;
        else if (result == -2)
            break;

        //Compute degradation score
    }
}

/**
 * @brief This function processes a site in the read and updates the mismatch dictionaries based on the provided parameters.
 * @param[in] start_distance The distance from the start of the read
 * @param[in] backStart_distance The distance from the end of the read
 * @param[in] real_ref_seq_pos The reference sequence base at the current position
 * @param[in] real_read_pos The read base at the current position
 * @param[in] addition The number of times this mismatch has been observed
 * @returns None
 */
void calPMD::platypus(const int &start_distance, const int &backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos,int addition)
{
    bool CpGcheck = false;
    std::string the_key = "";
    if (start_distance + 1 < real_ref_seq.length())
    {
        if (real_ref_seq.at(start_distance) == 'C' && real_ref_seq.at(start_distance + 1) == 'G')
        {
            CpGcheck = true;
        }
    }
    the_key.push_back(real_ref_seq_pos);
    the_key.push_back(real_read_pos);
    the_key += std::to_string(start_distance);

    {
        std::lock_guard<std::mutex> lock(statics_dict.dict_mutex); //only lock when modifying the dictionary
        if (CpGcheck == true)
        {
            if (mismatch_dict_CpG.find(the_key) != mismatch_dict_CpG.end())
            {
                addition = mismatch_dict_CpG[the_key];
                addition += 1;
                mismatch_dict_CpG[the_key] = addition;
            }
            else
            {
                mismatch_dict_CpG[the_key] = 1;
            }
        }
        else
        {
            if (mismatch_dict.find(the_key) != mismatch_dict.end())
            {
                addition = mismatch_dict[the_key];
                addition += 1;
                mismatch_dict[the_key] = addition;
            }
            else
            {
                mismatch_dict[the_key] = 1;
            }
        }
    }

    // count 3' end
    CpGcheck = false;
    if (start_distance > 0) // make sure not out of range
    {
        if (real_ref_seq.at(start_distance - 1) == 'C' && real_ref_seq.at(start_distance) == 'G')
        {
            CpGcheck = true;
        }
    }
    the_key.clear();
    the_key.push_back(real_ref_seq_pos);
    the_key.push_back(real_read_pos);
    the_key += std::to_string(backStart_distance);

    // set lock only when modifying the dictionary
    {
        std::lock_guard<std::mutex> lock(statics_dict.dict_mutex);
        if (CpGcheck == true)
        {
            if (mismatch_dict_CpG_rev.find(the_key) != mismatch_dict_CpG_rev.end())
            {
                addition = mismatch_dict_CpG_rev[the_key];
                addition += 1;
                mismatch_dict_CpG_rev[the_key] = addition;
            }
            else
            {
                mismatch_dict_CpG_rev[the_key] = 1;
            }
        }
        else
        {
            if (mismatch_dict_rev.find(the_key) != mismatch_dict_rev.end())
            {
                addition = mismatch_dict_rev[the_key];
                addition += 1;
                mismatch_dict_rev[the_key] = addition;
            }
            else
            {
                mismatch_dict_rev[the_key] = 1;
            }
        }
    }
}

/**
* @returns -1 skip current loop
*        -2 break the loop
*        0 like normal
*/
int calPMD::computeDegradationScore(int start_distance,int backStart_distance,const char &real_ref_seq_pos,const char &real_read_pos, std::string &qualsRev)
{
    if (start_distance >= real_read_length)
        return -1;
    if(real_ref_seq_pos=='C')
    {
        if(FLAGS_CpG)
        {
            if (start_distance + 1 >= real_ref_seq.length())
                return -2;
            if (real_ref_seq.at(start_distance + 1) != 'G')
                return -1;
        }
        ///@todo implement else if UDGhalf
        ///@todo implement else if EcoliCpG
        ///@todo implement else if Ecoli
        
        if(real_read_pos=='T')
        {
            L_MD.L_D = L_MD.L_D * L_mismatch(start_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_mismatch(start_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);

            if(FLAGS_ss)
            {
                L_MD.L_D = L_MD.L_D * L_mismatch_ss(start_distance, backStart_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
                L_MD.L_M = L_MD.L_M * L_mismatch_ss(start_distance, backStart_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);
            }
            ///@todo if options.adjustbaseq
            function_maskterminaldeam_init_maskedseq(start_distance, backStart_distance, false);
        }
        else if (real_read_pos=='C')
        {
            L_MD.L_D = L_MD.L_D * L_match(start_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_match(start_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);

            if(FLAGS_ss)
            {
                L_MD.L_D = L_MD.L_D * L_match_ss(start_distance, backStart_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
                L_MD.L_M = L_MD.L_M * L_match_ss(start_distance, backStart_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);
            }
        }
        //if options.PMDSprim
    }
    if (real_ref_seq_pos == 'G' && FLAGS_ss == false)
    {
        if (FLAGS_CpG)
        {
            if (start_distance == 0)
                return -1;
            if (real_ref_seq[start_distance - 1] != 'C')
                return -1;
        }
        //if options.UDGhalf
        if (real_read_pos == 'A')
        {
            L_MD.L_D = L_MD.L_D * L_mismatch(backStart_distance, ancient_model_deam, qualsRev, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_mismatch(backStart_distance, modern_model_deam, qualsRev, FLAGS_polymorphism_contamination);
            function_maskterminaldeam_init_maskedseq(start_distance, backStart_distance, true);
        }
        //if options.maskterminaldeaminations != False and options.ss ==False:
        //else
        else if (real_read_pos == 'G')
        {
            L_MD.L_D = L_MD.L_D * L_match(backStart_distance, ancient_model_deam, qualsRev, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_match(backStart_distance, modern_model_deam, qualsRev, FLAGS_polymorphism_contamination);
        }
    }

    LR = std::log(L_MD.L_D / L_MD.L_M);

    quals = temp_quals;
    return 0;
}

/**
 * @brief This function checks if the calculated likelihood ratio (LR) falls within the specified threshold range. It returns true if the LR is greater than or equal to the lower threshold and less than the upper threshold, indicating that the PMD score meets the filtering criteria.
 * @returns true if the LR is within the threshold range, false otherwise
 */
bool calPMD::threshold_filter()
{
    if (LR >= FLAGS_threshold && LR < FLAGS_upperthreshold)
    {
        return true;
    }
    return false;
}

/**
 * @brief This function initializes the masked sequence based on the provided start and back start distances, as well as the reverse context flag. It modifies the masked sequence by replacing bases with 'N' at specified positions if certain conditions are met, such as being within the threshold for masking terminal deaminations.
 * @param[in] start_distance The distance from the start of the read
 * @param[in] backstart_distance The distance from the end of the read
 * @param[in] is_reverse_context A boolean indicating if the sequence is in reverse orientation
 * @note The function checks the FLAGS_maskterminaldeaminations and FLAGS_ss flags to determine
 */
void calPMD::function_maskterminaldeam_init_maskedseq(int start_distance, int backstart_distance, bool is_reverse_context)
{
    if (!IS_USED_maskterminaldeaminations)
    {
        return;
    }

    if (start_distance <= FLAGS_maskterminaldeaminations)
    {
        maskedseq = real_read.substr(0, start_distance) + "N" + real_read.substr(start_distance + 1);
    }
    else if (backstart_distance <= FLAGS_maskterminaldeaminations && (FLAGS_ss || is_reverse_context))
    {
        maskedseq = real_read.substr(0, start_distance) + "N" + real_read.substr(start_distance + 1);
    }
}