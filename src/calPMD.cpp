#include <cassert>
#include "calPMD.hpp"
#include "arguments.hpp"
#include <algorithm>

/***
 * @details The constructor takes in a real_data_t object containing the real read and reference sequence,
 *          vectors for modern and ancient deamination models, quality scores, a masked sequence, and a statics_dicts_t object for managing mismatch dictionaries.
 *          It initializes the member variables and calls the calPMD_loop() function to perform the PMD calculation.
 */
calPMD::calPMD(real_data_t &&real_data, const std::vector<double> &modern_model_deam, const std::vector<double> &ancient_model_deam, std::string_view quals, const std::string &maskedseq_input, statics_dicts_t &statics_dict, statics_denominator_table_t &denominator_table) : real_read(std::move(real_data.real_read)),
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
                                                                                                                                                                                                                                  statics_dict(statics_dict),
                                                                                                                                                                                                                                  statics_denominator_table(denominator_table)
{
    assert(quals.size() >= real_read.size());
    const bool masking_enabled =
        IS_USED_maskterminaldeaminations || IS_USED_maskterminalbases;
    assert(!masking_enabled || maskedseq.size() == real_read.size());

    real_read_length = std::min({real_read.size(), real_ref_seq.size()});
    if (ancient_model_deam.size() < real_read.size())
    {
        throw std::invalid_argument(
            "ancient deamination model is shorter than the read");
    }

    if (modern_model_deam.size() < real_read.size())
    {
        throw std::invalid_argument(
            "modern deamination model is shorter than the read");
    }
    start_pos = 0;
    backStart_pos = real_read.length() > 0 ? real_read.length() - 1 : 0;

    if (FLAGS_CpG)
    {
        if (real_ref_seq.find('C') == std::string::npos && real_ref_seq.find('G') == std::string::npos)
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

    LR = std::log(L_MD.L_D) - std::log(L_MD.L_M);
    assert(std::isfinite(LR));
    ///@todo implement if options.PMDSprim
    quals = temp_quals;
}

/**
 * @brief This function iterates through the real read and reference sequence, calculating the PMD score based on the provided ancient and modern deamination models.
 *        It updates mismatch dictionaries accordingly and computes the degradation score for each position in the read.
 * @details The function loops through each position in the real read, checking for valid bases and
 *          applying the ancient and modern deamination models to compute the degradation score.
 *          It also updates mismatch dictionaries based on the observed mismatches and their positions.
 *          The function handles both forward and reverse sequences, taking into account the specified flags for CpG context and other options.
 * @note The function assumes that the input sequences and quality scores are valid and properly formatted.
 */
void calPMD::calPMD_loop()
{
    // a=real_read_pos
    // b=real_ref_seq_pos
    std::string qualsRev(quals.data(), quals.size());
    std::reverse(qualsRev.begin(), qualsRev.end());

    for (size_t site = 0; site < real_read_length; ++site)
    {
        if (real_read[site] == 'N' || real_ref_seq[site] == 'N' || real_ref_seq[site] == '-')
            continue;
        size_t start_distance = site - start_pos;
        size_t backStart_distance = backStart_pos - site;

        ///@todo implement if FLAGS_adjustbaseq_all
        ///@todo implement if FLAGS_adjustss:

        if (quals[site] - 33 < FLAGS_requirebaseq)
        {
            ///@todo implement if FLAGS_adjustbaseq
            continue;
        }

        if (FLAGS_platypus)
        {
            platypus(start_distance, backStart_distance, real_ref_seq[site], real_read[site]);
        }
        ///@todo implement options.deamination

        const int result = computeDegradationScore(start_distance, backStart_distance, real_ref_seq[site], real_read[site], qualsRev);
        if (result == -1)
            continue;
        else if (result == -2)
            break;

        // Compute degradation score
    }
}

/**
 * @brief This function processes a site in the read and updates the mismatch dictionaries based on the provided parameters.
 * @param[in] start_distance The distance from the start of the read
 * @param[in] backStart_distance The distance from the end of the read
 * @param[in] real_ref_seq_pos The reference sequence base at the current position
 * @param[in] real_read_pos The read base at the current position
 */
void calPMD::platypus(const int &start_distance, const int &backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos)
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

    //count 5' end
    the_key.push_back(real_ref_seq_pos);
    the_key.push_back(real_read_pos);
    the_key += std::to_string(start_distance);
    
    {
        std::lock_guard<std::mutex> lock(statics_dict.dict_mutex); // only lock when modifying the dictionary
        if (CpGcheck == true)
        {
            ++mismatch_dict_CpG[the_key];
            std::vector<double> *nucleo_total_table_vector_ptr = choose_nucleo_total_table_vector(real_ref_seq_pos, statics_denominator_table.forward_CpG);
            if (nucleo_total_table_vector_ptr != nullptr && static_cast<size_t>(start_distance) < nucleo_total_table_vector_ptr->size())
            {
                nucleo_total_table_vector_ptr->at(start_distance) += 1.0;
            }
        }
        else
        {
            ++mismatch_dict[the_key];
            std::vector<double> *nucleo_total_table_vector_ptr = choose_nucleo_total_table_vector(real_ref_seq_pos, statics_denominator_table.forward);
            if (nucleo_total_table_vector_ptr != nullptr && static_cast<size_t>(start_distance) < nucleo_total_table_vector_ptr->size())
            {
                nucleo_total_table_vector_ptr->at(start_distance) += 1.0;
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

    //set lock only when modifying the dictionary
    {
        std::lock_guard<std::mutex> lock(statics_dict.dict_mutex);
        if (CpGcheck == true)
        {
            ++mismatch_dict_CpG_rev[the_key];
            std::vector<double> *nucleo_total_table_vector_ptr = choose_nucleo_total_table_vector(real_ref_seq_pos, statics_denominator_table.reverse_CpG);
            if (nucleo_total_table_vector_ptr != nullptr && static_cast<size_t>(backStart_distance) < nucleo_total_table_vector_ptr->size())
            {
                nucleo_total_table_vector_ptr->at(backStart_distance) += 1.0;
            }
        }
        else
        {
            ++mismatch_dict_rev[the_key];
            std::vector<double> *nucleo_total_table_vector_ptr = choose_nucleo_total_table_vector(real_ref_seq_pos, statics_denominator_table.reverse);
            if (nucleo_total_table_vector_ptr != nullptr && static_cast<size_t>(backStart_distance) < nucleo_total_table_vector_ptr->size())
            {
                nucleo_total_table_vector_ptr->at(backStart_distance) += 1.0;
            }
        }
    }
}

/**
 * @brief Compute the degradation score for a single read position.
 *
 * This function evaluates the current base and reference context, applies PMD
 * model likelihoods, updates the internal damage/match likelihood state, and
 * determines whether the current read position should be skipped or whether the
 * loop should be terminated.
 *
 * @param[in] start_distance Distance from the 5' end of the read.
 * @param[in] backStart_distance Distance from the 3' end of the read.
 * @param[in] real_ref_seq_pos Reference base at the current position.
 * @param[in] real_read_pos Read base at the current position.
 * @param[in,out] qualsRev Reversed quality string used for reverse-strand scoring.
 *
 * @return int
 *         -1 when the current position should be skipped and processing continues.
 *         -2 when processing should stop early.
 *          0 when the position was processed normally.
 */
int calPMD::computeDegradationScore(size_t start_distance, size_t backStart_distance, const char &real_ref_seq_pos, const char &real_read_pos, std::string &qualsRev)
{
    if (start_distance >= real_read_length)
        return -1;
    if (real_ref_seq_pos == 'C')
    {
        if (FLAGS_CpG)
        {
            if (start_distance + 1 >= real_read.size() ||
                start_distance + 1 >= real_ref_seq.size())
                return -2;
            if (real_ref_seq.at(start_distance + 1) != 'G')
                return -1;
        }
        ///@todo implement else if UDGhalf
        ///@todo implement else if EcoliCpG
        ///@todo implement else if Ecoli

        if (real_read_pos == 'T')
        {
            L_MD.L_D = L_MD.L_D * L_mismatch(start_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_mismatch(start_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);

            if (FLAGS_ss)
            {
                L_MD.L_D = L_MD.L_D * L_mismatch_ss(start_distance, backStart_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
                L_MD.L_M = L_MD.L_M * L_mismatch_ss(start_distance, backStart_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);
            }
            ///@todo if options.adjustbaseq
            function_maskterminaldeam_init_maskedseq(start_distance, backStart_distance, false);
        }
        else if (real_read_pos == 'C')
        {
            L_MD.L_D = L_MD.L_D * L_match(start_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_match(start_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);

            if (FLAGS_ss)
            {
                L_MD.L_D = L_MD.L_D * L_match_ss(start_distance, backStart_distance, ancient_model_deam, quals, FLAGS_polymorphism_ancient);
                L_MD.L_M = L_MD.L_M * L_match_ss(start_distance, backStart_distance, modern_model_deam, quals, FLAGS_polymorphism_contamination);
            }
        }
        // if options.PMDSprim
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
        // if options.UDGhalf
        if (real_read_pos == 'A')
        {
            L_MD.L_D = L_MD.L_D * L_mismatch(backStart_distance, ancient_model_deam, qualsRev, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_mismatch(backStart_distance, modern_model_deam, qualsRev, FLAGS_polymorphism_contamination);
            function_maskterminaldeam_init_maskedseq(start_distance, backStart_distance, true);
        }
        // if options.maskterminaldeaminations != False and options.ss ==False:
        // else
        else if (real_read_pos == 'G')
        {
            L_MD.L_D = L_MD.L_D * L_match(backStart_distance, ancient_model_deam, qualsRev, FLAGS_polymorphism_ancient);
            L_MD.L_M = L_MD.L_M * L_match(backStart_distance, modern_model_deam, qualsRev, FLAGS_polymorphism_contamination);
        }
    }


    return 0;
}

bool calPMD::threshold_filter()
{
    if (LR >= FLAGS_threshold && LR < FLAGS_upperthreshold)
    {
        return true;
    }
    return false;
}

/**
 * @brief This function initializes the masked sequence based on the provided start and back start distances, as well as the reverse context flag.
 *        It modifies the masked sequence by replacing bases with 'N' at specified positions if certain conditions are met,
 *        such as being within the threshold for masking terminal deaminations.
 * @param[in] start_distance The distance from the start of the read
 * @param[in] backstart_distance The distance from the end of the read
 * @param[in] is_reverse_context A boolean indicating if the sequence is in reverse orientation
 * @note The function checks the FLAGS_maskterminaldeaminations and FLAGS_ss flags to determine
 */
void calPMD::function_maskterminaldeam_init_maskedseq(size_t start_distance, size_t backstart_distance, bool is_reverse_context)
{
    if (!IS_USED_maskterminaldeaminations)
    {
        return;
    }

    const bool should_mask = start_distance <= FLAGS_maskterminaldeaminations ||
                             (backstart_distance <= FLAGS_maskterminaldeaminations &&
                              (FLAGS_ss || is_reverse_context));

    if(!should_mask)
        return;

    assert(start_distance < maskedseq.size());
    maskedseq[start_distance] = 'N';
}

std::vector<double>* calPMD::choose_nucleo_total_table_vector(const char &base,statics_nucleo_total_table_t &nucleo_total_table)
{
    switch (base)
    {
    case 'A':
        return &nucleo_total_table.A;
    case 'C':
        return &nucleo_total_table.C;
    case 'G':
        return &nucleo_total_table.G;
    case 'T':
        return &nucleo_total_table.T;
    default:
        return nullptr; // Return nullptr for invalid bases
    }
}
