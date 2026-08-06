#include <iostream>
#include <ranges>
#include "arguments.hpp"
#include "statics.hpp"

std::vector<string> pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};
std::vector<string> CpG_pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};



/**
 * @brief Splits a key into name and index components
 * 
 * @param[in] key The input key to split
 * @param[out] name The resulting name component
 * @param[out] index The resulting index component
 */
void split_key(const string &key, string &name,int &index)
{
    name = key.substr(0, 2);
    auto result=std::from_chars(key.data()+2, key.data() + key.size(), index);
    if(result.ec != std::errc())
    {
        throw std::runtime_error("Failed to parse index from key: " + key);
    }
}

void statics(statics_dicts_t &statics_dict)
{
    match_dict_t &mismatch_dict = statics_dict.mismatch_dict;
    match_dict_t &mismatch_dict_rev = statics_dict.mismatch_dict_rev;
    match_dict_t &mismatch_dict_CpG = statics_dict.mismatch_dict_CpG;
    match_dict_t &mismatch_dict_CpG_rev = statics_dict.mismatch_dict_CpG_rev;

    unordered_map<string, int> forward_total_dict;
    unordered_map<string, int> CpG_forward_total_dict;
    unordered_map<string, int> reverse_total_dict;
    unordered_map<string, int> CpG_reverse_total_dict;

    for (const char &key : "CTGA")
        for (int i = 0; i < FLAGS_range; i++)
        {
            int forwardTotal = 0; //forward strand total count
            int reverseTotal = 0; //reverse strand total count
            int cpg_forwardTotal = 0;
            int cpg_reverseTotal = 0;
            for (const char &key2 : "CTGA")
            {
                string thekey;
                thekey.push_back(key);
                thekey.push_back(key2);
                thekey += std::to_string(i);
                if (mismatch_dict.find(thekey) != mismatch_dict.end())
                {
                    forwardTotal += mismatch_dict[thekey];
                }
                if (mismatch_dict_CpG.find(thekey) != mismatch_dict_CpG.end())
                {
                    cpg_forwardTotal += mismatch_dict_CpG[thekey];
                }
                if (mismatch_dict_rev.find(thekey) != mismatch_dict_rev.end())
                {
                    reverseTotal += mismatch_dict_rev[thekey];
                }
                if (mismatch_dict_CpG_rev.find(thekey) != mismatch_dict_CpG_rev.end())
                {
                    cpg_reverseTotal += mismatch_dict_CpG_rev[thekey];
                }
            }
            string dict_key = std::to_string(i) + key;
            forward_total_dict[dict_key] = forwardTotal;
            CpG_forward_total_dict[dict_key] = cpg_forwardTotal;
            reverse_total_dict[dict_key] = reverseTotal;
            CpG_reverse_total_dict[dict_key] = cpg_reverseTotal;
        }

    std::cout << "z\t ";
    for (const auto &pair : pairs)
    {
        std::cout << pair << "5\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "3\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "_CpG_5\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "_CpG_3\t ";
    }
    std::cout << std::endl;

    double thecount = 0.0f;
    double thetotal = 0.0f;
    double frac = 0.0f;

    for (int i = 0; i < FLAGS_range; i++)
    {
        std::cout << i << "\t ";
        for (const auto &pair : pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict.find(temp_key) != mismatch_dict.end())
            {
                thecount = mismatch_dict[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = forward_total_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            #ifdef __DEBUG__
            if (frac == 1.0) {
                std::cout << "1.0\t ";
            } else {
                std::cout << py_round(frac, 4) << "\t ";
            }
            #else
                std::cout << py_round(frac, 4) << "\t ";
            #endif
        }

        for (const auto &pair : pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_rev.find(temp_key) != mismatch_dict_rev.end())
            {
                thecount = mismatch_dict_rev[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = reverse_total_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            #ifdef __DEBUG__
            if (frac == 1.0) {
                std::cout << "1.0\t ";
            } else {
                std::cout << py_round(frac, 4) << "\t ";
            }
            #else
                std::cout << py_round(frac, 4) << "\t ";
            #endif
        }

        for (const auto &pair : CpG_pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_CpG.find(temp_key) != mismatch_dict_CpG.end())
            {
                thecount = mismatch_dict_CpG[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = CpG_forward_total_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            #ifdef __DEBUG__
            if (frac == 1.0){
                std::cout << "1.0\t ";
            }
            else{
                std::cout << py_round(frac, 4) << "\t ";
            }
            #else
            std::cout << py_round(frac, 4) << "\t ";
            #endif
        }

        for (const auto &pair : CpG_pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_CpG_rev.find(temp_key) != mismatch_dict_CpG_rev.end())
            {
                thecount = mismatch_dict_CpG_rev[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = CpG_reverse_total_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            #ifdef __DEBUG__
            if (frac == 1.0){
                std::cout << "1.0\t ";
            }
            else{
                std::cout << py_round(frac, 4) << "\t ";
            }
            #else
            std::cout << py_round(frac, 4) << "\t ";
            #endif
        }
        std::cout << std::endl;
    }
}

void init_platypus_result_struct(platypus_result_struct &platypus_result,platypus_result_struct &denominator_result)
{
    for (const auto &pair : pairs)
    {
        platypus_result[pair + "5"] = std::vector<double>(FLAGS_range, 0.0);
        denominator_result[pair + "5"] = std::vector<double>(FLAGS_range, 0.0);
        platypus_result[pair + "3"] = std::vector<double>(FLAGS_range, 0.0);
        denominator_result[pair + "3"] = std::vector<double>(FLAGS_range, 0.0);
    }
    for (const auto &pair : CpG_pairs)
    {
        platypus_result[pair + "_CpG_5"] = std::vector<double>(FLAGS_range, 0.0);
        platypus_result[pair + "_CpG_3"] = std::vector<double>(FLAGS_range, 0.0);
        denominator_result[pair + "_CpG_5"] = std::vector<double>(FLAGS_range, 0.0);
        denominator_result[pair + "_CpG_3"] = std::vector<double>(FLAGS_range, 0.0);
    }
}

static const std::vector<double>* get_denominator_column(const std::string &result_key, const statics_denominator_table_t &denominator_table)
{
    const char ref_base = result_key[0];
    const bool is_CpG = result_key.find("_CpG_") != std::string::npos;
    const bool is_5 = result_key.back() == '5';

    const statics_nucleo_total_table_t *table = nullptr;
    if (is_CpG)
    {
        table = is_5 ? &denominator_table.forward_CpG : &denominator_table.reverse_CpG;
    }
    else
    {
        table = is_5 ? &denominator_table.forward : &denominator_table.reverse;
    }

    if (!table)
        return nullptr;

    switch (ref_base)
    {
    case 'A':
        return &table->A;
    case 'C':
        return &table->C;
    case 'G':
        return &table->G;
    case 'T':
        return &table->T;
    default:
        return nullptr;
    }
}

void statics(statics_dicts_t &statics_dict, platypus_result_struct &platypus_result, const statics_denominator_table_t &denominator_table)
{
    string name;
    int index;
    for (const auto &[key, value] : statics_dict.mismatch_dict)
    {
        split_key(key, name, index);
        if(index>=FLAGS_range)
        {
            continue;
        }
        const string result_key = name + "5";
        platypus_result.at(result_key).at(static_cast<int>(index)) = static_cast<double>(value);
    }
    for(const auto& [key,value]:statics_dict.mismatch_dict_rev)
    {
        split_key(key, name, index);
        if(index>=FLAGS_range)
        {
            continue;
        }
        const string result_key = name + "3";
        platypus_result.at(result_key).at(static_cast<int>(index)) = static_cast<double>(value);
    }
    for(const auto& [key,value]:statics_dict.mismatch_dict_CpG)
    {
        split_key(key, name, index);
        if(index>=FLAGS_range)
        {
            continue;
        }
        const string result_key = name + "_CpG_5";
        platypus_result.at(result_key).at(static_cast<int>(index)) = static_cast<double>(value);
    }
    for(const auto& [key,value]:statics_dict.mismatch_dict_CpG_rev)
    {
        split_key(key, name, index);
        if(index>=FLAGS_range)
        {
            continue;
        }
        const string result_key = name + "_CpG_3";
        platypus_result.at(result_key).at(static_cast<int>(index)) = static_cast<double>(value);
    }

    for (auto &[key, result_column] : platypus_result)
    {
        const std::vector<double> *denominator_column = get_denominator_column(key, denominator_table);
        if (!denominator_column)
            continue;

        std::transform(
            result_column.begin(),
            result_column.end(),
            denominator_column->begin(),
            result_column.begin(),
            [](double count, double denominator)
            {
                const double safe_denominator = (denominator == 0.0) ? 1.0 : denominator; // Avoid division by zero
                return denominator == 0.0 ? 0.0 : count / safe_denominator;
            });
    }
}

void print_statics_result(const platypus_result_struct &platypus_result)
{
    constexpr std::array output_header = {"CT5","CA5","CG5","CC5","GA5","GT5", "GC5", "GG5", "AA5", "AT5", "AC5", "AG5", "TA5", "TT5", "TC5", "TG5", "CT3", "CA3", "CG3", "CC3", "GA3", "GT3", "GC3", "GG3", "AA3", "AT3", "AC3", "AG3", "TA3", "TT3", "TC3", "TG3", "CT_CpG_5", "CA_CpG_5", "CG_CpG_5", "CC_CpG_5", "GA_CpG_5", "GT_CpG_5", "GC_CpG_5", "GG_CpG_5", "AA_CpG_5", "AT_CpG_5", "AC_CpG_5", "AG_CpG_5", "TA_CpG_5", "TT_CpG_5", "TC_CpG_5", "TG_CpG_5", "CT_CpG_3", "CA_CpG_3", "CG_CpG_3", "CC_CpG_3", "GA_CpG_3","GT_CpG_3", "GC_CpG_3", "GG_CpG_3", "AA_CpG_3", "AT_CpG_3", "AC_CpG_3", "AG_CpG_3", "TA_CpG_3", "TT_CpG_3", "TC_CpG_3", "TG_CpG_3"};
    std::cout << "z\t ";
    for (const auto &header : output_header)
    {
        std::cout << header << "\t ";
    }
    std::cout << '\n';

    std::array<const std::vector<double> *, output_header.size()> output_columns;
    std::transform(
        output_header.begin(),
        output_header.end(),
        output_columns.begin(),
        [&platypus_result](const auto &header)
        {
            return &platypus_result.at(header);
        });

    for (int i = 0; i < FLAGS_range; i++)
    {
        std::cout << i << "\t ";
        for (const auto *column : output_columns)
        {
            const long double rounded_value = py_round(column->at(i), 4);
            if (rounded_value == 1.0L)
            {
                std::cout << "1.0\t ";
            }
            else
            {
                std::cout << rounded_value << "\t ";
            }
        }
        std::cout << '\n';
    }
}
