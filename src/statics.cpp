#include <iostream>
#include "utility.hpp"
#include "arguments.hpp"
#include "statics.hpp"

using std::string;
void statics(statics_dicts_t &statics_dict)
{
    match_dict_t &mismatch_dict = statics_dict.mismatch_dict;
    match_dict_t &mismatch_dict_rev = statics_dict.mismatch_dict_rev;
    match_dict_t &mismatch_dict_CpG = statics_dict.mismatch_dict_CpG;
    match_dict_t &mismatch_dict_CpG_rev = statics_dict.mismatch_dict_CpG_rev;

    std::vector<string> pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};
    std::vector<string> CpG_pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};

    std::unordered_map<string, int> itotal_dict;
    std::unordered_map<string, int> CpG_itotal_dict;
    std::unordered_map<string, int> ztotal_dict;
    std::unordered_map<string, int> CpG_ztotal_dict;

    for (const char &key : "CTGA")
        for (int i = 0; i < FLAGS_range; i++)
        {
            int itotal = 0;
            int ztotal = 0;
            int cpg_itotal = 0;
            int cpg_ztotal = 0;
            for (const char &key2 : "CTGA")
            {
                string thekey;
                thekey.push_back(key);
                thekey.push_back(key2);
                thekey += std::to_string(i);
                if (mismatch_dict.find(thekey) != mismatch_dict.end())
                {
                    itotal += mismatch_dict[thekey];
                    if (mismatch_dict_CpG.find(thekey) != mismatch_dict_CpG.end())
                    {
                        cpg_itotal += mismatch_dict_CpG[thekey];
                    }
                }
                if (mismatch_dict_rev.find(thekey) != mismatch_dict_rev.end())
                {
                    ztotal += mismatch_dict_rev[thekey];
                    if (mismatch_dict_CpG_rev.find(thekey) != mismatch_dict_CpG_rev.end())
                    {
                        cpg_ztotal += mismatch_dict_CpG_rev[thekey];
                    }
                }
            }
            string dict_key = std::to_string(i) + key;
            itotal_dict[dict_key] = itotal;
            CpG_itotal_dict[dict_key] = cpg_itotal;
            ztotal_dict[dict_key] = ztotal;
            CpG_ztotal_dict[dict_key] = cpg_ztotal;
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
            thetotal = itotal_dict[std::to_string(i) + pair[0]];
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
            thetotal = ztotal_dict[std::to_string(i) + pair[0]];
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
            thetotal = CpG_itotal_dict[std::to_string(i) + pair[0]];
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
            thetotal = CpG_ztotal_dict[std::to_string(i) + pair[0]];
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