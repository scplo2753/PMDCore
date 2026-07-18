#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <mutex>
#include "utilities_math.hpp"

enum badRecordType
{
    BAD_RECORD_HAVE_EMPTY_FIELD,
    BAD_RECORD_ILLEGAL_FLAG,
    BAD_RECORD_ILLEGAL_POS,
    BAD_RECORD_ILLEGAL_MAPQ
};

inline std::unordered_map<char, char> Nmap{{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}, {'Y', 'Y'}, {'R', 'R'}, {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'}, {'N', 'N'}, {'-', '-'}};

enum flagBits
{
    READ_PAIRED = 0x1,
    PROPER_PAIR = 0x2,
    READ_UNMAPPED = 0x4,
    MATE_UNMAPPED = 0x8,
    READ_REVERSE_STRAND = 0x10,
    MATE_REVERSE_STRAND = 0x20,
    FIRST_IN_PAIR = 0x40,
    SECOND_IN_PAIR = 0x80,
    NOT_PRIMARY_ALIGNMENT = 0x100,
    READ_FAILS_QUALITY_CHECKS = 0x200,
    DUPLICATE_READ = 0x400,
    SUPPLEMENTARY_ALIGNMENT = 0x800
};

struct options_value
{
    char type;
    std::string value;
};

using match_dict_t = std::unordered_map<std::string, int>;
using options_map_t = std::unordered_map<std::string,options_value>;

/*
* @todo review sam/bam optional field, if value can be multiple
*/
struct AlignLine_Data_t
{
    std::string QNAME;
    std::string FLAG;
    std::string RNAME;
    uint POS;
    uint MAPQ;
    std::string cigar;
    std::string read_seq;
    std::string quality_scores;
    options_map_t options_map;
};

/****** 
 * @brief struct of reconstructed reference and alignment data
 * @details first var is reference, the second one is read's alignment data
*/
struct alignnmentData_t
{
    std::string ref_seq; ///<reconstructed reference sequence
    std::string alignment; ///<read's alignment data
};

struct real_data_t
{
    std::string real_read;
    std::string real_ref_seq;
};

struct statics_dicts_t
{
    match_dict_t match_dict;
    match_dict_t match_dict_CpG;
    match_dict_t match_dict_rev;
    match_dict_t match_dict_CpG_rev;
    match_dict_t mismatch_dict;
    match_dict_t mismatch_dict_CpG;
    match_dict_t mismatch_dict_rev;
    match_dict_t mismatch_dict_CpG_rev;
    
    // 互斥锁保护所有字典的并发访问
    mutable std::mutex dict_mutex;
};

static std::string revcomp(const std::string &input)
{
    auto temp = input;
    std::reverse(temp.begin(), temp.end());
    std::string output={};
    for (auto i : temp)
    {
        output += Nmap.at(i);
    }
    return output;
}

/*****
 * @todo need to implement the function to get the reference sequence from fasta file
 */
static void fa_get(std::string ref_seq, std::string chromosome, int start, int end)
{
}

/**
 * @brief split record line and save fields to the vector<string>
 * 
 * @param[in] str record line
 * @return std::vector<std::string> vector of fields
 */
static std::vector<std::string> split(const std::string &str)
{
    std::vector<std::string> fields;
    std::string field;
    char del = '\x09'; //aka Tab button or \t
    for (char ch : str)
    {
        if (ch == del)
        {
            fields.push_back(field);
            field.clear();
        }
        else
        {
            field += ch;
        }
    }
    fields.push_back(field);
    return fields;
}

/**
 * @brief To check is the string only include + - . and number
 *
 * @param[in] str string that need to check
 * @return true if is legal number;
 * @return false
 */
static bool isStringDigit(const std::string &str)
{
    bool already_have_point = false;
    bool already_have_plusNeg = false;
    for (int i = 0; i < str.size(); ++i)
    {
        if (str[i] == '+' || str[i] == '-')
            if (i == 0)
                already_have_plusNeg = true;
            else
            {
                return false;
            }
        else if (str[i] == '.')
            if (already_have_point == false)
            {
                already_have_point = true;
            }
            else
            {
                return false;
            }
        else
        {
            if(!isdigit(str[i]))
                return false;
        }
    }
    return true;
}

/**
 * @brief To check is the string only include alphabet
 * 
 * @param[in] str string that need to check
 * @return true if only have alphabet
 * @return false if have digit
 */
static bool isStringAlphabet(const std::string &str)
{
    return std::all_of(str.begin(),str.end(),::isalpha);
}

static std::string strip(const std::string &str, const char &delimiter)
{
    size_t start = str.find_first_not_of(delimiter);
    size_t end = str.find_last_not_of(delimiter);
    if (start == std::string::npos || end == std::string::npos)
    {
        return "";
    }
    return str.substr(start, end - start + 1);
}

static std::string lstrip(const std::string &str, const char &delimiter)
{
    size_t start = str.find_first_not_of(delimiter);
    if (start == std::string::npos)
    {
        return "";
    }
    return str.substr(start);
}

static std::string rstrip(const std::string &str, const char &delimiter)
{
    size_t end = str.find_last_not_of(delimiter);
    if (end == std::string::npos)
    {
        return "";
    }
    return str.substr(0, end + 1);
}

static double L_match(int fpos, const std::vector<double> &fmodel, std::string &fqual, const double &polymorphism)
{
    double P_damage = fmodel[fpos];
    double P_error = phred2prob(fqual[fpos] - 33) / 3.0;
    const double &P_poly = polymorphism;
    double P_match = (1.0 - P_damage) * (1.0 - P_error) * (1.0 - P_poly) + (P_damage * P_error * (1.0 - P_poly)) + (P_error * P_poly * (1.0 - P_damage));
    return P_match;
}

static double L_mismatch(int fpos, const std::vector<double> &fmodel, std::string &fqual, const double &polymorphism)
{
    return 1.0f - L_match(fpos, fmodel, fqual, polymorphism);
}

static double L_match_ss(int fpos, int zpos, const std::vector<double> &fmodel,std::string_view fqual,const double &polymorphism)
{
    double P_damage = fmodel[fpos];
    double P_damage2 = fmodel[zpos];
    double P_error = phred2prob(fqual[fpos] - 33) / 3.0;
    const double &P_poly = polymorphism;
    double P_match = (1.0 - P_damage) * (1.0 - P_error) * (1.0 - P_poly) * (1.0 - P_damage2) + (P_damage * P_error * (1.0 - P_poly) * (1.0 - P_damage2)) + (P_damage2 * P_error * (1.0 - P_poly) * (1.0 - P_damage)) + (P_error * P_poly * (1.0 - P_damage2) * (1.0 - P_damage));
    return P_match;
}

static double L_mismatch_ss(int fpos, int zpos, const std::vector<double> &fmodel,std::string_view fqual,const double &polymorphism)
{
    return 1.0f - L_match_ss(fpos, zpos, fmodel, fqual, polymorphism);
}