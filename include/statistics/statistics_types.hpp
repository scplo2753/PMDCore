#pragma once
#include <vector>
#include <unordered_map>
#include <string>
#include <mutex>

using match_dict_t = std::unordered_map<std::string, int>;

struct platypus_statics_dicts_t
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

struct statics_nucleo_total_table_t
{
    std::vector<double> C;
    std::vector<double> A;
    std::vector<double> G;
    std::vector<double> T;

    explicit statics_nucleo_total_table_t(size_t range = 0){
        A.assign(range, 0.0);
        C.assign(range, 0.0);
        G.assign(range, 0.0);
        T.assign(range, 0.0);
    }
};

struct platypus_denominator_table_t
{
    statics_nucleo_total_table_t forward;
    statics_nucleo_total_table_t reverse;
    statics_nucleo_total_table_t forward_CpG;
    statics_nucleo_total_table_t reverse_CpG;

    explicit platypus_denominator_table_t(size_t range=0):
        forward(range),reverse(range),
        forward_CpG(range),reverse_CpG(range)
    {}
};
