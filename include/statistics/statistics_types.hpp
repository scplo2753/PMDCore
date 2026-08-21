#pragma once
#include <vector>
#include <stdexcept>
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


    platypus_statics_dicts_t &operator+=(const platypus_statics_dicts_t &right)
    {
        merge_match_dicts(match_dict, right.match_dict);
        merge_match_dicts(match_dict_CpG, right.match_dict_CpG);
        merge_match_dicts(match_dict_rev, right.match_dict_rev);
        merge_match_dicts(match_dict_CpG_rev, right.match_dict_CpG_rev);
        merge_match_dicts(mismatch_dict, right.mismatch_dict);
        merge_match_dicts(mismatch_dict_CpG, right.mismatch_dict_CpG);
        merge_match_dicts(mismatch_dict_rev, right.mismatch_dict_rev);
        merge_match_dicts(mismatch_dict_CpG_rev, right.mismatch_dict_CpG_rev);

        return *this;
    }

    private:
        static void merge_match_dicts(match_dict_t &dest, const match_dict_t &src)
        {
            for (const auto &[key, count] : src)
            {
                dest[key] += count;
            }
        }
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

    [[nodiscard]] bool is_compatible_with(const statics_nucleo_total_table_t &right) const noexcept
    {
        const bool left_valid = A.size() == C.size() &&
                                A.size() == G.size() &&
                                A.size() == T.size();

        const bool right_valid = right.A.size() == right.C.size() &&
                                 right.A.size() == right.G.size() &&
                                 right.A.size() == right.T.size();

        return left_valid && right_valid && A.size() == right.A.size();
    }

    statics_nucleo_total_table_t &operator+=(const statics_nucleo_total_table_t &right)
    {
        if (!is_compatible_with(right))
        {
            throw std::invalid_argument("cannot merge nucleo tables with different ranges");
        }

        for (std::size_t row = 0; row < A.size(); ++row)
        {
            A[row] += right.A[row];
            C[row] += right.C[row];
            G[row] += right.G[row];
            T[row] += right.T[row];
        }
        return *this;
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

    platypus_denominator_table_t &operator+=(const platypus_denominator_table_t &right)
    {
        if (!forward.is_compatible_with(right.forward) ||
            !reverse.is_compatible_with(right.reverse) ||
            !forward_CpG.is_compatible_with(right.forward_CpG) ||
            !reverse_CpG.is_compatible_with(right.reverse_CpG))
        {
            throw std::invalid_argument("cannot merge denominator tables with different ranges");
        }

        forward += right.forward;
        reverse += right.reverse;
        forward_CpG += right.forward_CpG;
        reverse_CpG += right.reverse_CpG;

        return *this;
    }
};
