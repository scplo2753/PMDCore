#pragma once
#include <vector>

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

struct statics_denominator_table_t
{
    statics_nucleo_total_table_t forward;
    statics_nucleo_total_table_t reverse;
    statics_nucleo_total_table_t forward_CpG;
    statics_nucleo_total_table_t reverse_CpG;

    explicit statics_denominator_table_t(size_t range=0):
        forward(range),reverse(range),
        forward_CpG(range),reverse_CpG(range)
    {}
};