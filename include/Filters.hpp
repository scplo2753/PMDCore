#pragma once
#include <string_view>
#include "arguments.hpp"
#include "sam/struct_record.hpp"

inline bool lengthFilter(const recordLine_struct_t &data);

bool cigarFilter(const recordLine_struct_t &data);

bool ArgsFilter(const recordLine_struct_t &data);

//if G C T A not in ref_seq, skip this line
bool badRefSeq_Vailder(std::string_view real_ref_seq, std::string_view line);
