#pragma once
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include "statistics/statistics_types.hpp"

using std::string;
using std::unordered_map;
using std::vector;
using platypus_result_struct = unordered_map<string, vector<double>>;

void split_key(const string &key, string &name, int &index);
void init_platypus_result_struct(platypus_result_struct &platypus_result, platypus_result_struct &denominator_result);
void statics(platypus_statics_dicts_t &platypus_statics_dict, platypus_result_struct &platypus_result, const platypus_denominator_table_t &denominator_table);
void statics(platypus_statics_dicts_t &platypus_statics_dict);

void print_statics_result(const platypus_result_struct &platypus_result);
