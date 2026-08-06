#pragma once
#include <map>
#include <string>
#include <vector>
#include "statics_types.hpp"
#include "utility.hpp"

using std::string;
using std::unordered_map;
using std::vector;
using platypus_result_struct = unordered_map<string, vector<double>>;

void split_key(const string &key, string &name, int &index);
void init_platypus_result_struct(platypus_result_struct &platypus_result, platypus_result_struct &denominator_result);
void statics(statics_dicts_t &statics_dict, platypus_result_struct &platypus_result, const statics_denominator_table_t &denominator_table);
void statics(statics_dicts_t &statics_dict);

void print_statics_result(const platypus_result_struct &platypus_result);