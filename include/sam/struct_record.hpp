#pragma once
#include <string>
#include <cstdint>
#include <unordered_map>
#include "option_field.hpp"

using options_map_t = std::unordered_map<std::string,options_value>;

/**
* @todo review sam/bam optional field, if value can be multiple
*/
struct recordLine_struct_t
{
    std::string QNAME;
    std::string FLAG;
    std::string RNAME;
    std::uint32_t POS;
    std::uint16_t MAPQ;
    std::string cigar;
    std::string read_seq;
    std::string quality_scores;
    options_map_t options_map;
};