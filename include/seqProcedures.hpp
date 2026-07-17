#pragma once
#include "utility.hpp"
#include "parsedData.hpp"

bool isLegalData(const AlignLine_Data_t &data);

bool validAndParse(const std::string &Line,std::vector<std::string> vector_splited_record, AlignLine_Data_t &data);

int ReconstructAlignmentAndRefSeq(parsedData &data, alignnmentData_t &alignnmentData);

bool isGCcontentInRange(alignnmentData_t& alignmentData);

int splitOneLine(const std::string &line, std::vector<std::string> vector_splited_record, AlignLine_Data_t &data);
