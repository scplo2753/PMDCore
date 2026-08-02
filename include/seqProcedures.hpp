#pragma once
#include "utility.hpp"
#include "parsedData.hpp"

/******
 * @todo: this function is not used yet, but it is a good idea to implement it for future use. It can be used to validate the data before processing.
 * @brief check the data is legal or not, if illegal, return false, otherwise return
 */
bool isLegalData(const recordLine_struct_t &data);

/**
 * @brief Check the record is legal or not, if legal, output the splited to data(for calculate PMD score), and splited_record
 *
 * @param[in] Line
 * @param[out] vector_splited_record
 * @param[out] data
 * @return true legal record
 * @return false bad record
 */
bool validAndParse(const std::string &Line,std::vector<std::string> vector_splited_record, recordLine_struct_t &data);

/**
 * @brief if MD list have nothing, then skip this line. alignment string shows where don't match the reference sequence
 * @param[in] &data: parsedData
 * @param[out] &alignmentData: the struct that need to output
 * @return -1 if the MD:Z is not exists, skip this line.
 *         0 like normal
 */
int ReconstructAlignmentAndRefSeq(parsedData &data, alignnmentData_t &alignnmentData);

/**
 * @brief Check whether the GC content of a reconstructed reference sequence is within configured bounds.
 *
 * @param[in] alignmentData The alignment data containing the reconstructed ref_seq.
 * @return true if the GC fraction is between FLAGS_minGC and FLAGS_maxGC, inclusive.
 * @return false if the ref_seq GC content is outside the allowed range.
 */
bool isGCcontentInRange(alignnmentData_t& alignmentData);

/****
* @brief parse the raw string to the struct of fields
* @related recordLine_struct_t
* @return return -1 if line is illegal, otherwise return 0
*****/
int splitOneLine(const std::string &line, std::vector<std::string> vector_splited_record, recordLine_struct_t &data);
