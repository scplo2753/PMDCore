#pragma once
#include <cassert>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
#include "arguments.hpp"
#include "parsedData.hpp"
#include "calPMD.hpp"
#include "alignment/alignment.hpp"
#include "sam/struct_record.hpp"
#include "statistics/statistics_types.hpp"

inline std::mutex output_mutex;
inline constexpr size_t OUTPUT_BUFFER_FLUSH_SIZE = 1 << 20;

inline thread_local platypus_statics_dicts_t *tls_statics_dict = nullptr;
inline thread_local platypus_denominator_table_t *tls_platypus_denominator_table = nullptr;
inline thread_local deamination_statics_t *tls_deamination_statics_table = nullptr;

// 数据处理结构体
struct WorkItem
{
    recordLine_struct_t raw_data;
    parsedData parsed_data;
    alignnmentData_t alignment_data;
};

struct TLS_output_buffer_holder_struct{
    std::string buffer;
    ~TLS_output_buffer_holder_struct(){
        if(!buffer.empty())
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << buffer;
            std::cout.flush();
        }
    }
};

inline thread_local TLS_output_buffer_holder_struct tls_output_buffer;

void buffered_output_line(const std::string &line);

void process_single_line(
    const WorkItem &work_item,
    const std::string &line,
    const std::vector<double> &modern_model,
    const std::vector<double> &ancient_model,
    std::string maskedseq,
    const std::vector<std::string> &split_record,
    bool is_reverse);
