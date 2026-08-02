#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <thread>

#include "utility.hpp"
#include "arguments.hpp"
#include "parsedData.hpp"
#include "seqProcedures.hpp"
#include "Filters.hpp"
#include "calPMD.hpp"
#include "ThreadPool.hpp"
#include "argparse/argparse.hpp"
#include "ThreadPool_entry.hpp"
#include "utilities_wrappers.hpp"
#include "statics.hpp"

#define __DEBUG__ ///Enable debug module for compare result with origin program output
//#define __VERBOSE__

static void merge_match_dicts(match_dict_t &dst, const match_dict_t &src)
{
    for (const auto &entry : src)
    {
        dst[entry.first] += entry.second;
    }
}

static void merge_statics_dicts(statics_dicts_t &dst, const statics_dicts_t &src)
{
    merge_match_dicts(dst.match_dict, src.match_dict);
    merge_match_dicts(dst.match_dict_CpG, src.match_dict_CpG);
    merge_match_dicts(dst.match_dict_rev, src.match_dict_rev);
    merge_match_dicts(dst.match_dict_CpG_rev, src.match_dict_CpG_rev);
    merge_match_dicts(dst.mismatch_dict, src.mismatch_dict);
    merge_match_dicts(dst.mismatch_dict_CpG, src.mismatch_dict_CpG);
    merge_match_dicts(dst.mismatch_dict_rev, src.mismatch_dict_rev);
    merge_match_dicts(dst.mismatch_dict_CpG_rev, src.mismatch_dict_CpG_rev);
}

int main(int argc, char *argv[])
{
    initCMDParse(argc, argv);

    constexpr size_t BUFFER_SIZE = 1024 * 1024;
    char *buffer = new char[BUFFER_SIZE];
    std::cin.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    #ifdef __VERBOSE__
    std::cout << "Processing the input file..." << std::endl;
    #endif

    std::vector<double> ancient_model_deam = ancientModelDeam_wrapper();
    std::vector<double> modern_model_deam(1000, 0.001);

    // ============ multi thread config ============
    size_t thread_count = std::thread::hardware_concurrency();
    if (thread_count == 0)
        thread_count = 1;

    std::vector<statics_dicts_t> thread_statics(thread_count);
    std::vector<std::string> thread_output_buffers; // removed per-thread external buffers; kept empty for compatibility

    ThreadPool thread_pool(
        thread_count, [&](size_t index)
        { 
            tls_statics_dict = &thread_statics[index];
            tls_output_buffer.buffer.reserve(OUTPUT_BUFFER_FLUSH_SIZE); });
    std::cout << "Using " << thread_pool.get_thread_count() << " threads" << std::endl;
    // ============= end ======================

    // std::ios::sync_with_stdio(false);
    // std::cin.tie(nullptr);

    std::string line;
    std::vector<std::string> splited_record;
    size_t line_count = 0;
    size_t task_count = 0;

    while (std::getline(std::cin, line))
    {
        recordLine_struct_t raw_data;

        if (!validAndParse(line, splited_record, raw_data))
            continue;

        line_count += 1;

        // 数据预处理
        if (!ArgsFilter(raw_data))
            continue;

        // is line reversed
        bool isReverse=false;
        if ((std::stoi(raw_data.FLAG) & 16) == 1)
        {
            isReverse = true;
        }

        if (FLAGS_noreverse)
        {
            if (isReverse == true)
            {
                continue;
            }
        }

        // end
        
        double LR=0.0f;
        bool DSfield = false;
        if (raw_data.options_map.find("DS") != raw_data.options_map.end())
        {
            DSfield = true;
            std::string DSvalue = raw_data.options_map.at("DS").value;
            LR = std::stod(DSvalue);
        }

        parsedData data_ptr(raw_data);
        alignnmentData_t alignnmentData;
        if (DSfield==false||FLAGS_writesamfield || IS_USED_basic && FLAGS_basic > 0||FLAGS_terminal) // not complete yet
        {
            int result = ReconstructAlignmentAndRefSeq(data_ptr, alignnmentData);
            if (result == -1)
                continue;
        }

        string maskedseq{};
        if (ParamChecks::isUsing_maskTerminalBases()|| ParamChecks::isUsing_maskTerminalDeaminations())
            maskedseq = data_ptr.getReadSeq();
        if (!isGCcontentInRange(alignnmentData))
            continue;
        if (!badRefSeq_Vailder(alignnmentData.ref_seq, line))
            continue;
        /// @todo imple basecomposition param
        /// @todo imple basic param
        /// @todo imple terminal param
        if (IS_USED_basic && FLAGS_basic > 0)
        {
            if (function_basicFilter(alignnmentData.ref_seq, data_ptr.getReadSeq(), data_ptr.getReadSeq().size(), data_ptr.getQualityScores()))
            {
                std::cout << line << std::endl;
            }
        }
        if (FLAGS_terminal)
        {
            if (function_basicTerminal(data_ptr.getReadSeq(), alignnmentData.ref_seq, data_ptr.getQualityScores()))
            {
                std::cout << line << std::endl;
                continue;
            }
        }
        /// @todo imple first param
        /// @todo imple Leipzigsimple
        /// @todo imple customterminus
        /// @todo imple if options.perc_identity > 0.01 or options.printalignments:

        // 注意：需要复制 raw_data 和 alignnmentData，避免栈空间问题
        WorkItem work_item = {raw_data, data_ptr, alignnmentData};

        thread_pool.enqueue(
            process_single_line,
            work_item,
            line,
            std::ref(modern_model_deam),
            std::ref(ancient_model_deam),
            maskedseq,
            splited_record,
            isReverse);

        task_count++;
    }

    #ifdef __VERBOSE__
    std::cout << "Waiting for all tasks to complete..." << std::endl;
    thread_pool.wait(); // 等待所有线程完成

    std::cout << "Processed " << task_count << " lines from " << line_count << " total lines" << std::endl;
    #else
    thread_pool.wait(); //wait for all thread complete
    #endif

    statics_dicts_t merged_statics;

    for (auto &local_statics : thread_statics)
    {
        merge_statics_dicts(merged_statics, local_statics);
    }

    /**
     * @todo imple param first
     * @todo imple param stats
     * @todo imple param deamination
     * @todo imple param basecomposition
     */

    if (FLAGS_platypus)
    {
        statics(merged_statics);
    }

    return 0;
}
