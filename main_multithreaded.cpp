#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <thread>

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
#include "statistics/statistics_types.hpp"
#include "statistics/deam_types.hpp"
#include "deamination_statics.hpp"
#include "pmd/deamination_model.hpp"

//#define __DEBUG__ ///Enable debug module for compare result with origin program output
//#define __VERBOSE__

static void merge_match_dicts(match_dict_t &dst, const match_dict_t &src)
{
    for (const auto &entry : src)
    {
        dst[entry.first] += entry.second;
    }
}

static void merge_statics_dicts(platypus_statics_dicts_t &dst, const platypus_statics_dicts_t &src)
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

static void merge_denominator_tables(platypus_denominator_table_t &dst, const platypus_denominator_table_t &src)
{
    for (size_t i = 0; i < FLAGS_range; ++i)
    {
        dst.forward.C[i] += src.forward.C[i];
        dst.forward.A[i] += src.forward.A[i];
        dst.forward.G[i] += src.forward.G[i];
        dst.forward.T[i] += src.forward.T[i];

        dst.reverse.C[i] += src.reverse.C[i];
        dst.reverse.A[i] += src.reverse.A[i];
        dst.reverse.G[i] += src.reverse.G[i];
        dst.reverse.T[i] += src.reverse.T[i];

        dst.forward_CpG.C[i] += src.forward_CpG.C[i];
        dst.forward_CpG.A[i] += src.forward_CpG.A[i];
        dst.forward_CpG.G[i] += src.forward_CpG.G[i];
        dst.forward_CpG.T[i] += src.forward_CpG.T[i];

        dst.reverse_CpG.C[i] += src.reverse_CpG.C[i];
        dst.reverse_CpG.A[i] += src.reverse_CpG.A[i];
        dst.reverse_CpG.G[i] += src.reverse_CpG.G[i];
        dst.reverse_CpG.T[i] += src.reverse_CpG.T[i];
    }
}

int main(int argc, char *argv[])
{
    initCMDParse(argc, argv);
    size_t range = FLAGS_range;

    constexpr size_t BUFFER_SIZE = 1024 * 1024;
    char *buffer = new char[BUFFER_SIZE];
    std::cin.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    #ifdef __VERBOSE__
    std::cout << "Processing the input file..." << std::endl;
    #endif

    // ======Initialize the ancient and modern deamination models=======
    std::vector<double> ancient_model_deam = ancientModelDeam_wrapper();
    std::vector<double> modern_model_deam(1000, 0.001);
    //==========================end=========================

    // ============ multi thread config ============
    size_t thread_count = std::thread::hardware_concurrency();
    if (thread_count == 0)
        thread_count = 1;

    std::vector<platypus_statics_dicts_t> thread_statics(thread_count);
    std::vector<std::string> thread_output_buffers; // removed per-thread external buffers; kept empty for compatibility
    std::vector<platypus_denominator_table_t> platypus_denominator_tables(thread_count, platypus_denominator_table_t(range));
    std::vector<deamination_statics_t> deamination_statics_tables(thread_count, deamination_statics_t(range));

    ThreadPool thread_pool(
        thread_count, [&](size_t index)
        { 
            tls_statics_dict = &thread_statics[index];
            tls_platypus_denominator_table = &platypus_denominator_tables[index];
            tls_deamination_statics_table = &deamination_statics_tables[index];
            tls_output_buffer.buffer.reserve(OUTPUT_BUFFER_FLUSH_SIZE); 
        });

    #ifdef __VERBOSE__
    std::cout << "Using " << thread_pool.get_thread_count() << " threads" << std::endl;
    #endif
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

        // data pre-process
        if (!ArgsFilter(raw_data))
            continue;

        // is line reversed
        bool isReverse=false;
        if ((std::stoi(raw_data.FLAG) & 16) != 0)
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

        //Temporarily disable reuse of DS:Z in optional field,
        //It would skip reference-sequence reconstruction, so the reference sequence var is empty.
        /**
        double LR=0.0f;
        bool DSfield = false;
        if (raw_data.options_map.find("DS") != raw_data.options_map.end())
        {
            DSfield = true;
            std::string DSvalue = raw_data.options_map.at("DS").value;
            LR = std::stod(DSvalue);
        }
        */

        parsedData data_ptr(raw_data);
        alignnmentData_t alignnmentData;
        int result = ReconstructAlignmentAndRefSeq(data_ptr, alignnmentData);
        if (result == -1)
            continue;

        string maskedseq{};
        if (ParamChecks::isUsing_maskTerminalBases()|| ParamChecks::isUsing_maskTerminalDeaminations())
            maskedseq = data_ptr.getReadSeq();
        if (!isGCcontentInRange(alignnmentData))
            continue;
        if (!badRefSeq_Vailder(alignnmentData.ref_seq, line))
            continue;
        /// @todo imple basecomposition param
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



    /**
     * @todo imple param first
     * @todo imple param stats
     * @todo imple param deamination
     * @todo imple param basecomposition
     */

    if(FLAGS_deamination)
    {
        deamination_statics_t merged_deamination_statistics(range);
        for (const auto &local_statistics : deamination_statics_tables)
        {
            merged_deamination_statistics += local_statistics;
        }

        print_deamination_statistics(merged_deamination_statistics);
    }

    if (FLAGS_platypus)
    {
        platypus_statics_dicts_t merged_statics;
        platypus_denominator_table_t merged_denominator_table(range);

        for (auto &local_statics : thread_statics)
        {
            merge_statics_dicts(merged_statics, local_statics);
        }
        for (auto &local_denominator_table : platypus_denominator_tables)
        {
            merge_denominator_tables(merged_denominator_table, local_denominator_table);
        }

        platypus_result_struct platypus_result;
        platypus_result_struct denominator_result;

        init_platypus_result_struct(platypus_result, denominator_result);
        statics(merged_statics, platypus_result, merged_denominator_table);
        print_statics_result(platypus_result);
    }

    return 0;
}
