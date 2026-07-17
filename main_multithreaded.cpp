
#include <cassert>
#include <iostream>
#include <fstream>
#include <iterator>
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

static thread_local statics_dicts_t *tls_statics_dict = nullptr;
static thread_local std::string *tls_output_buffer = nullptr;

static std::mutex output_mutex;
static constexpr size_t OUTPUT_BUFFER_FLUSH_SIZE = 1 << 20;

// 数据处理结构体
struct WorkItem
{
    AlignLine_Data_t raw_data;
    parsedData parsed_data;
    alignnmentData_t alignment_data;
};

std::vector<double> get_ancient_model_deam()
{
    std::vector<double> ancient_model_deam(1000);
    for (int pos = 1; pos < 1000; ++pos)
    {
        ancient_model_deam.at(pos - 1) = geometric(FLAGS_PMDpparam, pos, FLAGS_PMDconstant);
    }
    return ancient_model_deam;
}

void statics(statics_dicts_t &statics_dict)
{
    match_dict_t &mismatch_dict = statics_dict.mismatch_dict;
    match_dict_t &mismatch_dict_rev = statics_dict.mismatch_dict_rev;
    match_dict_t &mismatch_dict_CpG = statics_dict.mismatch_dict_CpG;
    match_dict_t &mismatch_dict_CpG_rev = statics_dict.mismatch_dict_CpG_rev;

    std::vector<string> pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};
    std::vector<string> CpG_pairs = {"CT", "CA", "CG", "CC", "GA", "GT", "GC", "GG", "AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG"};

    std::unordered_map<string, int> itotal_dict;
    std::unordered_map<string, int> CpG_itotal_dict;
    std::unordered_map<string, int> ztotal_dict;
    std::unordered_map<string, int> CpG_ztotal_dict;

    for (const char &key : "CTGA")
        for (int i = 0; i < FLAGS_range; i++)
        {
            int itotal = 0;
            int ztotal = 0;
            int cpg_itotal = 0;
            int cpg_ztotal = 0;
            for (const char &key2 : "CTGA")
            {
                string thekey;
                thekey.push_back(key);
                thekey.push_back(key2);
                thekey += std::to_string(i);
                if (mismatch_dict.find(thekey) != mismatch_dict.end())
                {
                    itotal += mismatch_dict[thekey];
                    if (mismatch_dict_CpG.find(thekey) != mismatch_dict_CpG.end())
                    {
                        cpg_itotal += mismatch_dict_CpG[thekey];
                    }
                }
                if (mismatch_dict_rev.find(thekey) != mismatch_dict_rev.end())
                {
                    ztotal += mismatch_dict_rev[thekey];
                    if (mismatch_dict_CpG_rev.find(thekey) != mismatch_dict_CpG_rev.end())
                    {
                        cpg_ztotal += mismatch_dict_CpG_rev[thekey];
                    }
                }
            }
            string dict_key = std::to_string(i) + key;
            itotal_dict[dict_key] = itotal;
            CpG_itotal_dict[dict_key] = cpg_itotal;
            ztotal_dict[dict_key] = ztotal;
            CpG_ztotal_dict[dict_key] = cpg_ztotal;
        }

    std::cout << "z\t ";
    for (const auto &pair : pairs)
    {
        std::cout << pair << "5\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "3\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "_CpG_5\t ";
    }
    for (const auto &pair : pairs)
    {
        std::cout << pair << "_CpG_3\t ";
    }
    std::cout << std::endl;

    double thecount = 0.0f;
    double thetotal = 0.0f;
    double frac = 0.0f;

    for (int i = 0; i < FLAGS_range; i++)
    {
        std::cout << i << "\t ";
        for (const auto &pair : pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict.find(temp_key) != mismatch_dict.end())
            {
                thecount = mismatch_dict[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = itotal_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            std::cout << py_round(frac, 4) << "\t ";
        }

        for (const auto &pair : pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_rev.find(temp_key) != mismatch_dict_rev.end())
            {
                thecount = mismatch_dict_rev[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = ztotal_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            std::cout << py_round(frac, 4) << "\t ";
        }

        for (const auto &pair : CpG_pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_CpG.find(temp_key) != mismatch_dict_CpG.end())
            {
                thecount = mismatch_dict_CpG[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = CpG_itotal_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            std::cout << py_round(frac, 4) << "\t ";
        }

        for (const auto &pair : CpG_pairs)
        {
            std::string temp_key = pair + std::to_string(i);
            if (mismatch_dict_CpG_rev.find(temp_key) != mismatch_dict_CpG_rev.end())
            {
                thecount = mismatch_dict_CpG_rev[temp_key];
            }
            else
            {
                std::cout << "0\t ";
                continue;
            }
            thetotal = CpG_ztotal_dict[std::to_string(i) + pair[0]];
            frac = 1.0 * thecount / thetotal;
            std::cout << py_round(frac, 4) << "\t ";
        }
        std::cout << std::endl;
    }
}

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

static void buffered_output_line(const std::string &line)
{
    assert(tls_output_buffer != nullptr);

    tls_output_buffer->append(line);
    tls_output_buffer->push_back('\n');

    if (tls_output_buffer->size() >= OUTPUT_BUFFER_FLUSH_SIZE)
    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << *tls_output_buffer;
        tls_output_buffer->clear();
    }
}

/**
 * @brief entry of the thread pool, entry function is calPMD
 * 
 * @param work_item struct of Read and reference sequence
 * @param line raw record line
 * @param modern_model modern deamination model
 * @param ancient_model ancient deamination model
 */
void process_single_line(
    const WorkItem &work_item,
    const std::string &line,
    const std::vector<double> &modern_model,
    const std::vector<double> &ancient_model,
    std::string maskedseq,
    const std::vector<std::string> &split_record,
    bool is_reverse)
{
    assert(tls_statics_dict != nullptr);
    calPMD calPMD_instance(
        real_data_t{work_item.parsed_data.getReadSeq(), work_item.alignment_data.ref_seq},
        modern_model,
        ancient_model,
        work_item.parsed_data.getQualityScores(),
        maskedseq,
        *tls_statics_dict);

    std::string output_line = line;
    if (IS_USED_maskterminaldeaminations || IS_USED_maskterminalbases)
    {
        std::string real_read = maskedseq;
        function_in_thread_pool_maskterminaldeam_or_maskterminalbases(maskedseq, real_read, is_reverse, split_record, output_line);
    }

    if (IS_USED_threshold)
    {
        if (calPMD_instance.threshold_filter())
        {
            buffered_output_line(output_line);
        }
    }
}



int main(int argc, char *argv[])
{
    initCMDParse(argc, argv);

    constexpr size_t BUFFER_SIZE = 1024 * 1024;
    char *buffer = new char[BUFFER_SIZE];
    std::cin.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    std::cout << "Processing the input file..." << std::endl;

    std::vector<double> ancient_model_deam = get_ancient_model_deam();
    std::vector<double> modern_model_deam(1000, 0.001);

    // ============ 多线程处理 ============
    size_t thread_count = std::thread::hardware_concurrency();
    if (thread_count == 0)
        thread_count = 1;

    std::vector<statics_dicts_t> thread_statics(thread_count);
    std::vector<std::string> thread_output_buffers(thread_count);

    ThreadPool thread_pool(
        thread_count, [&](size_t index)
        { 
            tls_statics_dict = &thread_statics[index]; 
            tls_output_buffer=&thread_output_buffers[index];
            tls_output_buffer->reserve(OUTPUT_BUFFER_FLUSH_SIZE); });
    std::cout << "Using " << thread_pool.get_thread_count() << " threads" << std::endl;

    // std::ios::sync_with_stdio(false);
    // std::cin.tie(nullptr);

    std::string line;
    std::vector<std::string> splited_record;
    size_t line_count = 0;
    size_t task_count = 0;

    while (std::getline(std::cin, line))
    {
        AlignLine_Data_t raw_data;

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

        /// @todo imple if(options.maskterminaldeaminations or options.maskterminalbases)
        string maskedseq{};
        if (IS_USED_maskterminalbases || IS_USED_maskterminaldeaminations)
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

    std::cout << "Waiting for all tasks to complete..." << std::endl;
    thread_pool.wait(); // 等待所有线程完FLAGS_basic成

    std::cout << "Processed " << task_count << " lines from " << line_count << " total lines" << std::endl;

    statics_dicts_t merged_statics;

    for (auto &local_statics : thread_statics)
    {
        merge_statics_dicts(merged_statics, local_statics);
    }

    if (FLAGS_platypus)
    {
        statics(merged_statics);
    }

    return 0;
}
