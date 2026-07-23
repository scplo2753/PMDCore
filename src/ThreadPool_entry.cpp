#include "ThreadPool_entry.hpp"
/**
 * @brief buffer the output line, and flush to stdout when the buffer size exceeds a threshold
 *
 * @param line the sam record line to be output
 */
void buffered_output_line(const std::string &line)
{
    // append to thread-local buffer
    tls_output_buffer.buffer.append(line);
    tls_output_buffer.buffer.push_back('\n');

    if (tls_output_buffer.buffer.size() >= OUTPUT_BUFFER_FLUSH_SIZE)
    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << tls_output_buffer.buffer;
        tls_output_buffer.buffer.clear();
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

    /***
     * @todo imple functions about add PMDS tag
     * @todo imple writesamfield
     * @todo imple flagss
     * @todo imple maskss
     * @todo printDS
     * @todo imple dry
     * @todo imple printalignments
     */
}
