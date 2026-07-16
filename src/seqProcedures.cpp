#include "seqProcedures.hpp"
#include "arguments.hpp"

/******
 * @todo: this function is not used yet, but it is a good idea to implement it for future use. It can be used to validate the data before processing.
 * @brief check the data is legal or not, if illegal, return false, otherwise return
 */
bool isLegalData(const AlignLine_Data_t &data)
{
    if (data.QNAME.empty() || data.FLAG.empty() || data.RNAME.empty() || data.cigar.empty() || data.read_seq.empty() || data.quality_scores.empty())
    {
        std::cout << "Warning: Missing required fields in the input data. This line will be skipped." << std::endl;
        return false;
    }
    return true;
}

/****
* @brief parse the raw string to the struct of tokens
* @related AlignLine_data_t
* @return return -1 if line is illegal, otherwise return 0
*****/
int splitOneLine(const std::string &line, AlignLine_Data_t &data)
{
    std::vector<std::string> tokens(split(line));
    
    if (tokens.size() < 11)
        return BAD_RECORD_HAVE_EMPTY_FIELD;
    
    for(auto &token : tokens)
    {
        if(token.empty())
        {
            std::cout << "Warning: Empty field detected in the input data. This line will be skipped." << std::endl;
            return BAD_RECORD_HAVE_EMPTY_FIELD;
        }
    }

    //begin to parse the tokens into the struct
    data.QNAME= tokens[0];

    if (tokens[1].find_first_not_of("0123456789") != std::string::npos)
        return BAD_RECORD_ILLEGAL_FLAG;
    data.FLAG = tokens[1];

    data.RNAME = tokens[2];

    if (!isStringDigit(tokens[3])|| std::stoi(tokens[3]) < 0)
        return BAD_RECORD_ILLEGAL_POS;
    data.POS = std::stoi(tokens[3]);

    if (!isStringDigit(tokens[4]) || std::stoi(tokens[4]) < 0 || std::stoi(tokens[4]) > 255)
        return BAD_RECORD_ILLEGAL_MAPQ;
    data.MAPQ = std::stoi(tokens[4]);

    data.cigar = tokens[5];
    data.read_seq = tokens[9];
    data.quality_scores = tokens[10];
    for (size_t i = 11; i < tokens.size(); ++i)
    {
        const auto &field=tokens[i];
        if (field.size() < 5 || field[2] != ':' || field[4] != ':') //aviod illegal field format
        {
            continue;
        }
        data.options_map.insert_or_assign(field.substr(0, 2), options_value{field[3], field.substr(5)});
    }
    return 0;
}

/*
* @brief valid the line is legal, and split line into struct
*/
bool validAndParse(const std::string &Line,AlignLine_Data_t &data)
{
    if(Line.empty())
        return false;
    if(Line[0]=='@')
        return false;
    splitOneLine(Line, data);
    if(data.quality_scores.size()<2)
        return false;
    return true;
}

/*******
 * @brief if MD list have nothing, then skip this line. alignment string shows where don't match the reference sequence
 * @param[in] &data: parsedData
 * @param[out] &alignmentData: the struct that need to output
 * @return -1 if the MD:Z is not exists, skip this line
 * @return 0 like normal
 *******/
int ReconstructAlignmentAndRefSeq(parsedData &data, alignnmentData_t &alignnmentData)
{
    std::vector<std::string> md_list = data.getMDList();
    if(md_list.empty())
    {
        std::cout << "Warning: MD tag is missing. This line will be skipped." << std::endl;
        return -1;
    }
    int current_pos = 0;
    std::vector<uint> insertion_list = data.getInsertionList();
    std::vector<uint> softclip_list = data.getSoftClipList();
    const std::string &read = data.getReadSeq();
    for (auto &md : md_list)
    {
        if (std::isdigit(md[0]))
        {
            int match_length = std::stoi(md);
            alignnmentData.alignment += std::string(match_length, '.');
            alignnmentData.ref_seq += read.substr(current_pos, match_length);
            current_pos += match_length;
        }
        else if (md[0] == ('^'))
        {
            alignnmentData.alignment += md.substr(1);
            continue;
        }
        else // if only have alphabet
        {
            alignnmentData.alignment += md;
            alignnmentData.ref_seq += md;
            current_pos += md.length();
        }
    }

    if (!insertion_list.empty() || !softclip_list.empty())
    {
        int alignment_ptr = 0;
        alignnmentData.ref_seq.clear();

        for (uint x = 0; x < read.length(); x++)
        {
            bool is_insertion = std::binary_search(insertion_list.begin(), insertion_list.end(), x);
            bool is_softclip = std::binary_search(softclip_list.begin(), softclip_list.end(), x);

            if (is_insertion || is_softclip)
            {
                alignnmentData.ref_seq += '-';
            }
            else
            {
                if (alignment_ptr >= static_cast<int>(alignnmentData.alignment.size()))
                {
                    return -1;
                }

                char ref_base = alignnmentData.alignment[alignment_ptr];
                if (ref_base == '.')
                {
                    alignnmentData.ref_seq += read[x];
                }
                else
                {
                    alignnmentData.ref_seq += ref_base;
                }
                alignment_ptr++;
            }
        }
    }
    if (data.isReversedSeq())
    {
        data.set_ReadSeq_reverseSeq(); //data.seq=revcmp() and reverse quality_score
        alignnmentData.ref_seq = revcomp(alignnmentData.ref_seq);
    }
    return 0;
}

//options: get masked seq
std::string isMaskOptionExists(std::string& real_read)
{
    if (FLAGS_maskterminalbases || FLAGS_maskterminaldeaminations)
    {
        return real_read;
    }
    return std::string();
}

//step2
bool isGCcontentInRange(alignnmentData_t& alignmentData)
{
    int C_content = std::count(alignmentData.ref_seq.cbegin(), alignmentData.ref_seq.cend(), 'C');
    int G_content = std::count(alignmentData.ref_seq.cbegin(), alignmentData.ref_seq.cend(), 'G');
    float GC_content = 1.0 * (C_content + G_content) / alignmentData.ref_seq.length();
    if (GC_content > FLAGS_maxGC|| GC_content < FLAGS_minGC)
    {
        return false;
    }
    else
        return true;
}