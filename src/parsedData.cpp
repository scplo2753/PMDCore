#include "parsedData.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <charconv>
#include <limits>
#include <string>
#include "sam/struct_record.hpp"
#include "utilities/sequence_utils.hpp"

static bool isCigarOp(char c)
{
    switch (c)
    {
        case 'M':
        case 'I':
        case 'D':
        case 'N':
        case 'S':
        case 'H':
        case 'P':
        case '=':
        case 'X':
            return true;
        default:
            return false;
    }
}

parsedRecordResult parsedData::parseRawData(const recordLine_struct_t& raw_data)
{
    CIGARList_t cigar_list{};
    auto parse_status = parseCIGAR(raw_data,cigar_list);
    if(parse_status!=parsedRecordError::SUCCESS)
    {
        return parse_status;
    }

    if (factory_isTagNotExists(raw_data, "MD"))
    {
        std::cerr<<"Warning: "<<raw_data.QNAME<<" has no MD tag, this line will be skipped.\n";
        return parsedRecordError::MD_TAG_NOT_EXIST;
    }
    else if (factory_isTagValueEmpty(raw_data, "MD"))
    {
        std::cerr<<"Warning: "<<raw_data.QNAME<<" has an empty MD tag, this line will be skipped.\n";
        return parsedRecordError::MD_TAG_EMPTY;
    }

    return parsedData(raw_data,std::move(cigar_list));
}

void parsedData::set_ReadSeq_reverseSeq()
{
    if (!flag_isReadReversed)
    {
        data.read_seq = revcomp(data.read_seq);
        std::reverse(data.quality_scores.begin(), data.quality_scores.end());
    }
    flag_isReadReversed = true;
}

bool parsedData::isStatusFlagSet(SamFlag flag) const
{
    const auto flag_value = static_cast<int>(flag);
    try
    {
        return std::stoi(data.FLAG) & flag_value;
    }
    catch (const std::invalid_argument& e)
    {
        if (flag == SamFlag::READ_REVERSE_STRAND || data.FLAG.find('r') != std::string::npos)
        {
            return true;
        }
        else
            return false;
    }
}

bool parsedData::factory_isTagNotExists(const recordLine_struct_t &data, const std::string &tag)
{
    return !data.options_map.contains(tag);
}

bool parsedData::isTagExists(const std::string& tag) const
{
    if (data.options_map.find(tag) != data.options_map.end())
    {
        return true;
    }
    return false;
}

bool parsedData::factory_isTagValueEmpty(const recordLine_struct_t& raw_data, const std::string& tag)
{
    return raw_data.options_map.at(tag).value.empty();
}

std::string parsedData::getTagValue(const std::string& tag) const
{
    if (isTagExists(tag))
    {
        std::string values=data.options_map.at(tag).value;
        return values;
    }
    return "";
}

/******
* @brief return splited MD tag, The MD tag is split into numbers and non-numbers, for example, "10A5^C3" will be split into ["10", "A", "5", "^C", "3"].
* @returns return splitd MD tags if exists, else return empty vector. 
*******/
std::vector<std::string> parsedData::getMDList() const
{
    std::vector<std::string> md_list;

    const std::string& md_value = getMDTagValue();

        size_t pos = 0;
        const size_t md_value_size = md_value.size();
        while (pos < md_value_size)
        {
            char c = md_value[pos];
            if (c >= '0' && c <= '9')
            {
            const auto start = pos;
                while (pos < md_value_size && md_value[pos] >= '0' && md_value[pos] <= '9')
                {
                    ++pos;
                }
            md_list.emplace_back(md_value.substr(start, pos - start));
            }
            else if (c == '^')
            {
                size_t start = pos;
                ++pos;
                while (pos < md_value_size && md_value[pos] >= 'A' && md_value[pos] <= 'Z')
                {
                    ++pos;
                }
                md_list.emplace_back(std::string(md_value.substr(start, pos - start)));
            }
            else if (c >= 'A' && c <= 'Z')
            {
                md_list.emplace_back(1,c);
                ++pos;
            }
            else
                ++pos;
    }
    return md_list;
}

parsedRecordError parsedData::parseCIGAR(const recordLine_struct_t& raw_data, CIGARList_t& output)
{
    if (raw_data.cigar.empty() || raw_data.cigar == "*")
    {
        std::cerr <<"Warning: "<<raw_data.QNAME << " has no CIGAR data\n";
        return parsedRecordError::CIGAR_EMPTY;
    }

    std::vector<std::pair<char,std::size_t>> temp_cigar_list{};
    std::string_view cigar = raw_data.cigar;
    const size_t n = cigar.size();

    size_t pos = 0;

    while (pos < n)
    {
        size_t num_start = pos;
        while (pos < n && cigar[pos] >= '0' && cigar[pos] <= '9')
            ++pos;

        if (pos == num_start)
        {
            std::cerr << "Warning: " << raw_data.QNAME << " has an invalid CIGAR format. This line will be skipped\n";
            return parsedRecordError::BAD_CIGAR_FORMAT;
        }

        if (pos >= n)
        {
            std::cerr << "Warning: " << raw_data.QNAME << " has an invalid CIGAR format. This line will be skipped\n";
            return parsedRecordError::BAD_CIGAR_FORMAT;
        }

        char op = cigar[pos];

        if (!isCigarOp(op))
        {
            std::cerr << "Warning: " << raw_data.QNAME << " has an invalid CIGAR operation: " << op
                      << ". This line will be skipped." << "\n";
            return parsedRecordError::INVALID_CIGAR_OPERATION;
        }

        std::string_view str_cigar_step{cigar.data()+num_start,pos-num_start};
        std::size_t step{};
        auto [ptr, ec] = std::from_chars(str_cigar_step.data(), str_cigar_step.data() + str_cigar_step.size(), step);
        if(ec!=std::errc{}||ptr!=str_cigar_step.data() + str_cigar_step.size())
        {
            std::cerr << "Warning: " << raw_data.QNAME << " has an invalid CIGAR operation length " << str_cigar_step << "\n";
            return parsedRecordError::INVALID_CIGAR_STEP;
        }
        temp_cigar_list.emplace_back(op, step);
        ++pos;
    }

    output = std::move(temp_cigar_list);
    return parsedRecordError::SUCCESS;
}


std::vector<std::size_t> parsedData::getOpListInCIGAR(char Op)
{
    std::vector<std::size_t> Op_list{};
    std::size_t current_pos = 0;
    for (auto &line : cigar_list)
    {
        std::size_t step_length = line.second;

        if(step_length>std::numeric_limits<std::size_t>::max() - current_pos)
        {
            throw std::overflow_error("CIGAR operation length overflow: " + std::to_string(line.second));
        }

        if (line.first == Op)
        {
            for (std::size_t temp = current_pos; temp<(current_pos + step_length);temp++)
                Op_list.emplace_back(temp);
        }
        current_pos += step_length;
    }
    return Op_list;
}
