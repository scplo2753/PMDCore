#include "parsedData.hpp"
#include <algorithm>

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

parsedData::parsedData(const AlignLine_Data_t data) : data(data), flag_isReadReversed(false)
{
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

bool parsedData::isStatusFlagSet(int flag) const
{
    try
    {
        return std::stoi(data.FLAG) & flag;
    }
    catch (const std::invalid_argument& e)
    {
        if (flag == 0x10 || data.FLAG.find('r') != std::string::npos)
        {
            return true;
        }
        else
            return false;
    }
}

bool parsedData::isTagExists(const std::string& tag) const
{
    if (data.options_map.find(tag) != data.options_map.end())
    {
        return true;
    }
    return false;
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
* @brief return splited MD tag, if the tag does not exist, return an empty vector. The MD tag is split into numbers and non-numbers, for example, "10A5^C3" will be split into ["10", "A", "5", "^C", "3"].
* @returns return splitd MD tags if exists, else return empty vector. 
*******/
std::vector<std::string> parsedData::getMDList() const
{
    std::vector<std::string> md_list;
    if(isTagExists("MD"))
    {
        const std::string &md_value = getMDTagValue();

        size_t pos = 0;
        const size_t md_value_size = md_value.size();
        while (pos < md_value_size)
        {
            char c = md_value[pos];
            if (c >= '0' && c <= '9')
            {
                int value = 0;
                while (pos < md_value_size && md_value[pos] >= '0' && md_value[pos] <= '9')
                {
                    value = value * 10 + (md_value[pos] - '0');
                    ++pos;
                }
                md_list.emplace_back(std::to_string(value));
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
    }
    return md_list;
}

/******
 * @brief return splited CIGAR tag, if the tag does not exist, return an empty vector. The CIGAR tag is split into numbers and non-numbers, for example, "10M5S" will be split into [("M", "10"), ("S", "5")].
 * @returns return splitd CIGAR tags if exists, else return empty vector.
 * @todo consider to delete return empty vector if the CIGAR tag does not exist, since the CIGAR tag is a required field in SAM format.
 */
std::vector<std::pair<char,std::string>> parsedData::getCIGARList()
{
    if (!cigar_list.empty())
    {
        return cigar_list;
    }

    const std::string &cigar = data.cigar;
    const size_t n = cigar.size();

    size_t pos = 0;

    while (pos < n)
    {
        size_t num_start = pos;
        while (pos < n && cigar[pos] >= '0' && cigar[pos] <= '9')
            ++pos;

        if (pos == num_start)
        {
            ++pos;
            continue; 
        }

        if (pos >= n)
            break; 

        char op = cigar[pos];

        if (!isCigarOp(op))
        {
            ++pos;
            continue; 
        }

        cigar_list.emplace_back(op, cigar.substr(num_start, pos - num_start));
        ++pos;
    }

    return cigar_list;
}

std::vector<uint> parsedData::getOpListInCIGAR(char Op)
{
    std::vector<uint> Op_list{};
    uint current_pos = 0;
    for (auto &line : getCIGARList())
    {
        uint step_long = std::stoi(line.second);
        if (line.first == Op)
        {
            for (int temp = current_pos; temp<(current_pos + step_long);temp++)
                Op_list.emplace_back(temp);
        }
        current_pos += std::stoi(line.second);
    }
    return Op_list;
}