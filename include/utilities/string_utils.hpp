#pragma once
#include <charconv>
#include <optional>
#include <string_view>
#include <string>
#include <algorithm>
#include <vector>
#include <cctype>

inline std::string strip(const std::string &str, const char &delimiter)
{
    size_t start = str.find_first_not_of(delimiter);
    size_t end = str.find_last_not_of(delimiter);
    if (start == std::string::npos || end == std::string::npos)
    {
        return "";
    }
    return str.substr(start, end - start + 1);
}

inline std::string lstrip(const std::string &str, const char &delimiter)
{
    size_t start = str.find_first_not_of(delimiter);
    if (start == std::string::npos)
    {
        return "";
    }
    return str.substr(start);
}

inline std::string rstrip(const std::string &str, const char &delimiter)
{
    size_t end = str.find_last_not_of(delimiter);
    if (end == std::string::npos)
    {
        return "";
    }
    return str.substr(0, end + 1);
}

inline std::optional<unsigned int> parse_unsigned_integer(std::string_view str)
{
    if(str.empty())
    {
        return std::nullopt;
    }

    unsigned int value{};
    const char *first = str.data();
    const char *last = first + str.size();
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc{} || ptr != last)
    {
        return std::nullopt;
    }
    return value;
}

/**
 * @brief To check is the string only include alphabet
 *
 * @param[in] str string that need to check
 * @return true if only have alphabet
 * @return false if have digit
 */
inline bool isStringAlphabet(std::string_view str)
{
    return !str.empty() && std::all_of(str.begin(), str.end(), [](unsigned char ch) {
        return std::isalpha(ch) != 0;
    });
}

/**
 * @brief Split a string using a single-character delimiter.
 *
 * Empty fields, including trailing empty fields, are preserved.
 *
 * @param[in] str String to split.
 * @param[in] delimiter Field delimiter.
 * @return Split fields.
 */
inline std::vector<std::string> split(const std::string &str,const char delimiter)
{
    std::vector<std::string> fields;
    std::string field;
    for (char ch : str)
    {
        if (ch == delimiter)
        {
            fields.push_back(field);
            field.clear();
        }
        else
        {
            field += ch;
        }
    }
    fields.push_back(field);
    return fields;
}

inline std::vector<std::string> splitTable(const std::string &str)
{
    return split(str, '\x09');
}
