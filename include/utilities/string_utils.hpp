#pragma once
#include <optional>
#include <string_view>
#include <string>
#include <algorithm>
#include <vector>
#include <concepts>
#include <type_traits>
#include <cctype>
#include <charconv>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <utility>

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
    if (str.empty())
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
    return !str.empty() && std::all_of(str.begin(), str.end(), [](unsigned char ch)
                                       { return std::isalpha(ch) != 0; });
}

enum class StringToIntStatus
{
    SUCCESS = 0,
    EMPTY_INPUT,
    INVALID_CHARACTER,
    OUT_OF_RANGE,
    TRAILING_CHARACTER
};

[[nodiscard]]
inline StringToIntStatus StringToInt(const std::string_view str, int &output) noexcept
{
    if (str.empty())
    {
        return StringToIntStatus::EMPTY_INPUT;
    }

    int result{};
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), result);
    switch (ec)
    {
    case std::errc::invalid_argument:
        return StringToIntStatus::INVALID_CHARACTER;
    case std::errc::result_out_of_range:
        return StringToIntStatus::OUT_OF_RANGE;
    default:
        break;
    }

    if (ptr != str.data() + str.size())
    {
        return StringToIntStatus::TRAILING_CHARACTER;
    }
    output = result;
    return StringToIntStatus::SUCCESS;
}

template <typename T, typename F>
using split_transform_result_t =
    std::invoke_result_t<F &, std::string_view, T &>;

template <typename T, typename F>
concept StatusSplitTransformer =
    std::default_initializable<T> &&
    std::move_constructible<T> &&
    std::invocable<F &, std::string_view, T &> &&
    std::equality_comparable<split_transform_result_t<T, F>>;

template <typename Status>
struct transformWhileSplitResult
{
    std::size_t field_index;
    std::size_t offset;
    Status status;
};

template <typename T, typename F>
    requires StatusSplitTransformer<T, F>
[[nodiscard]]
inline auto transformWhileSplit(std::string_view input, char delimiter, std::vector<T> &output,
                                F &&func, std::invoke_result_t<F &, std::string_view, T &> success_status)
    -> transformWhileSplitResult<split_transform_result_t<T, F>>
{
    using Status = split_transform_result_t<T, F>;

    std::vector<T> transformed;
    std::size_t field_begin = 0;
    std::size_t field_index = 0;

    for (std::size_t pos = 0; pos <= input.size(); ++pos)
    {
        if (pos != input.size() && input[pos] != delimiter)
        {
            continue;
        }
        const std::string_view field = input.substr(field_begin, pos - field_begin);

        T converted{};

        const Status status = std::invoke(func, field, converted);

        if (status != success_status)
        {
            return {
                field_index, field_begin, status};
        }

        transformed.push_back(std::move(converted));
        field_begin = pos + 1;
        field_index += 1;
    }
    output = std::move(transformed);

    return {
        field_index,
        input.size(),
        success_status};
}

/**
 * @brief split record line and save fields to the vector<string>
 *
 * @param[in] str record line
 * @return std::vector<std::string> vector of fields
 */
inline std::vector<std::string> split(const std::string &str)
{
    std::vector<std::string> fields;
    std::string field;
    char del = '\x09'; // aka Tab button or \t
    for (char ch : str)
    {
        if (ch == del)
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
