#pragma once

#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include "parsedData.hpp"

// Use inside test bodies, not in global initialization. GoogleTest reports an
// unexpected factory error as a failed test instead of terminating at startup.
inline parsedData requireParsedRecord(const recordLine_struct_t& record)
{
    auto result = parsedData::parseRawData(record);
    if (const auto* error = std::get_if<parsedRecordError>(&result))
    {
        throw std::runtime_error(
            "Failed to create test record " + record.QNAME +
            " (CIGAR=" + record.cigar + ", status=" +
            std::to_string(static_cast<int>(*error)) + ")");
    }
    return std::get<parsedData>(std::move(result));
}
