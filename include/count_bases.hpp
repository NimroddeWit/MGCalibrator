#ifndef COUNT_BASES_HPP
#define COUNT_BASES_HPP

#include <array>
#include <string>
#include <cstdint>
#include <filesystem>
#include <unordered_map>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/views/zip.hpp>

// ===================================
// 🔹 Function Prototypes
// ===================================

// ✅ Count total number of bases in single files
std::uint64_t count_total_bases_SE(const std::string& filePath);
std::uint64_t count_total_bases_PE(const std::string& file1, const std::string& file2);

// ✅ Batch processing for total base count from a folder
std::unordered_map<std::string, std::uint64_t> count_total_bases_SE_from_folder(
    const std::string& directory, const std::string& extension);

std::unordered_map<std::string, std::uint64_t> count_total_bases_PE_from_folder(
    const std::string& directory, const std::string& extension, const std::vector<std::string>& suffixes);

#endif // COUNT_BASES_HPP