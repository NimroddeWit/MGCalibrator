#ifndef COUNT_BASES_HPP
#define COUNT_BASES_HPP

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/utility/views/zip.hpp>

// **Function Prototypes**
std::array<size_t, seqan3::dna5::alphabet_size> count_bases_SE(const std::string& filePath);
std::array<size_t, seqan3::dna5::alphabet_size> count_bases_PE(const std::string& file1, const std::string& file2);
std::vector<std::pair<std::string, std::array<size_t, 5>>> count_bases_SE_from_folder(
    const std::string& directory, const std::string& extension);
std::vector<std::pair<std::string, std::array<size_t, 5>>> count_bases_PE_from_folder(
    const std::string& directory, const std::string& extension, const std::vector<std::string>& suffixes);

#endif


