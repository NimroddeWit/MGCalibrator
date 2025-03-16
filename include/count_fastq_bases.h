#ifndef COUNT_FASTQ_BASES_H
#define COUNT_FASTQ_BASES_H

#include <unordered_map>
#include <fstream>

std::unordered_map<char, int> count_fastq_bases(const std::string& filePath);

#endif // COUNT_FASTQ_BASES_H