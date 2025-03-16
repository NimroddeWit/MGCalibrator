// count_fastq_bases.cpp
#include "count_fastq_bases.h"

std::unordered_map<char, int> count_fastq_bases(const std::string& filePath) {

    std::ifstream file(filePath);

    std::unordered_map<char, int> nucleotide_counts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'N', 0}};
    std::string line;
    int line_count = 0;

    while (std::getline(file, line)) {
        if (++line_count == 2) {
            for (char nt : line) {
                if (nucleotide_counts.find(nt) != nucleotide_counts.end()) {
                    nucleotide_counts[nt]++;
                }
            }
            std::getline(file, line); // Skip '+'
            std::getline(file, line); // Skip quality scores
            line_count = 0; // Reset line counter 
        }
    }
    
    file.close();

    return nucleotide_counts;
}