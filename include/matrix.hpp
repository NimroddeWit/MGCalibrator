#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <string>
#include <unordered_map>

// ✅ Holds a matrix of numeric values with sample and taxa metadata
struct Matrix {
    std::vector<std::string> taxa_names;   // Column headers
    std::vector<std::string> sample_order; // Row order
    std::unordered_map<std::string, std::vector<double>> data; // sample_id -> abundances
};

#endif 
