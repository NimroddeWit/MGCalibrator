#ifndef IO_MATRIX_HPP
#define IO_MATRIX_HPP

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

// Function to read a matrix from a CSV file
std::pair<std::vector<std::string>, std::unordered_map<std::string, std::vector<double>>> 
read_matrix(const std::string& filePath);

// Function to write a matrix to a CSV file
void write_matrix(
    const std::string& output_file,
    const std::vector<std::string>& taxa_names,
    const std::unordered_map<std::string, std::vector<double>>& absolute_matrix
);

#endif