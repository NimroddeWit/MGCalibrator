#ifndef IO_MATRIX_HPP
#define IO_MATRIX_HPP

#include "matrix.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_set>
#include <stdexcept>

// ✅ Read matrix from file into Matrix struct
Matrix read_matrix(const std::string& filePath);

// ✅ Write Matrix struct to CSV file
void write_matrix(const std::string& output_file, const Matrix& matrix);

#endif 