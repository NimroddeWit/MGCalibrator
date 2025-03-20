#ifndef READ_METADATA_HPP
#define READ_METADATA_HPP

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

std::unordered_map<std::string, double> read_metadata(const std::string& filePath);

#endif 