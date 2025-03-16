#ifndef READ_METADATA_H
#define READ_METADATA_H

#include <string>
#include <unordered_map>

std::unordered_map<std::string, double> read_metadata(const std::string& filePath);

#endif // READ_METADATA_H