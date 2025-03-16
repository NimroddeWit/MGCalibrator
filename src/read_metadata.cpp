#include "read_metadata.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

std::unordered_map<std::string, double> read_metadata(const std::string& filePath) {
    std::unordered_map<std::string, double> metadata;
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the metadata file: " + filePath);
    }

    std::string header;
    getline(file, header);
    std::istringstream headerStream(header);
    std::string header1, header2;

    getline(headerStream, header1, ',');
    getline(headerStream, header2, ',');

    if (header1 != "sample_id" || header2 != "DNA_prop") {
        throw std::runtime_error("Invalid header format in metadata file: Expected 'sample_id,DNA_conc'");
    }

    std::string line;
    while (getline(file, line)) {
        std::istringstream lineStream(line);
        std::string sample_id;
        double concentration;
        if (!(getline(lineStream, sample_id, ',') && (lineStream >> concentration))) {
            throw std::runtime_error("Failed to parse line: " + line);
        }

        // Check if the concentration value is within the expected range [0, 1]
        if (concentration < 0.0 || concentration > 1.0) {
            std::cerr << "Invalid DNA concentration value: " << concentration << " in sample " << sample_id << std::endl;
            throw std::runtime_error("DNA concentration value out of bounds [0,1]");
        }

        metadata[sample_id] = concentration;
    }

    return metadata;
}