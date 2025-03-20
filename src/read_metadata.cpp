#include "read_metadata.hpp"

std::unordered_map<std::string, double> read_metadata(const std::string& filePath) {
    std::unordered_map<std::string, double> metadata;
    std::ifstream file(filePath);

    if (!file) {
        throw std::runtime_error("❌ Error: Unable to open metadata file: " + filePath);
    }

    std::string header, col1, col2;
    if (!std::getline(file, header)) {
        throw std::runtime_error("❌ Error: Metadata file is empty: " + filePath);
    }

    std::istringstream headerStream(header);
    if (!(std::getline(headerStream, col1, ',') && std::getline(headerStream, col2, ',')) ||
        col1 != "sample_id" || col2 != "DNA_prop") {
        throw std::runtime_error("❌ Error: Invalid header in file: " + filePath + ". Expected: 'sample_id,DNA_prop'");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string sample_id;
        double concentration;

        if (!(std::getline(lineStream, sample_id, ',') && (lineStream >> concentration))) {
            std::cerr << "⚠️ Warning: Skipping malformed line in " << filePath << " -> " << line << std::endl;
            continue;
        }

        if (concentration < 0.0 || concentration > 1.0) {
            std::cerr << "⚠️ Warning: DNA concentration out of range [0,1] in " << filePath << " -> "
                      << sample_id << ": " << concentration << std::endl;
            continue;
        }

        if (metadata.find(sample_id) != metadata.end()) {
            std::cerr << "⚠️ Warning: Duplicate entry found for sample_id: " << sample_id << " in " << filePath
                      << ". Keeping the first occurrence." << std::endl;
            continue;
        }

        metadata[sample_id] = concentration;
    }

    if (metadata.empty()) {
        throw std::runtime_error("❌ Error: No valid metadata entries found in file: " + filePath);
    }

    return metadata;
}