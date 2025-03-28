#include "read_metadata.hpp"

std::unordered_map<std::string, double> read_metadata(const std::string& filePath) {
    std::unordered_map<std::string, double> metadata;
    std::ifstream file(filePath);

    std::cout << "📂 Reading DNA concentration file: " << filePath << "\n";

    if (!file.is_open()) {
        throw std::runtime_error("❌ Error: Unable to open metadata file: " + filePath);
    }

    std::string header, col1, col2;
    if (!std::getline(file, header)) {
        throw std::runtime_error("❌ Error: Metadata file is empty: " + filePath);
    }

    std::istringstream headerStream(header);
    if (!(std::getline(headerStream, col1, ',') && std::getline(headerStream, col2, ',')) ||
        col1 != "sample_id" || col2 != "DNA_conc") {
        throw std::runtime_error("❌ Error: Invalid header. Expected 'sample_id,DNA_conc'");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string sample_id;
        double concentration;

        if (!(std::getline(lineStream, sample_id, ',') && (lineStream >> concentration))) {
            std::cerr << "⚠️ Warning: Skipping malformed line: " << line << std::endl;
            continue;
        }

        if (concentration < 0.0) {
            std::cerr << "⚠️ Warning: DNA concentration cannot be negative for sample " << sample_id << std::endl;
            continue;
        }

        if (metadata.find(sample_id) != metadata.end()) {
            throw std::runtime_error("❌ Error: Duplicate sample_id '" + sample_id + "' in metadata file.");
        }

        metadata[sample_id] = concentration;
    }

    if (metadata.empty()) {
        throw std::runtime_error("❌ Error: No valid metadata entries found.");
    }

    std::cout << "✅ DNA concentration file successfully loaded: " << filePath << "\n";

    return metadata;
}