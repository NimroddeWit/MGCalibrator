#include "io_matrix.hpp"

// ✅ Read matrix file into a Matrix struct (preserves sample order)
Matrix read_matrix(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("❌ Error: Failed to open matrix file: " + filePath);
    }

    Matrix matrix;
    std::unordered_set<std::string> seen_taxa;
    std::unordered_set<std::string> seen_samples;

    std::string line;
    bool first_row = true;
    size_t expected_columns = 0;

    std::cout << "📂 Reading matrix file: " << filePath << "\n";

    while (getline(file, line)) {
        std::istringstream lineStream(line);
        std::string token;

        if (first_row) {
            getline(lineStream, token, ',');  // skip "sample_id"
            while (getline(lineStream, token, ',')) {
                if (!seen_taxa.insert(token).second) {
                    throw std::runtime_error("❌ Error: Duplicate taxa name in header: " + token);
                }
                matrix.taxa_names.push_back(token);
            }
            expected_columns = matrix.taxa_names.size();
            first_row = false;
        } else {
            std::string sample_id;
            getline(lineStream, sample_id, ',');

            if (!seen_samples.insert(sample_id).second) {
                throw std::runtime_error("❌ Error: Duplicate sample_id: " + sample_id);
            }

            matrix.sample_order.push_back(sample_id);
            std::vector<double> abundances;
            while (getline(lineStream, token, ',')) {
                try {
                    double value = std::stod(token);
                    if (value < 0) {
                        throw std::runtime_error("❌ Error: Negative abundance value in " + filePath + " at sample: " + sample_id);
                    }
                    abundances.push_back(value);
                } catch (...) {
                    throw std::runtime_error("❌ Error: Invalid numeric in " + filePath + " at sample: " + sample_id);
                }
            }

            if (abundances.size() != expected_columns) {
                throw std::runtime_error("❌ Error: Column count mismatch at sample: " + sample_id);
            }

            matrix.data[sample_id] = abundances;
        }
    }

    if (matrix.data.empty()) {
        throw std::runtime_error("❌ Error: No matrix data found in file: " + filePath);
    }

    std::cout << "✅ Matrix file successfully loaded: " << filePath << "\n";
    return matrix;
}

// ✅ Write matrix from Matrix struct (preserves sample order)
void write_matrix(const std::string& output_file, const Matrix& matrix) {
    std::ofstream outFile(output_file);
    if (!outFile.is_open()) {
        std::cerr << "❌ Error: Could not open output file: " << output_file << "\n";
        return;
    }

    std::cout << "📝 Writing matrix to: " << output_file << "\n";

    outFile << "sample_id";
    for (const auto& taxa : matrix.taxa_names) {
        outFile << "," << taxa;
    }
    outFile << "\n";

    for (const auto& sample : matrix.sample_order) {
        outFile << sample;
        for (double value : matrix.data.at(sample)) {
            outFile << "," << value;
        }
        outFile << "\n";
    }

    std::cout << "✅ Matrix successfully saved: " << output_file << "\n";
    outFile.close();
}