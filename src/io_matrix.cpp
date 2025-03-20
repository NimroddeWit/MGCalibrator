#include "io_matrix.hpp"

std::pair<std::vector<std::string>, std::unordered_map<std::string, std::vector<double>>> 
read_matrix(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("❌ Error: Failed to open matrix file: " + filePath);
    }

    std::vector<std::string> taxa_names;
    std::unordered_map<std::string, std::vector<double>> matrix_data;

    std::string line;
    bool first_row = true;
    size_t expected_columns = 0;  // Track number of columns

    std::cout << "📂 Reading matrix file: " << filePath << "\n";

    while (getline(file, line)) {
        std::istringstream lineStream(line);
        std::string token;

        if (first_row) {
            // ✅ Read header (Extract taxa names, skipping first column)
            getline(lineStream, token, ',');  // Skip "sample_id"
            while (getline(lineStream, token, ',')) {
                taxa_names.push_back(token);
            }
            expected_columns = taxa_names.size();
            first_row = false;
        } else {
            // ✅ Read sample_id and taxa abundances
            std::string sample_id;
            getline(lineStream, sample_id, ',');

            std::vector<double> abundances;
            while (getline(lineStream, token, ',')) {
                try {
                    abundances.push_back(std::stod(token));  // Convert to double
                } catch (const std::exception &e) {
                    throw std::runtime_error("❌ Error: Invalid numeric value in " + filePath + " at sample: " + sample_id);
                }
            }

            // ✅ Validate row consistency
            if (abundances.size() != expected_columns) {
                throw std::runtime_error("❌ Error: Mismatched column count in " + filePath + " at sample: " + sample_id);
            }

            // ✅ Store in matrix
            matrix_data[sample_id] = abundances;
        }
    }

    if (matrix_data.empty()) {
        throw std::runtime_error("❌ Error: No valid matrix data found in file: " + filePath);
    }

    std::cout << "✅ Matrix file successfully loaded: " << filePath << "\n";
    return {taxa_names, matrix_data};
}

void write_matrix(
    const std::string& output_file,
    const std::vector<std::string>& taxa_names,
    const std::unordered_map<std::string, std::vector<double>>& absolute_matrix
) {
    std::ofstream outFile(output_file);
    if (!outFile.is_open()) {
        std::cerr << "❌ Error: Could not open output file: " << output_file << "\n";
        return;
    }

    std::cout << "📝 Writing matrix to: " << output_file << "\n";

    // ✅ Write header row (Sample ID + Taxa Names)
    outFile << "sample_id";
    for (const auto& taxa : taxa_names) {
        outFile << "," << taxa;
    }
    outFile << "\n";

    // ✅ Write data rows
    for (const auto& [sample_name, abundances] : absolute_matrix) {
        outFile << sample_name;
        for (double value : abundances) {
            outFile << "," << value;
        }
        outFile << "\n";
    }

    std::cout << "✅ Matrix successfully saved: " << output_file << "\n";
    outFile.close();
}