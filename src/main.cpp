#include <iostream>
#include <vector>
#include <unordered_map>
#include "count_bases.hpp"
#include "read_metadata.hpp"
#include "io_matrix.hpp"

void printHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "  --metagenomes <path>     Path to metagenomes directory\n"
              << "  --extension <ext>        File extension (default: .fastq)\n"
              << "  --suffixes <R1,R2>       Paired-end suffixes (e.g., _R1,_R2)\n"
              << "  --meta <file>            Metadata file path\n"
              << "  --counts <file>          Counts matrix file path\n"
              << "  --output <file>          Output file for absolute quantification\n"
              << "  --volume <mL>            Processed volume in uL (default: 500 umL)\n"
              << "  --help                   Display this help and exit\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printHelp(argv[0]);
        return 1;
    }

    // 🔹 Default values
    std::string file_extension = ".fastq"; 
    std::string metagenomes_dir, metadata_file, counts_file, output_file;
    std::vector<std::string> suffixes;
    double PROCESSED_VOLUME_mL = 500;  // Default: 500 mL 

    // 🔹 Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help") {
            printHelp(argv[0]);
            return 0;
        } else if (arg == "--metagenomes" && i + 1 < argc) {
            metagenomes_dir = argv[++i];
        } else if (arg == "--extension" && i + 1 < argc) {
            file_extension = argv[++i];
        } else if (arg == "--suffixes" && i + 1 < argc) {
            std::istringstream ss(argv[++i]);
            std::string token;
            while (std::getline(ss, token, ',')) {
                suffixes.push_back(token);
            }
        } else if (arg == "--meta" && i + 1 < argc) {
            metadata_file = argv[++i];
        } else if (arg == "--counts" && i + 1 < argc) {
            counts_file = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--volume" && i + 1 < argc) {
            PROCESSED_VOLUME_mL = std::stod(argv[++i]);
        } else {
            std::cerr << "❌ Invalid option: " << arg << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }

    // 🔹 Ensure required parameters are provided
    if (metagenomes_dir.empty() || metadata_file.empty() || counts_file.empty() || output_file.empty()) {
        std::cerr << "❌ Error: --metagenomes, --meta, --counts, and --output are REQUIRED.\n";
        printHelp(argv[0]);
        return 1;
    }

    // 🔹 Load metadata
    std::unordered_map<std::string, double> metadata;
    try {
        metadata = read_metadata(metadata_file);
    } catch (const std::exception &e) {
        std::cerr << "❌ Error: Failed to load metadata -> " << e.what() << "\n";
        return 1;
    }

    // 🔹 Load counts matrix
    std::vector<std::string> taxa_names;
    std::unordered_map<std::string, std::vector<double>> counts_matrix;
    try {
        std::tie(taxa_names, counts_matrix) = read_matrix(counts_file);
    } catch (const std::exception &e) {
        std::cerr << "❌ Error: Failed to load counts matrix -> " << e.what() << "\n";
        return 1;
    }

    // 🔹 Perform base counting (Single-End or Paired-End)
    auto sample_nucleotide_counts = suffixes.empty()
        ? count_bases_SE_from_folder(metagenomes_dir, file_extension)
        : count_bases_PE_from_folder(metagenomes_dir, file_extension, suffixes);

    // 🔹 Calculate absolute quantification
    std::unordered_map<std::string, std::vector<double>> absolute_quantification;
    for (const auto& [sample_name, counts] : sample_nucleotide_counts) {
        if (metadata.find(sample_name) == metadata.end()) {
            std::cerr << "❌ Error: Missing metadata for sample " << sample_name << "\n";
            return 1;
        }

        // Compute total DNA mass (grams) for the sample
        const double PROCESSED_VOLUME_UL = PROCESSED_VOLUME_mL * 1000;
        double sampleDNA_g = metadata[sample_name] * PROCESSED_VOLUME_UL * 1e-9;

        // Compute total base pairs
        double total_weight = (counts[0] * 313.2) + (counts[1] * 289.2) +
                              (counts[2] * 329.2) + (counts[3] * 304.2) + 79.0;
        double sampleDNA_bp = sampleDNA_g / (total_weight / (counts[0] + counts[1] + counts[2] + counts[3]));

        // Compute fraction sequenced
        double seqDNA_bp = counts_matrix[sample_name][0];  // Total sequenced DNA
        double fraction_sequenced = seqDNA_bp / sampleDNA_bp;

        // Compute absolute quantification for each taxon
        std::vector<double> absolute_counts;
        for (double basepair : counts_matrix[sample_name]) {
            absolute_counts.push_back(basepair / fraction_sequenced);
        }
        absolute_quantification[sample_name] = absolute_counts;
    }

    // 🔹 Save the results
    write_matrix(output_file, taxa_names, absolute_quantification);

    return 0;
}