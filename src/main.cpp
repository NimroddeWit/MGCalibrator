#include <iostream>
#include <vector>
#include <unordered_map>
#include <set> 
#include "count_bases.hpp"
#include "read_metadata.hpp"
#include "io_matrix.hpp"

void printHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "  --metagenomes <path>     Path to metagenomes directory (required)\n"
              << "  --extension <ext>        File extension (default: .fastq)\n"
              << "  --suffixes <s1,s2>       Paired-end suffixes, e.g., _R1,_R2. If omitted, samples are treated as single-end.\n"
              << "  --meta <file>            Metadata file path - Contains 'sample_id' and 'DNA_conc' in ng/μL (required)\n"
              << "  --counts <file>          Counts matrix file path (required)\n"
              << "  --output <file>          Output file for absolute quantification (required)\n"
              << "  --volume <μL>            Processed volume in μL (required)\n"
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
    double processed_volume_uL = -1;

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
            processed_volume_uL = std::stod(argv[++i]);
            if (processed_volume_uL <= 0) {
                std::cerr << "❌ Error: --volume must be greater than 0 μL.\n";
                return 1;
            }
        } else {
            std::cerr << "❌ Invalid option: " << arg << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }

    if (metagenomes_dir.empty() || metadata_file.empty() || counts_file.empty() || output_file.empty() || processed_volume_uL == -1) {
        std::cerr << "❌ Error: --metagenomes, --meta, --counts, --output, and --volume are REQUIRED.\n";
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

      // 🔍 Check sample ID consistency across all data sources
      std::set<std::string> base_ids, meta_ids, count_ids;
      for (const auto& [id, _] : sample_nucleotide_counts) base_ids.insert(id);
      for (const auto& [id, _] : metadata) meta_ids.insert(id);
      for (const auto& [id, _] : counts_matrix) count_ids.insert(id);
  
      if (base_ids != meta_ids || base_ids != count_ids) {
        std::cerr << "❌ Error: Sample IDs do not match across metagenomes, metadata, and counts. "
                  << "(Did you set --suffixes for paired-end metagenomes?)\n";
        return 1;
    }

    // 🔹 Calculate absolute quantification
    std::unordered_map<std::string, std::vector<double>> absolute_quantification;

    for (const auto& [sample_name, counts] : sample_nucleotide_counts) {

        // 1. 📦 DNA mass extracted from sample (in grams)
        // DNA_conc is in ng/µL, volume is in µL -> total DNA in ng -> convert to grams (×1e-9)
        double extracted_mass_g = metadata[sample_name] * processed_volume_uL * 1e-9;

        // 2. Total weight of sequenced DNA in daltons
        double sequenced_mass_daltons =
        (counts[0] * 313.2) +
        (counts[1] * 289.2) +
        (counts[2] * 329.2) +
        (counts[3] * 304.2) + 79.0;

        // 3. Convert to grams
        double sequenced_mass_g = sequenced_mass_daltons * 1.66054e-24;

        // 4. Compute ratio
        double ratio = sequenced_mass_g / extracted_mass_g;

        // 5. Scale each count by the inverse of the ratio
        std::vector<double> absolute_counts;
        for (double abun : counts_matrix[sample_name]) {
            absolute_counts.push_back(abun / ratio);  
        }

        absolute_quantification[sample_name] = absolute_counts;
    }

    // 🔹 Save the results
    write_matrix(output_file, taxa_names, absolute_quantification);

    return 0;
}
