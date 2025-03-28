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
              << "  --debug                  Enable debug mode and save intermediate files\n"
              << "  --help                   Display this help and exit\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printHelp(argv[0]);
        return 1;
    }

    std::string file_extension = ".fastq";
    std::string metagenomes_dir, metadata_file, counts_file, output_file;
    std::vector<std::string> suffixes;
    double processed_volume_uL = -1;
    bool debug = false;
    std::filesystem::path debug_dir = "debug";

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
                std::cerr << "\u274c Error: --volume must be greater than 0 \u00b5L.\n";
                return 1;
            }
        } else if (arg == "--debug") {
            debug = true;
            if (!std::filesystem::exists(debug_dir)) {
                std::filesystem::create_directory(debug_dir);
            }
        } else {
            std::cerr << "\u274c Invalid option: " << arg << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }

    if (metagenomes_dir.empty() || metadata_file.empty() || counts_file.empty() || output_file.empty() || processed_volume_uL == -1) {
        std::cerr << "\u274c Error: --metagenomes, --meta, --counts, --output, and --volume are REQUIRED.\n";
        printHelp(argv[0]);
        return 1;
    }

    Matrix abundance_matrix;
    try {
        abundance_matrix = read_matrix(counts_file);
    } catch (const std::exception& e) {
        std::cerr << "\u274c Error loading counts matrix: " << e.what() << "\n";
        return 1;
    }

    if (debug) {
        write_matrix((debug_dir / "abundance_matrix_debug.csv").string(), abundance_matrix);
    }

    std::unordered_map<std::string, double> metadata;
    try {
        metadata = read_metadata(metadata_file);
    } catch (const std::exception& e) {
        std::cerr << "\u274c Error loading metadata: " << e.what() << "\n";
        return 1;
    }

    if (debug) {
        std::ofstream meta_debug_file(debug_dir / "metadata_debug.csv");
        meta_debug_file << "sample_id,DNA_conc\n";
        for (const auto& [id, conc] : metadata) {
            meta_debug_file << id << "," << conc << "\n";
        }
    }

    std::unordered_map<std::string, std::uint64_t> sequenced_bp_map = suffixes.empty()
        ? count_total_bases_SE_from_folder(metagenomes_dir, file_extension)
        : count_total_bases_PE_from_folder(metagenomes_dir, file_extension, suffixes);

    if (debug) {
        std::ofstream base_debug_file(debug_dir / "base_counts_debug.csv");
        base_debug_file << "sample_id,total_bases\n";
        for (const auto& [id, count] : sequenced_bp_map) {
            base_debug_file << id << "," << count << "\n";
        }
    }

    const double weightPerBP = 660.0 / 6.022e23; // g/bp

    std::unordered_map<std::string, double> extracted_dna_mass;
    std::unordered_map<std::string, double> total_bp_dna;
    std::unordered_map<std::string, double> fraction_sequenced;

    Matrix absolute_matrix;
    absolute_matrix.taxa_names = abundance_matrix.taxa_names;
    absolute_matrix.sample_order = abundance_matrix.sample_order;

    for (const auto& sample_id : abundance_matrix.sample_order) {
        double extracted_g = metadata.at(sample_id) * processed_volume_uL * 1e-9;
        std::uint64_t theoretical_bp = static_cast<std::uint64_t>(extracted_g / weightPerBP);
        double sequenced_bp = static_cast<double>(sequenced_bp_map.at(sample_id));

        extracted_dna_mass[sample_id] = extracted_g;
        total_bp_dna[sample_id] = theoretical_bp;
        fraction_sequenced[sample_id] = sequenced_bp / theoretical_bp;

        std::vector<double> absolute_counts;
        for (double aligned_bp : abundance_matrix.data.at(sample_id)) {
            absolute_counts.push_back(aligned_bp / fraction_sequenced[sample_id]);
        }
        absolute_matrix.data[sample_id] = absolute_counts;
    }

    if (debug) {
        std::ofstream frac_debug_file(debug_dir / "fraction_sequenced_debug.csv");
        frac_debug_file << "sample_id,DNA_mass_g,total_bp,sequenced_bp,fraction_sequenced\n";
        for (const auto& sample_id : abundance_matrix.sample_order) {
            frac_debug_file << sample_id << ","
                            << extracted_dna_mass[sample_id] << ","
                            << total_bp_dna[sample_id] << ","
                            << sequenced_bp_map.at(sample_id) << ","
                            << fraction_sequenced[sample_id] << "\n";
        }
    }

    write_matrix(output_file, absolute_matrix);

    return 0;
}