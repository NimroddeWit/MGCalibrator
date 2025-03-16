#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <vector>
#include "count_fastq_bases.h"
#include "read_metadata.h"

void printHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "Options:\n"
              << "  --metagenomes <path>  Specify the path to the metagenomes directory\n"
              << "  --extension <ext>     Specify the file extension (default: .fastq)\n"
              << "  --meta <file>         Specify the metadata file\n"
              << "  --counts <file>       Specify the counts file\n"
              << "  --output <file>       Specify the output matrix file\n"
              << "  --help                Display this help and exit\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printHelp(argv[0]);
        return 1; // Exit if no arguments are provided
    }

    std::string file_extension = ".fastq"; 
    std::string metagenomes_dir;
    std::string metadata_file;
    std::string counts_file;
    std::string output_file;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help") {
            printHelp(argv[0]);
            return 0; 
        } else if (arg == "--metagenomes" && i + 1 < argc) {
            metagenomes_dir = argv[++i];
        } else if (arg == "--extension" && i + 1 < argc) {
            file_extension = argv[++i];
        } else if (arg == "--meta" && i + 1 < argc) {
            metadata_file = argv[++i];
        } else if (arg == "--counts" && i + 1 < argc) {
            counts_file = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else {
            std::cerr << "Invalid or unrecognized option: " << arg << std::endl;
            return 1;
        }
    }

    if (metagenomes_dir.empty() || metadata_file.empty() || counts_file.empty() || output_file.empty()) {
        std::cerr << "Error: Missing required parameters.\n";
        printHelp(argv[0]);
        return 1;
    }

    std::unordered_map<std::string, double> metadata;
    try {
        metadata = read_metadata(metadata_file);
    } catch (const std::exception& e) {
        std::cerr << "Error reading metadata: " << e.what() << std::endl;
        return 2; // Exit the program with an error code
    }

    std::vector<std::pair<std::string, std::unordered_map<char, int>>> sample_nucleotide_counts;
    try {
        for (const auto& entry : std::filesystem::directory_iterator(metagenomes_dir)) {
            if (entry.path().extension() == file_extension) {
                auto nucleotide_counts = count_fastq_bases(entry.path());
                sample_nucleotide_counts.push_back({entry.path().filename().stem(), nucleotide_counts});
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error processing files: " << e.what() << '\n';
        return 1;
    }

    // Output the nucleotide counts and metadata for each sample
    for (const auto& sample : sample_nucleotide_counts) {
        std::cout << "Sample: " << sample.first << std::endl;
        for (const auto& pair : sample.second) {
            std::cout << "  " << pair.first << ": " << pair.second << std::endl;
        }
        if (metadata.find(sample.first) != metadata.end()) {
            std::cout << "DNA Concentration: " << metadata[sample.first] << std::endl;
        } else {
            std::cout << "DNA Concentration: Not available" << std::endl;
        }
        std::cout << "-------------------------------------\n";
    }

    return 0;
}