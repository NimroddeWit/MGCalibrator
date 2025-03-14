#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <vector>

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
    std::string fileExtension = ".fastq"; // Default value for file extension
    std::string metagenomesDir;
    std::string metadataFile;
    std::string countsFile;
    std::string outputFile;

    if (argc == 1) {
        printHelp(argv[0]);
        return 1; // Exit if no arguments are provided
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help") {
            printHelp(argv[0]);
            return 0; // Exit after displaying help
        } else if (arg == "--metagenomes" && i + 1 < argc) {
            metagenomesDir = argv[++i];
        } else if (arg == "--extension" && i + 1 < argc) {
            fileExtension = argv[++i];
        } else if (arg == "--meta" && i + 1 < argc) {
            metadataFile = argv[++i];
        } else if (arg == "--counts" && i + 1 < argc) {
            countsFile = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            outputFile = argv[++i];
        } else {
            std::cerr << "Invalid or unrecognized option: " << arg << std::endl;
            return 1;  // exit if the option is unrecognized
        }
    }

    // Check if all required parameters are provided
    if (metagenomesDir.empty() || metadataFile.empty() || countsFile.empty() || outputFile.empty()) {
        std::cerr << "Error: Missing required parameters.\n";
        printHelp(argv[0]);
        return 1; // Exit if any required parameter is missing
    }

    std::vector<std::pair<std::string, std::unordered_map<char, int>>> sample_nucleotide_counts;

    try {
        for (const auto& entry : std::filesystem::directory_iterator(metagenomesDir)) {
            if (entry.path().extension() == fileExtension) {
                std::unordered_map<char, int> nucleotide_counts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'N', 0}};
                std::ifstream file_in(entry.path());
                std::string line;
                int line_count = 0;

                while (std::getline(file_in, line)) {
                    if (++line_count == 2) {
                        for (char nt : line) {
                            if (nucleotide_counts.find(nt) != nucleotide_counts.end()) {
                                nucleotide_counts[nt]++;
                            }
                        }
                        std::getline(file_in, line); // Skip '+'
                        std::getline(file_in, line); // Skip quality scores
                        line_count = 0; // Reset line counter for the next record
                    }
                }

                // Store the filename and the counts in the vector
                sample_nucleotide_counts.push_back({entry.path().filename(), nucleotide_counts});
            }
        }
        
        // Optionally, output the counts for each file
        for (const auto& [filename, counts] : sample_nucleotide_counts) {
            std::cout << "File: " << filename << "\n";
            for (const auto& [nt, count] : counts) {
                std::cout << nt << ": " << count << '\n';
            }
            std::cout << "-----------------\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error processing files: " << e.what() << '\n';
        return 1;
    }

    return 0;
}