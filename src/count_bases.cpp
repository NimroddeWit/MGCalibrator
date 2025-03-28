#include "count_bases.hpp"

// 🔹 Count total bases in a SINGLE-END FASTQ/FASTA file
// 📄 Each base is counted from the input file regardless of nucleotide type
std::uint64_t count_total_bases_SE(const std::string& filePath) {
    std::uint64_t total = 0;
    try {
        seqan3::sequence_file_input file{filePath};
        for (auto &record : file) {
            total += record.sequence().size();
        }
    } catch (const std::exception &e) {
        std::cerr << "❌ Error reading file " << filePath << ": " << e.what() << std::endl;
    }
    return total;
}

// 🔹 Count total bases in a PAIRED-END FASTQ file
// 📄 Counts total number of bases from both R1 and R2 reads
std::uint64_t count_total_bases_PE(const std::string& file1, const std::string& file2) {
    std::uint64_t total = 0;
    try {
        seqan3::sequence_file_input fin1{file1};
        seqan3::sequence_file_input fin2{file2};

        for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) {
            if (rec1.id() != rec2.id()) {
                throw std::runtime_error("Mismatched paired-end reads in: " + file1 + " and " + file2);
            }
            total += rec1.sequence().size();
            total += rec2.sequence().size();
        }
    } catch (const std::exception &e) {
        std::cerr << "❌ Error processing paired-end files: " << e.what() << std::endl;
    }
    return total;
}

// 🔹 Process ALL SINGLE-END FASTQ/FASTA FILES in a folder
// 📁 Iterates over all single-end files and sums total base counts
std::unordered_map<std::string, std::uint64_t> count_total_bases_SE_from_folder(const std::string& directory, const std::string& extension) {
    std::unordered_map<std::string, std::uint64_t> sample_counts;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        std::string ext = entry.path().extension().string();
        if (ext == extension || ext == (extension + ".gz")) {
            std::string sample_name = entry.path().stem().string();
            std::uint64_t count = count_total_bases_SE(entry.path().string());
            sample_counts[sample_name] = count;
        }
    }
    return sample_counts;
}

// 🔹 Process ALL PAIRED-END FASTQ FILES in a folder
// 📁 Matches R1 and R2 files based on suffixes and processes them
std::unordered_map<std::string, std::uint64_t> count_total_bases_PE_from_folder(
    const std::string& directory,
    const std::string& extension,
    const std::vector<std::string>& suffixes) {

    std::unordered_map<std::string, std::uint64_t> sample_counts;
    std::unordered_map<std::string, std::string> r1_files;
    std::unordered_map<std::string, std::string> r2_files;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        std::string filename = entry.path().stem().string();
        std::string ext = entry.path().extension().string();

        if (ext == extension || ext == (extension + ".gz")) {
            if (filename.ends_with(suffixes[0])) {
                std::string base_name = filename.substr(0, filename.size() - suffixes[0].size());
                r1_files[base_name] = entry.path().string();
            } else if (filename.ends_with(suffixes[1])) {
                std::string base_name = filename.substr(0, filename.size() - suffixes[1].size());
                r2_files[base_name] = entry.path().string();
            }
        }
    }

    for (const auto& [base_name, r1_path] : r1_files) {
        if (r2_files.contains(base_name)) {
            std::uint64_t count = count_total_bases_PE(r1_path, r2_files[base_name]);
            sample_counts[base_name] = count;
        }
    }

    return sample_counts;
}