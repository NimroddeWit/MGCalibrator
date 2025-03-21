#include "count_bases.hpp"

// 🔹 **Count bases in a SINGLE-END FASTQ/FASTA file**
std::array<size_t, seqan3::dna5::alphabet_size> count_bases_SE(const std::string& filePath) {
    std::array<size_t, seqan3::dna5::alphabet_size> nucleotide_counts{}; 

    try {
        seqan3::sequence_file_input file{filePath};
        for (auto &record : file) {
            for (seqan3::dna5 symbol : record.sequence()) {
                ++nucleotide_counts[symbol.to_rank()];
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error reading file " << filePath << ": " << e.what() << std::endl;
    }

    return nucleotide_counts;
}

// 🔹 **Count bases in a PAIRED-END FASTQ file**
std::array<size_t, seqan3::dna5::alphabet_size> count_bases_PE(const std::string& file1, const std::string& file2) {
    std::array<size_t, seqan3::dna5::alphabet_size> combined_counts{}; 

    try {
        seqan3::sequence_file_input fin1{file1};
        seqan3::sequence_file_input fin2{file2};

        for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) {
            if (rec1.id() != rec2.id()) {
                throw std::runtime_error("Mismatched paired-end reads in: " + file1 + " and " + file2);
            }

            for (seqan3::dna5 symbol : rec1.sequence()) {
                ++combined_counts[symbol.to_rank()];
            }
            for (seqan3::dna5 symbol : rec2.sequence()) {
                ++combined_counts[symbol.to_rank()];
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error processing paired-end files: " << e.what() << std::endl;
    }

    return combined_counts;
}

// 🔹 **Process ALL SINGLE-END FASTQ/FASTA FILES in a folder**
std::vector<std::pair<std::string, std::array<size_t, 5>>> count_bases_SE_from_folder(const std::string& directory, const std::string& extension) {
    
    std::vector<std::pair<std::string, std::array<size_t, 5>>> sample_nucleotide_counts;
    std::unordered_set<std::string> sample_names;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        std::string ext = entry.path().extension().string();
        if (ext == extension || ext == (extension + ".gz")) {
            std::string sample_name = entry.path().stem().string();
            if (!sample_names.insert(sample_name).second) {
                throw std::runtime_error("❌ Duplicate sample name detected: " + sample_name);
            }

            std::array<size_t, 5> counts = count_bases_SE(entry.path().string());
            sample_nucleotide_counts.emplace_back(sample_name, counts);
        }
    }

    return sample_nucleotide_counts;
}

// 🔹 **Process ALL PAIRED-END FASTQ FILES in a folder**
std::vector<std::pair<std::string, std::array<size_t, 5>>> count_bases_PE_from_folder(
    const std::string& directory,
    const std::string& extension,
    const std::vector<std::string>& suffixes)
{
    std::vector<std::pair<std::string, std::array<size_t, 5>>> sample_nucleotide_counts;
    std::unordered_map<std::string, std::string> r1_files;
    std::unordered_map<std::string, std::string> r2_files;

    // 📂 Scan directory to find `_R1` and `_R2` files
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        std::string filename = entry.path().stem().string();
        std::string ext = entry.path().extension().string();

        if (ext == extension || ext == (extension + ".gz")) {
            if (filename.size() >= suffixes[0].size() &&
                filename.compare(filename.size() - suffixes[0].size(), suffixes[0].size(), suffixes[0]) == 0)
            {
                std::string base_name = filename.substr(0, filename.size() - suffixes[0].size());
                r1_files[base_name] = entry.path().string();
            }
            else if (filename.size() >= suffixes[1].size() &&
                     filename.compare(filename.size() - suffixes[1].size(), suffixes[1].size(), suffixes[1]) == 0)
            {
                std::string base_name = filename.substr(0, filename.size() - suffixes[1].size());
                r2_files[base_name] = entry.path().string();
            }
        }
    }

    // ✅ Process only paired R1 & R2 files
    for (const auto& [base_name, r1_file] : r1_files) {
        if (r2_files.contains(base_name)) {
            std::array<size_t, 5> counts = count_bases_PE(r1_file, r2_files[base_name]);
            sample_nucleotide_counts.emplace_back(base_name, counts);
        } else {
            std::cerr << "⚠️ Warning: Missing paired file for base name: " << base_name << "\n";
        }
    }

    // ✅ Check for duplicate sample names
    std::unordered_set<std::string> sample_names;
    for (const auto& [sample, _] : sample_nucleotide_counts) {
        if (!sample_names.insert(sample).second) {
            throw std::runtime_error("❌ Error: Duplicate sample name found: " + sample);
        }
    }

    return sample_nucleotide_counts;
}
