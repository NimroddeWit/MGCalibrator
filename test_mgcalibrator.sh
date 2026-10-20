# Test BAM mode
mgcalibrator \
    --bam_folder data/input/ \
    --dna_mass data/input/dna_mass.csv \
    --cluster data/input/reference_clusters.csv \
    --bin data/input/reference_bins.csv \
    --output data/output/output_bam.csv \
    --scaling_factors data/output/scaling_factors.txt

# Test SingleM + FASTQ mode
mgcalibrator \
    --singlem_folder data/input/singlem/ \
    --fastq_folder data/input/fastq/ \
    --dna_mass data/input/dna_mass.csv \
    --output data/output/output_singlem.csv \
    --scaling_factors data/output/scaling_factors.txt

# Test Sylph + FASTQ mode (two databases, with confidence intervals)
mgcalibrator \
    --sylph_folder data/input/sylph/gtdb data/input/sylph/virus \
    --sylph_profile_folder data/input/sylph_profiles/ \
    --fastq_folder data/input/fastq/ \
    --dna_mass data/input/dna_mass.csv \
    --output data/output/output_sylph.csv \
    --scaling_factors data/output/scaling_factors.txt
