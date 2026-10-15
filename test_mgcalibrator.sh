# Test BAM mode
mgcalibrator \
    --bam_folder data/input/ \
    --dna_mass data/input/dna_mass.csv \
    --cluster data/input/reference_clusters.csv \
    --bin data/input/reference_bins.csv \
    --output data/output/output.csv \
    --scaling_factors data/output/scaling_factors.txt 

# Test SingleM + FASTQ mode
mgcalibrator \
    --fastq_folder data/input/fastq/ \
    --singlem_folder data/input/singlem/ \
    --dna_mass data/input/dna_mass.csv \
    --output data/output/output.csv \
    --scaling_factors data/output/scaling_factors.txt