import pandas as pd
import numpy as np
from .parser import calculate_total_base_pairs

def compute_absolute_abundance(counts_df, dna_conc, volume, fastq_files):
    """
    Compute absolute abundance by scaling raw counts with a scaling factor.
    
    Formula: absolute_abundance = raw_counts * scaling_factor
    Where: scaling_factor = initial_dna_weight / final_dna_weight
           initial_dna_weight = dna_conc * volume
           final_dna_weight = calculated from base pairs in FASTQ files
    """
    if not fastq_files:
        raise ValueError("FASTQ files are required to calculate DNA weight from base pairs")
    
    # Calculate initial DNA weight for each sample
    initial_dna_weight = {sample: dna_conc[sample] * volume[sample] for sample in counts_df.columns}
    
    # Calculate final DNA weight from FASTQ files
    sample_base_pairs = calculate_total_base_pairs(fastq_files)
    # Convert base pairs to DNA weight (assuming average molecular weight per base pair)
    # Using approximately 650 Da per base pair (average of A, T, G, C)
    # 1 Da = 1.66054e-15 ng, so 650 Da = 1.079e-12 ng per base pair
    final_dna_weight = {sample: sample_base_pairs.get(sample, 0) * 1.079e-12  # Convert to ng
                       for sample in counts_df.columns}
    
    # Check for samples with no FASTQ data
    for sample in counts_df.columns:
        if final_dna_weight[sample] == 0:
            raise ValueError(f"No base pairs found for sample {sample} in FASTQ files")
    
    # Calculate scaling factors
    scaling_factors = {sample: initial_dna_weight[sample] / final_dna_weight[sample] 
                      for sample in counts_df.columns}
    
    print("Scaling factors:", scaling_factors)
    
    absolute = counts_df.multiply(pd.Series(scaling_factors))
    return absolute, scaling_factors

def compute_absolute_abundance_with_error(counts_df, dna_conc, volume, fastq_files,
                                        n_monte_carlo=1000,
                                        alpha=0.5):
    """
    Compute absolute abundance with 95% confidence intervals using a Bayesian approach.

    The method adds a Dirichlet prior to the counts and then uses Monte Carlo sampling
    from the posterior distribution to estimate uncertainty.

    SCALING FACTOR ERROR PROPAGATION:
    The scaling factor directly affects the width of confidence intervals. Samples with
    higher scaling factors (i.e., higher initial DNA concentration relative to what was
    sequenced) will have proportionally larger confidence intervals, reflecting
    the increased uncertainty when extrapolating from a smaller sequenced fraction.

    Framework:
    1. Model counts with a Multinomial distribution.
    2. Use a Dirichlet prior on the proportions, which is equivalent to adding a pseudocount (alpha).
    3. The posterior is a Dirichlet distribution. We draw Monte Carlo samples of compositions from it.
    4. These are converted to counts and scaled to get absolute abundance samples.
    5. Confidence intervals are derived from the distribution of these absolute abundance samples.

    Parameters:
    -----------
    counts_df : pd.DataFrame
        DataFrame with features (genes/taxa) as rows and samples as columns
    dna_conc : dict
        Dictionary mapping sample names to DNA concentrations
    volume : dict
        Dictionary mapping sample names to volumes
    fastq_files : list
        List of FASTQ files to calculate final DNA weight (mandatory)
    n_monte_carlo : int
        Number of Monte Carlo samples for confidence interval estimation
    alpha : float
        Dirichlet prior parameter (pseudocount to add to each feature).

    Returns:
    --------
    tuple: (absolute_abundance, lower_ci, upper_ci, zero_replaced_counts, scaling_factors)
        - absolute_abundance: DataFrame with point estimates (posterior mean)
        - lower_ci: DataFrame with lower 95% confidence bounds
        - upper_ci: DataFrame with upper 95% confidence bounds
        - zero_replaced_counts: DataFrame with posterior mean counts for transparency
        - scaling_factors: Dictionary of scaling factors for each sample
    """
    if not fastq_files:
        raise ValueError("FASTQ files are required to calculate DNA weight from base pairs")

    # Calculate initial DNA weight for each sample
    initial_dna_weight = {sample: dna_conc[sample] * volume[sample] for sample in counts_df.columns}

    # Calculate final DNA weight from FASTQ files
    sample_base_pairs = calculate_total_base_pairs(fastq_files)
    final_dna_weight = {sample: sample_base_pairs.get(sample, 0) * 1.079e-12  # Convert to ng
                       for sample in counts_df.columns}

    # Check for samples with no FASTQ data
    for sample in counts_df.columns:
        if final_dna_weight[sample] == 0:
            raise ValueError(f"No base pairs found for sample {sample} in FASTQ files")

    # Calculate scaling factors
    scaling_factors = {sample: initial_dna_weight[sample] / final_dna_weight[sample]
                      for sample in counts_df.columns}

    print("Scaling factors:", scaling_factors)
    print(f"Note: Samples with higher scaling factors will have proportionally larger confidence intervals")
    print(f"Using Dirichlet prior for Monte Carlo sampling with alpha={alpha}")

    counts_matrix = counts_df.values.astype(float)
    prior = alpha

    # Initialize arrays to store Monte Carlo results
    mc_results = np.zeros((n_monte_carlo, counts_df.shape[0], counts_df.shape[1]))

    for sample_idx, sample in enumerate(counts_df.columns):
        sample_counts = counts_matrix[:, sample_idx]
        scaling_factor = scaling_factors[sample]
        total_counts = sample_counts.sum()

        if total_counts == 0:
            mc_results[:, :, sample_idx] = 0
            continue

        # Posterior is Dirichlet(counts + prior)
        dirichlet_params = sample_counts + prior

        #Monte Carlo samples from the posterior distribution
        #sample from a dirichlet distribution `n_monte_carlo` times
        sampled_proportions = np.random.dirichlet(dirichlet_params, size=n_monte_carlo)

        # For each set of proportions, sample from a multinomial to get counts
        #adds the sampling noise from the sequencing process
        sampled_counts = np.array([np.random.multinomial(int(total_counts), p) for p in sampled_proportions])

        #Apply scaling factor to get absolute abundance
        sampled_absolute = sampled_counts * scaling_factor
        mc_results[:, :, sample_idx] = sampled_absolute.T

    # Compute point estimate and CIs from the MC samples.
    absolute_point_values = np.mean(mc_results, axis=0)
    lower_ci_values = np.percentile(mc_results, 2.5, axis=0)
    upper_ci_values = np.percentile(mc_results, 97.5, axis=0)

    # Convert back to DataFrames
    absolute_point = pd.DataFrame(absolute_point_values, index=counts_df.index, columns=counts_df.columns)
    lower_ci = pd.DataFrame(lower_ci_values, index=counts_df.index, columns=counts_df.columns)
    upper_ci = pd.DataFrame(upper_ci_values, index=counts_df.index, columns=counts_df.columns)

    # For transparency, provide the posterior mean of the counts
    # This can be seen as the "zero-replaced" or "denoised" counts.
    posterior_mean_counts_matrix = np.zeros_like(counts_matrix)
    for sample_idx, sample in enumerate(counts_df.columns):
        sample_counts = counts_matrix[:, sample_idx]
        total_counts = sample_counts.sum()
        if total_counts > 0:
            dirichlet_params = sample_counts + prior
            post_mean_props = dirichlet_params / dirichlet_params.sum()
            posterior_mean_counts_matrix[:, sample_idx] = post_mean_props * total_counts

    zero_replaced_df = pd.DataFrame(posterior_mean_counts_matrix, index=counts_df.index, columns=counts_df.columns)

    print(f"Confidence intervals computed with proper scaling factor error propagation.")

    if len(scaling_factors) > 1:
        max_sf = max(scaling_factors.values())
        min_sf = min(scaling_factors.values())
        if min_sf > 0:
            print(f"Check: Samples with scaling factors {max_sf:.2f} vs {min_sf:.2f}")
            print(f"should have ~{max_sf/min_sf:.1f}x wider confidence intervals")

    return absolute_point, lower_ci, upper_ci, zero_replaced_df, scaling_factors
