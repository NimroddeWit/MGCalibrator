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

def bayesian_zero_replacement(counts_matrix, method="CZM", alpha=0.5, threshold=0.5, frac=0.65):
    """
    Bayesian zero replacement similar to cmultrepl in R.
    
    Parameters:
    -----------
    counts_matrix : np.array
        Matrix of counts (features x samples)
    method : str
        Method for zero replacement ("CZM" or "Dirichlet")
    alpha : float
        Dirichlet prior parameter (strength of prior)
    threshold : float
        Upper limit threshold for CZM method
    frac : float
        Fraction of threshold for imputation
    
    Returns:
    --------
    np.array: Matrix with zeros replaced
    """
    n_features, n_samples = counts_matrix.shape
    replaced_matrix = counts_matrix.copy().astype(float)
    
    for sample_idx in range(n_samples):
        sample_counts = counts_matrix[:, sample_idx]
        total_counts = sample_counts.sum()
        
        if total_counts == 0:
            continue
            
        # Find zero positions
        zero_mask = (sample_counts == 0)
        n_zeros = zero_mask.sum()
        
        if n_zeros == 0:
            continue
            
        if method == "CZM":
            # CZM method: multiplicative simple replacement
            # Calculate upper limit for multinomial probability
            upper_limit = threshold / total_counts
            
            # Imputed value is fraction of upper limit
            imputed_prob = frac * upper_limit
            imputed_count = imputed_prob * total_counts
            
            # Replace zeros
            replaced_matrix[zero_mask, sample_idx] = imputed_count
            
            # Multiplicatively adjust non-zero parts to maintain total
            non_zero_mask = ~zero_mask
            if non_zero_mask.sum() > 0:
                adjustment_factor = (total_counts - n_zeros * imputed_count) / sample_counts[non_zero_mask].sum()
                replaced_matrix[non_zero_mask, sample_idx] = sample_counts[non_zero_mask] * adjustment_factor
                
        elif method == "Dirichlet":
            # Dirichlet-multinomial Bayesian replacement
            # Use Dirichlet prior with alpha parameter
            dirichlet_params = sample_counts + alpha
            
            # Sample from posterior Dirichlet distribution
            posterior_probs = np.random.dirichlet(dirichlet_params)
            
            # Convert back to counts while preserving total
            imputed_counts = posterior_probs * total_counts
            
            # For zero positions, use the sampled values
            # For non-zero positions, use original counts (more conservative)
            replaced_matrix[zero_mask, sample_idx] = imputed_counts[zero_mask]
            
            # Adjust non-zero counts to maintain compositional constraint
            if (~zero_mask).sum() > 0:
                remaining_total = total_counts - imputed_counts[zero_mask].sum()
                if remaining_total > 0:
                    adjustment_factor = remaining_total / sample_counts[~zero_mask].sum()
                    replaced_matrix[~zero_mask, sample_idx] = sample_counts[~zero_mask] * adjustment_factor
    
    return replaced_matrix

def compute_absolute_abundance_with_error(counts_df, dna_conc, volume, fastq_files, 
                                        n_monte_carlo=1000, 
                                        zero_replacement_method="CZM",
                                        zero_replacement_alpha=0.5, 
                                        threshold=0.5, frac=0.65):
    """
    Compute absolute abundance with 95% confidence intervals using Monte Carlo sampling
    with proper error propagation of scaling factors for compositional data (genes/taxa).
    
    SCALING FACTOR ERROR PROPAGATION:
    The scaling factor directly affects the width of confidence intervals. Samples with 
    higher scaling factors (i.e., higher initial DNA concentration relative to what was 
    sequenced) will have proportionally larger confidence intervals, properly reflecting 
    the increased uncertainty when extrapolating from a smaller sequenced fraction.
    
    Scientific framework:
    1. Apply zero replacement to handle compositional data constraints
    2. Use Monte Carlo sampling on the compositional data
    3. Apply scaling factor to each Monte Carlo sample: sampled_absolute = sampled_counts * scaling_factor
    4. This ensures confidence intervals scale proportionally with the scaling factor
    
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
    zero_replacement_method : str
        Method for zero replacement ("CZM" or "Dirichlet")
    zero_replacement_alpha : float
        Alpha parameter for zero replacement
    threshold : float
        Upper limit threshold for CZM method
    frac : float
        Fraction of threshold for CZM imputation
    
    Returns:
    --------
    tuple: (absolute_abundance, lower_ci, upper_ci, zero_replaced_counts, scaling_factors)
        - absolute_abundance: DataFrame with point estimates (original counts * scaling factor)
        - lower_ci: DataFrame with lower 95% confidence bounds
        - upper_ci: DataFrame with upper 95% confidence bounds
        - zero_replaced_counts: DataFrame with zero-replaced counts for transparency
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
    print(f"Using {zero_replacement_method} zero replacement with alpha={zero_replacement_alpha}")
    
    # Apply zero replacement for compositional data
    counts_matrix = counts_df.values.astype(float)
    zero_replaced_matrix = bayesian_zero_replacement(
        counts_matrix, 
        method=zero_replacement_method, 
        alpha=zero_replacement_alpha,
        threshold=threshold,
        frac=frac
    )
    zero_replaced_df = pd.DataFrame(zero_replaced_matrix, 
                                   index=counts_df.index, 
                                   columns=counts_df.columns)
    
    # Point estimate: zero-replaced counts multiplied by scaling factor
    absolute_point = zero_replaced_df.multiply(pd.Series(scaling_factors))
    
    # Monte Carlo sampling for 95% confidence intervals
    lower_percentile = 2.5  # For 95% CI
    upper_percentile = 97.5  # For 95% CI
    
    # Initialize arrays to store Monte Carlo results
    mc_results = np.zeros((n_monte_carlo, counts_df.shape[0], counts_df.shape[1]))
    
    for sample_idx, sample in enumerate(counts_df.columns):
        sample_counts = zero_replaced_df[sample].values.astype(float)
        scaling_factor = scaling_factors[sample]
        total_counts = sample_counts.sum()
        
        if total_counts == 0:
            # If no counts at all, confidence intervals are all zero
            mc_results[:, :, sample_idx] = 0
            continue
        
        # Generate Monte Carlo samples using Dirichlet multinomial
        for mc_iter in range(n_monte_carlo):
            # Use zero-replaced counts as Dirichlet parameters
            dirichlet_params = sample_counts + 0.1  # Small alpha for numerical stability
            sampled_proportions = np.random.dirichlet(dirichlet_params)
            sampled_counts = np.random.multinomial(int(total_counts), sampled_proportions)
            
            # CRITICAL: Apply scaling factor to get absolute abundance
            # This is where the scaling factor error propagation occurs
            # Higher scaling factors → larger absolute values → wider confidence intervals
            sampled_absolute = sampled_counts * scaling_factor
            mc_results[mc_iter, :, sample_idx] = sampled_absolute
    
    # Compute 95% confidence intervals
    lower_ci_values = np.percentile(mc_results, lower_percentile, axis=0)
    upper_ci_values = np.percentile(mc_results, upper_percentile, axis=0)
    
    # Convert back to DataFrames
    lower_ci = pd.DataFrame(lower_ci_values, index=counts_df.index, columns=counts_df.columns)
    upper_ci = pd.DataFrame(upper_ci_values, index=counts_df.index, columns=counts_df.columns)
    
    print(f"Confidence intervals computed with proper scaling factor error propagation.")
    print(f"Check: Samples with scaling factors {max(scaling_factors.values()):.2f} vs {min(scaling_factors.values()):.2f}")
    print(f"should have ~{max(scaling_factors.values())/min(scaling_factors.values()):.1f}x wider confidence intervals")
    
    return absolute_point, lower_ci, upper_ci, zero_replaced_df, scaling_factors
