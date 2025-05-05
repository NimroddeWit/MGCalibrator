import pandas as pd

def compute_absolute_abundance(counts_df, dna_conc, volume):
    total_dna = {sample: dna_conc[sample] * volume[sample] for sample in counts_df.columns}
    relative_abundance = counts_df.div(counts_df.sum(axis=0), axis=1)
    absolute = relative_abundance.multiply(pd.Series(total_dna))
    return absolute
