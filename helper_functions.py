def process_experiment(df_full, experiment, protease='trypsin', max_missed_cleavages=2):
    """Process a single experiment"""
    # Filter
    mask = create_combined_mask(df_full, experiment['file_tags'])
    _df = df_full[mask]
    
    # Calculate CV statistics
    _df_pg_cv = _df.groupby('Protein.Group')['PG.MaxLFQ'].agg(['mean', 'std'])
    _df_pg_cv['cv'] = _df_pg_cv['std'] / _df_pg_cv['mean']
    
    _df_pr_cv = _df.groupby('Precursor.Id')['Precursor.Normalised'].agg(['mean', 'std'])
    _df_pr_cv['cv'] = _df_pr_cv['std'] / _df_pr_cv['mean']
    
    # Calculate missed cleavages PER RUN
    def calculate_mc_per_run(group):
        """Calculate missed cleavage distribution for a single run"""
        # Get unique peptides (drop duplicates, keep first occurrence)
        unique_peptides = group['Stripped.Sequence'].drop_duplicates(keep='first')
        
        mc_counts = {}
        
        for seq in unique_peptides:
            mc = count_missed_cleavages(seq, protease)
            mc_capped = min(mc, max_missed_cleavages)
            mc_counts[mc_capped] = mc_counts.get(mc_capped, 0) + 1
        
        total = len(unique_peptides)
        
        # Create result dict with relative proportions
        result = {}
        for i in range(max_missed_cleavages + 1):
            count = mc_counts.get(i, 0)
            result[f'MC{i}'] = count / total if total > 0 else 0
        
        return pd.Series(result)
    
    # Apply per run
    mc_per_run = _df.groupby('Run').apply(calculate_mc_per_run).reset_index()
    
    # Aggregate per Run
    _df_agg = _df.groupby('Run', as_index=False).agg({
        'Modified.Sequence': 'nunique',
        'Precursor.Id': 'nunique',
        'Protein.Group': 'nunique'
    }).rename(columns={
        'Modified.Sequence': 'peptide',
        'Precursor.Id': 'precursor',
        'Protein.Group': 'protein'
    })
    
    # Merge missed cleavage data
    _df_agg = _df_agg.merge(mc_per_run, on='Run', how='left')
    
    # Calculate average missed cleavages (weighted average)
    avg_mc = 0
    for i in range(max_missed_cleavages + 1):
        avg_mc += _df_agg[f'MC{i}'] * i
    _df_agg['avg_MC'] = avg_mc
    
    # Add statistics
    _df_agg = _df_agg.assign(
        PG20=(_df_pg_cv['cv'] < 0.2).sum(),
        Pr20=(_df_pr_cv['cv'] < 0.2).sum(),
        total_peptides=_df['Stripped.Sequence'].nunique(),
        total_protein_groups=_df['Protein.Group'].nunique(),
        total_precursors=_df['Precursor.Id'].nunique(),
        instrument=experiment['instrument'],
        method=experiment['method']
    )
    
    return _df_agg

def calculate_mc_per_run(group):
        """Calculate missed cleavage distribution for a single run"""
        # Get unique peptides (drop duplicates, keep first occurrence)
        unique_peptides = group['Stripped.Sequence'].drop_duplicates(keep='first')
        
        mc_counts = {}
        
        for seq in unique_peptides:
            mc = count_missed_cleavages(seq, protease)
            mc_capped = min(mc, max_missed_cleavages)
            mc_counts[mc_capped] = mc_counts.get(mc_capped, 0) + 1
        
        total = len(unique_peptides)
        
        # Create result dict with relative proportions
        result = {}
        for i in range(max_missed_cleavages + 1):
            count = mc_counts.get(i, 0)
            result[f'MC{i}'] = count / total if total > 0 else 0
        
        return pd.Series(result)

import pandas as pd
import re
from functools import lru_cache
from collections import defaultdict
@lru_cache(maxsize=None)
def load_parquet_cached(path):
    df = pd.read_parquet(path, 
                        columns=['Run', 'PG.Q.Value', 'PG.MaxLFQ', 
                                'Precursor.Normalised', 'Precursor.Id',
                                'Protein.Group', 'Stripped.Sequence','Modified.Sequence','Genes'])
    return df


def generate_pattern_list(prefixes, start, end, prefix_all=''):
    """
    Generate a list of strings following a pattern with optional prefix for all elements.
    
    Parameters:
    prefixes (str or list): Single prefix or list of prefixes (e.g., 'A' or ['A', 'B', 'C'])
    start (int): The starting number
    end (int): The ending number (inclusive)
    prefix_all (str): Optional prefix to add before every element (default: '')
    
    Returns:
    list: A list of strings following the pattern
    
    Examples:
    >>> generate_pattern_list(['A', 'B', 'C'], 1, 4, 'SPEC_')
    ['SPEC_A1', 'SPEC_A2', 'SPEC_A3', 'SPEC_A4', 'SPEC_B1', 'SPEC_B2', 'SPEC_B3', 'SPEC_B4', 'SPEC_C1', 'SPEC_C2', 'SPEC_C3', 'SPEC_C4']
    """
    if isinstance(prefixes, str):
        prefixes = [prefixes]
    
    return [f"{prefix_all}{prefix}{i}" for prefix in prefixes for i in range(start, end + 1)]
