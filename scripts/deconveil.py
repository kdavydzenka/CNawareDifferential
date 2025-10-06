import os
import pickle as pkl
import pandas as pd
import numpy as np

import deconveil
from deconveil.dds import deconveil_fit
from deconveil.inference import Inference
from deconveil.default_inference import DefInference
from deconveil.utils_fit import *
from deconveil.utils_processing import *
from deconveil import deconveil_fit
from deconveil.ds import deconveil_stats


def run_deconveil(rna_counts, metadata, cnv, output_path, design_factors="condition", alpha=0.05):
    """
    Runs DeConveil analysis and saves the results.

    Parameters:
        rna_counts (pd.DataFrame): Count matrix with genes as rows and samples as columns.
        metadata (pd.DataFrame): Metadata for the samples with design factors.
        cnv (pd.DataFrame): Copy number variation (CNV) data matrix  with genes as rows and samples as columns.
        output_path (str): Directory to save the results.
        design_factors (str): Column in metadata to use for design.
        alpha (float): Significance level for statistical tests.
    """
    os.makedirs(output_path, exist_ok=True)
    
    # Initialize DeConveil inference
    inference = DefInference(n_cpus=8)

    # Fit DeConveil model
    dds = deconveil_fit(
        counts=rna_counts,
        metadata=metadata,
        cnv=cnv,
        design_factors=design_factors,
        inference=inference,
        refit_cooks=True
    )
    dds.fit_size_factors()
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()

    if dds.refit_cooks:
        dds.refit()  # Replace outlier counts

    # Statistical analysis
    stat_res_deconveil = deconveil_stats(
        dds, 
        alpha=alpha, 
        independent_filter=True, 
        cooks_filter=True
    )
    stat_res_deconveil.run_wald_test()

    if stat_res_deconveil.independent_filter:
        stat_res_deconveil._independent_filtering()
    else:
        stat_res_deconveil._p_value_adjustment()

    # Log-fold change shrinkage
    stat_res_deconveil.lfc_shrink(coeff="condition_B_vs_A")
    stat_res_deconveil.summary()

    # Save results
    results_path = os.path.join(output_path, "res_CNaware.csv")
    stat_res_deconveil.results_df.to_csv(results_path)
    return(stat_res_deconveil.results_df)

### Load data ###

DATA_PATH = "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA_test/data/BRCA/deconveil_input"
rna_counts = pd.read_csv(os.path.join(DATA_PATH, "rna.csv"), index_col=0)
rna_counts = rna_counts.T
metadata = pd.read_csv(os.path.join(DATA_PATH, "metadata.csv"), index_col=0)
cnv = pd.read_csv(os.path.join(DATA_PATH, "cnv.csv"), index_col=0)
cnv = cnv.T

## Reorder cnv to match rna_counts ##

cnv = cnv.loc[rna_counts.index]
assert (rna_counts.index == metadata.index).all(), "Sample order mismatch between rna_counts and metadata"
assert (cnv.index == rna_counts.index).all(), "Sample order mismatch between cnv and rna_counts"

## Remove all zero genes ##

all_zero_mask = (rna_counts.sum(axis=0) == 0)
print("All-zero genes:", all_zero_mask.sum())
rna_counts = rna_counts.loc[:, ~all_zero_mask]
cnv = cnv.loc[:, ~all_zero_mask]  

# Filter RNA + CNV together

res = filter_low_count_genes(rna_counts, other_dfs=[cnv], min_count=300, min_samples=50)
rna_counts = res["filtered_df"]
cnv = res["other_filtered"][0]
print("After low-count filtering:", rna_counts.shape, cnv.shape)


### Run DeConveil ###

deconveil_output_path = "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA_test/results/BRCA/"
run_deconveil(rna_counts, metadata, cnv, deconveil_output_path)  
