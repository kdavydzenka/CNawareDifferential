import os
import pickle as pkl
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import deconveil
from deconveil.utils_processing import *


### Define gene groups ###

pydeseq2_output_path = "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA_test/results/BRCA/res_CNnaive.csv"
deconveil_output_path = "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA_test/results/BRCA/res_CNaware.csv"

# Preprocess PyDESeq2 results
res_pydeseq = process_results(pydeseq2_output_path, "CN naive", lfc_cut = 1.0, pval_cut = 0.05)
res_pydeseq.columns = ["logFC", "padj", "isDE", "DEtype", "method"]

# Preprocess DeConveil results
res_deconveil = process_results(deconveil_output_path, "CN aware", lfc_cut = 1.0, pval_cut = 0.05)
res_deconveil.columns = ["logFC", "padj", "isDE", "DEtype", "method"]

# Join results
res_joint = pd.concat([
    res_pydeseq.add_suffix("_naive"),
    res_deconveil.add_suffix("_aware")
], axis=1)

gene_groups = define_gene_groups(res_joint)





