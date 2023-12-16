#!/usr/bin/env julia

"""
Provides `scSigScores` and functions to perform single-cell mean rank signature scoring and
between-group statistical inference
"""
module scRankSigScores

using Statistics: mean, median
using HypothesisTests: pvalue, MannWhitneyUTest, ApproximatePermutationTest
using MultipleTesting: adjust, Bonferroni, BenjaminiHochberg, Hochberg, Holm, Hommel

using SparseArrays: SparseMatrixCSC, sparse, dropzeros!
using HDF5: h5read
using DataFrames: DataFrame, innerjoin

using CairoMakie
using Distances: pairwise, Euclidean
using Clustering: hclust
using Random: seed!, randperm

import CSV

export scSigScores, compute_sig_scores, compute_pvalues, compute_scores_pvalues, calc_fold_changes,
diff_sig_scores, prop_diff_scores, filter_kw, access_result, access_result_k, mwu_test,
scSigScores_pair_from_pf, scSigScores_tri_from_pf, scSigScores_sig_gmt, prefilter,
read_gmt, get_subset_h5ad, get_subset_h5, get_subset_loom, create_filters, write_csv, filter_cells,
write_csv_sig_scores_p, write_csv_p_fc,
heatmap_df

include("sc_sig_scoring.jl")
include("read_write.jl")
include("utilities.jl")
include("plot.jl")

end
