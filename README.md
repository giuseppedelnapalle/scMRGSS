# scRankSigScores.jl

## Overview

scRankSigScores.jl is a Julia module that provides functionality for mean rank signature scoring and between-group statistical inference on single-cell RNA-seq data. The module includes functions for computing signature scores, p-values, and adjusted p-values, as well as identifying differential signature scores based on fold changes and adjusted p-values.

## Installation

To use the scRankSigScores.jl module, include it in your Julia environment by running:

`using scRankSigScores`

## Dependencies

The module depends on the following Julia packages:

    Statistics: for mean and median calculations.
    SparseArrays: for handling sparse matrices.
    HypothesisTests: for hypothesis testing (Mann-Whitney U test and approximate permutation test).
    MultipleTesting: for multiple comparisons correction.
    HDF5: for reading HDF5 files.
    DataFrames: for working with tabular data.
    CSV: for reading and writing CSV files.

## Usage

Here is a brief overview of the main functionalities provided by the module:

### Signature Scoring

compute_sig_scores: Compute signature scores based on gene ranks across cells for each signature (gene set).

`compute_sig_scores(object::scSigScores; kw...)`

### Hypothesis Testing

compute_pvalues: Compute p-values for the signature scores over signatures.

`compute_pvalues(object::scSigScores, ref_object::scSigScores, proc::String; kw...)`

### Differential Signature Scores

diff_sig_scores: Identify differential signature scores based on fold changes and adjusted p-values.

`diff_sig_scores(object::scSigScores, object_2::scSigScores; cut_off::Float64=1.2, adj_p_thresh::Float64=0.01, filter_scores::Bool=true, modify_obj::Bool=false)`

### Fold Change Calculation

calc_fold_changes: Calculate fold changes of signature scores between two objects.

`calc_fold_changes(object::scSigScores, object_2::scSigScores; filter_scores::Bool=true, modify_obj::Bool=false)`

### Accessing Results

access_result: Access the result using the key of the results field.

`access_result(object::scSigScores, key::AbstractString) -> DataFrame`

## Example

```julia
using scRankSigScores

# Load your single-cell RNA-seq data
data = ...

# Load your reference single-cell RNA-seq data
ref_data = ...

# Define signature genes
signature_genes = [...]

# Create scSigScores objects
sc_object = scSigScores("dataset1", cell_annot, data, genes, Dict("Signature1" => signature_genes))
ref_object = scSigScores("ref_dataset", ref_cell_annot, ref_data, ref_genes, Dict("Signature1" => signature_genes))

# Compute signature scores and adjust for a covariate
compute_adj_sig_scores_p(sc_object, ref_object, "n_counts")

# Access the results
sig_scores_df = access_result(sc_object, "sig_scores")
pvalues_df = access_result(sc_object, "two_sample_test")

# Identify differential signature scores
diff_scores_df = diff_sig_scores(sc_object, ref_object)
```

For more details, refer to the function documentation and examples provided in the module.

## Contributors

    Nikola (Nov 2023)

## Version

    1.6.6

## Notes

    The module does not adjust scores for covariates.
    If the signature test fails, scores are set to missing.

## License

This module is distributed under the MIT License.