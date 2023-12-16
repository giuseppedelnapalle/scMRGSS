#!/usr/bin/env julia
#= 
functions to perform comparison of signature scores of Jurkat or 293T data
nikola
Jan 2023
version 1.2
=#

# import modules and functions
using Base.Threads: @threads
using Combinatorics: combinations
using DataFrames: DataFrame

const dir_pkg = "/home/nikola/Project_Data/Julia_data/packages/scRankSigScores/src"
push!(LOAD_PATH, dir_pkg)
using scRankSigScores

# jurkat
"""
    Read data from jurkat dataset

Parameters
----------
...
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell annotations
"""
function read_jurkat(label::String="jurkat",
    filename_dt::String="/home/nikola/Project_Data/R_data/Jurkat_293T/Jurkat/output/mat_jurkat.h5ad",
    dataset_id::String="jurkat", path_exp_dataset::String="X/data",
    path_gene_names::String="var/_index", path_cell_annot::String="obs/label/codes";
    full_path_cell_cat::String="obs/label/categories", kw...)
    f = path_cell_annot => ["jurkat"]
    filters = create_filters(f)
    mat, genes, cell_ann = read_jk_293t_h5ad(filename_dt, dataset_id, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    full_path_cell_cat=full_path_cell_cat, kw...)
    return (mat, genes, cell_ann) 
end

function read_jk_293t_h5ad(filename_dt::String, dataset_id::String,
    path_exp_dataset::String, path_gene_names::String, path_cell_annot::String, filters::Dict{String, Vector{String}};
    kw...)
    kw_gsh = filter_kw(get_subset_h5ad; kw...)
    (mat, genes, cell_ann) = get_subset_h5ad(filename_dt, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    kw_gsh...)
    kw_pf = filter_kw(prefilter; kw...)
    mat, genes, cell_ann = prefilter(mat, genes, cell_ann; kw_pf...)
    return (mat, genes, cell_ann) 
end

# 293t
"""
    Read data from 293t dataset

Parameters
----------
...
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell annotations
"""
function read_293t(label::String="293t",
    filename_dt::String="/home/nikola/Project_Data/R_data/Jurkat_293T/293T/output/mat_293t.h5ad",
    dataset_id::String="293t", path_exp_dataset::String="X/data",
    path_gene_names::String="var/_index", path_cell_annot::String="obs/label/codes";
    full_path_cell_cat::String="obs/label/categories", kw...)
    f = path_cell_annot => ["293t"]
    filters = create_filters(f)
    mat, genes, cell_ann = read_jk_293t_h5ad(filename_dt, dataset_id, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    full_path_cell_cat=full_path_cell_cat, kw...)
    return (mat, genes, cell_ann) 
end

# jurkat_293t_50_50
"""
    Read data from jurkat_293t_50_50 dataset

Parameters
----------
...
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell annotations
"""
function read_jurkat_293t_50_50(label::String,
    filename_dt::String="/home/nikola/Project_Data/R_data/Jurkat_293T/Jurkat_293T_50_50/output/mat_jurkat_293t_50_50.h5ad",
    dataset_id::String="jurkat_293t_50_50", path_exp_dataset::String="X/data",
    path_gene_names::String="var/_index", path_cell_annot::String="obs/label/codes";
    full_path_cell_cat::String="obs/label/categories", kw...)
    f = path_cell_annot => [label]
    filters = create_filters(f)
    mat, genes, cell_ann = read_jk_293t_h5ad(filename_dt, dataset_id, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    full_path_cell_cat=full_path_cell_cat, kw...)
    return (mat, genes, cell_ann) 
end

# jurkat_293t_99_1
"""
    Read data from jurkat_293t_99_1 dataset

Parameters
----------
...
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell annotations
"""
function read_jurkat_293t_99_1(label::String="jurkat",
    filename_dt::String="/home/nikola/Project_Data/R_data/Jurkat_293T/Jurkat_293T_99_1/output/mat_jurkat_293t_99_1.h5ad",
    dataset_id::String="jurkat_293t_99_1", path_exp_dataset::String="X/data",
    path_gene_names::String="var/_index", path_cell_annot::String="obs/label/codes";
    full_path_cell_cat::String="obs/label/categories", kw...)
    f = path_cell_annot => ["jurkat"]
    filters = create_filters(f)
    mat, genes, cell_ann = read_jk_293t_h5ad(filename_dt, dataset_id, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    full_path_cell_cat=full_path_cell_cat, kw...)
    return (mat, genes, cell_ann) 
end

"""
    Read data from Jurkat, 293T, or mixture dataset

Parameters
----------
label
    Cell label (default values set for 293t, jurkat, and jurkat_293t_99_1 datasets)
...
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell annotations
"""
function read_jk_293t(label::String, dataset_id::String;
    kw...)
    if dataset_id=="jurkat"
        mat, genes, cell_ann = read_jurkat("jurkat"; kw...)
    elseif dataset_id=="293t"
        mat, genes, cell_ann = read_293t("293t"; kw...)
    elseif dataset_id=="jurkat_293t_50_50"
        mat, genes, cell_ann = read_jurkat_293t_50_50(label; kw...)
    elseif dataset_id=="jurkat_293t_99_1"
        mat, genes, cell_ann = read_jurkat_293t_99_1("jurkat"; kw...)
    else
        throw(ArgumentError(("""
        Incorrect argument for dataset_id. Please choose from "jurkat", "293t",
        "jurkat_293t_50_50", or "jurkat_293t_99_1".
        """)))
    end
    return (mat, genes, cell_ann)
end

"""
    Workflow to calculate proportions of differential signature scores over combinations of cell labels and datasets.
"""
function prop_diff_scores_wf(labels::Vector{String}, dataset_ids::Vector{String}, filename_sig::String,
    dataset_dict::Dict{String, Vector{String}};
    min_features::Int64=2000, min_cells::Int64=2, sig_thresh::Float64=.6, tail::String="both",
    fc_thresh::Float64=1.15, adj_p_thresh::Float64=.01, n_c_thresh::Int64=100)
    println("Cell labels: ", labels)
    println("Dataset IDs: ", dataset_ids)

    lb_ds = [(i, j) for i in labels for j in dataset_ids if i in dataset_dict[j]]
    r_c_nm = [join([lb_ds[i][1], lb_ds[i][2]], "-") for i in eachindex(lb_ds)]
    comb = collect(combinations(lb_ds, 2)) # label_dataset pairs
    len = length(r_c_nm)

    # coordinate dict
    r = vcat([repeat([i], len-i) for i in 1:(len-1)]...) # Vector{Vector{Int64}} collapsed to Vector{Int64}
    c = vcat([[(j+1):len;] for j in 1:(len-1)]...)
    ind_p = collect(zip(r, c))
    pos = Dict(zip(comb, ind_p)) # {Vector{Tuple{String, String}}, Tuple{Int64, Int64}}

    prop = Matrix{Union{Float64, Missing}}(undef, len, len)
    # upper triangular matrix excluding diagonal elements
    println("Begin calculating proportions of differential signature scores...")
    @threads for p in comb
        label, ds_id = p[1]
        label_2, ds_id_2 = p[2]
        i, j = pos[p]
        key_s = join([label, ds_id], "-")
        key_s_2 = join([label_2, ds_id_2], "-")

        mat, genes, cell_ann = read_jk_293t(label, ds_id; min_features=min_features, min_cells=min_cells)
        mat_2, genes_2, cell_ann_2 = read_jk_293t(label_2, ds_id_2; min_features=min_features, min_cells=min_cells)
        
        if size(mat)[2] < n_c_thresh || size(mat_2)[2] < n_c_thresh
            if size(mat)[2] < n_c_thresh
                print("""
                # cells with label "$label" less than $n_c_thresh in $(ds_id). Proportion of differential signature scores coerced to `missing`.
                """)
            end
            if size(mat_2)[2] < n_c_thresh
                print("""
                # cells with label "$label_2" less than $n_c_thresh in $(ds_id_2). Proportion of differential signature scores coerced to `missing`.
                """)
            end
            prop[i, j] = missing
        else
            obj, obj_2 = scSigScores_pair_from_pf(ds_id, mat, genes, cell_ann, ds_id_2, mat_2, genes_2, cell_ann_2, filename_sig)
            # compute_adj_sig_scores_p(obj, obj_2; sig_thresh=sig_thresh)
            compute_sig_scores(obj; sig_thresh=sig_thresh)
            compute_sig_scores(obj_2; sig_thresh=sig_thresh)
            compute_pvalues(obj, obj_2, "Mann_Whitney_U"; tail=tail)
            diff_scores = diff_sig_scores(obj, obj_2; cut_off=fc_thresh, adj_p_thresh=adj_p_thresh)
            prop[i, j] = prop_diff_scores(diff_scores)
        end
    end

    # complete symmetric elements
    for j in 1:(len-1)
        for i in (j+1):len
            prop[i, j] = prop[j, i]
        end
    end

    # complete diagonal elements
    for i in 1:len
        prop[i, i] = 0.0
    end
    
    df_p = hcat(DataFrame(row_name=r_c_nm), DataFrame(prop, r_c_nm))
end

"""
    Workflow to calculate proportions of differential signature scores over combinations of cell labels and datasets
    and export results of signature scores, p-values, and fold changes.
"""
function prop_diff_scores_export_res_wf(labels::Vector{String}, dataset_ids::Vector{String}, filename_sig::String,
    dataset_dict::Dict{String, Vector{String}}, res_dir::String, gene_set_id::T;
    min_features::Int64=2000, min_cells::Int64=2, sig_thresh::Float64=.6, tail::String="both",
    fc_thresh::Float64=1.15, adj_p_thresh::Float64=.01, n_c_thresh::Int64=100) where {T<:AbstractString}
    println("Cell labels: ", labels)
    println("Dataset IDs: ", dataset_ids)

    lb_ds = [(i, j) for i in labels for j in dataset_ids if i in dataset_dict[j]]
    r_c_nm = [join([lb_ds[i][1], lb_ds[i][2]], "-") for i in eachindex(lb_ds)]
    comb = collect(combinations(lb_ds, 2)) # label_dataset pairs
    len = length(r_c_nm)

    # coordinate dict
    r = vcat([repeat([i], len-i) for i in 1:(len-1)]...) # Vector{Vector{Int64}} collapsed to Vector{Int64}
    c = vcat([[(j+1):len;] for j in 1:(len-1)]...)
    ind_p = collect(zip(r, c))
    pos = Dict(zip(comb, ind_p)) # {Vector{Tuple{String, String}}, Tuple{Int64, Int64}}

    prop = Matrix{Union{Float64, Missing}}(undef, len, len)
    # upper triangular matrix excluding diagonal elements
    println("Begin calculating proportions of differential signature scores...")
    @threads for p in comb
        label, ds_id = p[1]
        label_2, ds_id_2 = p[2]
        i, j = pos[p]
        key_s = join([label, ds_id], "-")
        key_s_2 = join([label_2, ds_id_2], "-")

        mat, genes, cell_ann = read_jk_293t(label, ds_id; min_features=min_features, min_cells=min_cells)
        mat_2, genes_2, cell_ann_2 = read_jk_293t(label_2, ds_id_2; min_features=min_features, min_cells=min_cells)
        
        if size(mat)[2] < n_c_thresh || size(mat_2)[2] < n_c_thresh
            if size(mat)[2] < n_c_thresh
                print("""
                # cells with label "$label" less than $n_c_thresh in $(ds_id). Proportion of differential signature scores coerced to `missing`.
                """)
            end
            if size(mat_2)[2] < n_c_thresh
                print("""
                # cells with label "$label_2" less than $n_c_thresh in $(ds_id_2). Proportion of differential signature scores coerced to `missing`.
                """)
            end
            prop[i, j] = missing
        else
            obj, obj_2 = scSigScores_pair_from_pf(ds_id, mat, genes, cell_ann, ds_id_2, mat_2, genes_2, cell_ann_2, filename_sig)
            # compute_adj_sig_scores_p(obj, obj_2; sig_thresh=sig_thresh)
            compute_sig_scores(obj; sig_thresh=sig_thresh)
            compute_sig_scores(obj_2; sig_thresh=sig_thresh)
            compute_pvalues(obj, obj_2, "Mann_Whitney_U"; tail=tail)
            diff_scores = diff_sig_scores(obj, obj_2; cut_off=fc_thresh, adj_p_thresh=adj_p_thresh)
            prop[i, j] = prop_diff_scores(diff_scores)

            calc_fold_changes(obj, obj_2; modify_obj=true)
            write_csv_sig_scores_p(obj, ds_id, label, obj_2, ds_id_2, label_2, res_dir, gene_set_id)
            write_csv_p_fc(obj, ds_id, label, ds_id_2, label_2, res_dir, gene_set_id)
        end
    end

    # complete symmetric elements
    for j in 1:(len-1)
        for i in (j+1):len
            prop[i, j] = prop[j, i]
        end
    end

    # complete diagonal elements
    for i in 1:len
        prop[i, i] = 0.0
    end
    
    df_p = hcat(DataFrame(row_name=r_c_nm), DataFrame(prop, r_c_nm))
end
