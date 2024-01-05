#!/usr/bin/env julia
#= 
functions to perform comparison of signature scores of PBMC data
nikola
Jan 2023
version 1.1
=#

# import modules and functions
using Base.Threads: @threads
import CSV
using DataFrames: DataFrame, innerjoin
# using Base.Threads: @spawn

const dir_pkg = "/home/nikola/Project_Data/Julia_data/packages/scRankSigScores/src"
push!(LOAD_PATH, dir_pkg)
using scRankSigScores

# pbmc_10k_3p
"""
    Read data from pbmc_10k_3p dataset

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
function read_pbmc_10k_3p(label::String,
    filename_dt::String="/home/nikola/Documents/Rsch/resources/scRNA_seq/blood/PBMC_10k/3p_pbmc10k_filt.h5",
    fn_cell_annot::String="/home/nikola/Project_Data/R_data/pbmc_10k_integration/output/cell_annotations/pbmc_10k_3p/c_ann_pbmc10k_3p.csv",
    path_gene_names::String="matrix/features/name", 
    path_exp_dataset::String="matrix/data", path_cell_barcodes::String="matrix/barcodes";
    kw...)
    mat, genes, cell_ann = read_pbmc_10k(filename_dt, fn_cell_annot, label, 
    path_gene_names, path_exp_dataset, path_cell_barcodes; kw...)
    return (mat, genes, cell_ann)
end

function read_pbmc_10k(filename_dt::String, fn_cell_annot::String, label::String,
    path_gene_names::String="matrix/features/name", 
    path_exp_dataset::String="matrix/data", path_cell_barcodes::String="matrix/barcodes";
    kw...)
    kw_gsh = filter_kw(get_subset_h5; kw...)
    (mat, genes, cell_ann) = get_subset_h5(filename_dt, fn_cell_annot, [label], path_gene_names, 
    path_exp_dataset, path_cell_barcodes; kw_gsh...)
    kw_pf = filter_kw(prefilter; kw...)
    (mat, genes, cell_ann) = prefilter(mat, genes, cell_ann; kw_pf...)
    return (mat, genes, cell_ann)
end

# pbmc_10k_5p
"""
    Read data from pbmc_10k_5p dataset

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
function read_pbmc_10k_5p(label::String,
    filename_dt::String="/home/nikola/Documents/Rsch/resources/scRNA_seq/blood/PBMC_10k/5p_pbmc10k_filt.h5",
    fn_cell_annot::String="/home/nikola/Project_Data/R_data/pbmc_10k_integration/output/cell_annotations/pbmc_10k_5p/c_ann_pbmc10k_5p.csv",
    path_gene_names::String="matrix/features/name", path_exp_dataset::String="matrix/data",
    path_cell_barcodes::String="matrix/barcodes";
    kw...)
    mat, genes, cell_ann = read_pbmc_10k(filename_dt, fn_cell_annot, label, 
    path_gene_names, path_exp_dataset, path_cell_barcodes; kw...)
    return (mat, genes, cell_ann)
end

# Tabula_Sapiens
const c_ann_ts = Dict{String, Vector{String}}(
    "mp" => ["macrophage"],
    "cd14_mono" => ["classical monocyte"],
    "cd16_mono" => ["non-classical monocyte"],
    "neu" => ["nampt neutrophil", "neutrophil"],
    "nk" => ["nk cell"],
    "cd8_t" => ["cd8-positive, alpha-beta t cell"],
    "cd4_t" => ["cd4-positive, alpha-beta t cell"],
    "naive_cd4_t" => ["naive thymus-derived cd4-positive, alpha-beta t cell"],
    "mem_b" => ["memory b cell"],
    "naive_b" => ["naive b cell"],
    "plsm" => ["plasma cell"],
    "b" => ["memory b cell", "naive b cell", "plasma cell"],
    "hsc" => ["hematopoietic stem cell"]
)

"""
    Read data from the TabulaSapiens dataset

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
function read_blood_ts(label::String,
    filename_dt::String="/home/nikola/Documents/Rsch/resources/omics_data_projects/HCA/Tabula_Sapiens/version_5/TabulaSapiens.h5ad",
    path_exp_dataset::String="layers/raw_counts/data", path_gene_names::String="var/gene_symbol/codes",
    path_cell_annot::String="obs/cell_ontology_class/codes", path_cell_grp::String="obs/compartment/codes";
    path_organ_orig::String="obs/organ_tissue/codes", full_path_gene_cat::String="var/gene_symbol/categories",
    full_path_cell_cat::String="obs/cell_ontology_class/categories", kw...)
    f = path_cell_annot => c_ann_ts[label]
    f_2 = path_cell_grp => ["immune"]
    f_3 = path_organ_orig => ["Blood"]
    filters = create_filters(f, f_2, f_3)
    mat, genes, cell_ann = read_blood_h5ad(filename_dt, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    full_path_gene_cat=full_path_gene_cat, full_path_cell_cat=full_path_cell_cat, kw...)
    return (mat, genes, cell_ann) 
end

function read_blood_h5ad(filename_dt::String, path_exp_dataset::String,
    path_gene_names::String, path_cell_annot::String, filters::Dict{String, Vector{String}};
    kw...)
    kw_gsh = filter_kw(get_subset_h5ad; kw...)
    (mat, genes, cell_ann) = get_subset_h5ad(filename_dt, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    kw_gsh...)
    kw_pf = filter_kw(prefilter; kw...)
    mat, genes, cell_ann = prefilter(mat, genes, cell_ann; kw_pf...)
    return (mat, genes, cell_ann) 
end

# HematopoieticImmuneCellAtlas
const c_ann_hica = Dict{String, Vector{String}}(
    "cd14_mono" => ["CD14+ monocytes"],
    "cd16_mono" => ["CD16+ monocytes"],
    "nk" => ["NK cells"],
    "cd8_t" => ["Cytotoxic T cells"],
    "cd4_t" => ["T helper cells"],
    "naive_cd4_t" => ["CD4+ naive T cells"],
    "mem_b" => ["Memory B cells"],
    "naive_b" => ["Naive B cells"],
    "plsm" => ["Plasma cells"],
    "b" => ["Memory B cells", "Naive B cells", "Plasma cells"]
)

"""
    Read data from the HematopoieticImmuneCellAtlas dataset

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
function read_blood_hica(label::String,
    filename_dt::String="/home/nikola/Documents/Rsch/resources/HCA/blood/HematopoieticImmuneCellAtlas/BL_standard_design.h5ad",
    path_exp_dataset::String="raw/X/data", path_gene_names::String="var/featurekey", path_cell_annot::String="obs/anno";
    kw...)
    f = path_cell_annot => c_ann_hica[label]
    filters = create_filters(f)
    mat, genes, cell_ann = read_blood_h5ad(filename_dt, path_exp_dataset, path_gene_names, path_cell_annot, filters;
    kw...)
    return (mat, genes, cell_ann) 
end

"""
    Read data from pbmc or blood dataset

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
function read_pbmc(label::String, dataset_id::String;
    kw...)
    if dataset_id=="pbmc_10k_3p"
        mat, genes, cell_ann = read_pbmc_10k_3p(label; kw...)
    elseif dataset_id=="pbmc_10k_5p"
        mat, genes, cell_ann = read_pbmc_10k_5p(label; kw...)
    elseif dataset_id=="TabulaSapiens"
        mat, genes, cell_ann = read_blood_ts(label; kw...)
    elseif dataset_id=="HematopoieticImmuneCellAtlas"
        mat, genes, cell_ann = read_blood_hica(label; kw...)
    else
        throw(ArgumentError(("""
        Incorrect argument for dataset_id. Please choose from "pbmc_10k_3p", "pbmc_10k_5p", 
        "TabulaSapiens", or "HematopoieticImmuneCellAtlas".
        """)))
    end
    return (mat, genes, cell_ann)
end

"""
    Workflow to calculate proportion of differential signature scores over combinations of cell labels and datasets,
    and export the results of signature scores, p-values, and fold changes.
"""
function prop_diff_scores_export_res_wf(labels::Vector{String}, dataset_ids::Vector{String},
    filename_sig::String, res_dir::String, gene_set_id::T;
    min_features::Int64=2000, min_cells::Int64=2, sig_thresh::Float64=.6, tail::String="both",
    fc_thresh::Float64=1.15, adj_p_thresh::Float64=.01, n_c_thresh::Int64=100) where {T<:AbstractString}
    println("Cell labels: ", labels)
    println("Dataset IDs: ", dataset_ids)
    n_lb = length(labels)
    n_ds = length(dataset_ids)
    len = n_lb * n_ds
    lb_ds = [(i, j) for i in labels for j in dataset_ids]
    r_c_nm = [join([lb_ds[i][1], lb_ds[i][2]], "-") for i in eachindex(lb_ds)]

    prop = Matrix{Union{Float64, Missing}}(undef, len, len)
    r_ind = Vector{Vector{Int32}}(undef, len-1)
    c_ind = Vector{Vector{Int32}}(undef, len-1)
    arg_p = Vector{Vector{Vector{Tuple{String, String}}}}(undef, len-1)

    # lower triangular matrix excluding diagonal elements
    println("Begin calculating proportions of differential signature scores...")
    for (c, r) in enumerate([2:len;])
        r_ind[c] = repeat([r], c)
        c_ind[c] = [1:c;]
        arg_p[c] = [[lb_ds[i], lb_ds[j]] for i in r for j in 1:c]
    end

    @threads for i in eachindex(arg_p) # @threads
        @threads for j in eachindex(arg_p[i])
            label, ds_id = arg_p[i][j][1]
            label_2, ds_id_2 = arg_p[i][j][2]
            m = r_ind[i][j]
            n = c_ind[i][j]

            mat, genes, cell_ann = read_pbmc(label, ds_id; min_features=min_features, min_cells=min_cells)
            mat_2, genes_2, cell_ann_2 = read_pbmc(label_2, ds_id_2; min_features=min_features, min_cells=min_cells)

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
                prop[m, n] = missing
            else
                obj, obj_2 = scSigScores_pair_from_pf(ds_id, mat, genes, cell_ann, ds_id_2, mat_2, genes_2, cell_ann_2, filename_sig)
                compute_sig_scores(obj; sig_thresh=sig_thresh)
                compute_sig_scores(obj_2; sig_thresh=sig_thresh)
                compute_pvalues(obj, obj_2, "Mann_Whitney_U"; tail=tail)
                diff_scores = diff_sig_scores(obj, obj_2; cut_off=fc_thresh, adj_p_thresh=adj_p_thresh)
                prop[m, n] = prop_diff_scores(diff_scores)

                calc_fold_changes(obj, obj_2; modify_obj=true)
                write_csv_sig_scores_p(obj, ds_id, label, obj_2, ds_id_2, label_2, res_dir, gene_set_id)
                write_csv_p_fc(obj, ds_id, label, ds_id_2, label_2, res_dir, gene_set_id)
            end
        end
    end
    
    # complete symmetric elements
    for r in 1:(len-1)
        for c in (r+1):len
            prop[r, c] = prop[c, r]
        end
    end

    # complete diagonal elements
    for i in 1:len
        prop[i, i] = 0.0
    end

    df_p = hcat(DataFrame(row_name=r_c_nm), DataFrame(prop, r_c_nm))
end

"""
    Workflow to calculate proportion of differential signature scores over combinations of cell labels and datasets.
"""
function prop_diff_scores_wf(labels::Vector{String}, dataset_ids::Vector{String}, filename_sig::String;
    min_features::Int64=2000, min_cells::Int64=2, sig_thresh::Float64=.6, tail::String="both",
    fc_thresh::Float64=1.15, adj_p_thresh::Float64=.01, n_c_thresh::Int64=100)
    println("Cell labels: ", labels)
    println("Dataset IDs: ", dataset_ids)
    n_lb = length(labels)
    n_ds = length(dataset_ids)
    len = n_lb * n_ds
    lb_ds = [(i, j) for i in labels for j in dataset_ids]
    r_c_nm = [join([lb_ds[i][1], lb_ds[i][2]], "-") for i in eachindex(lb_ds)]

    prop = Matrix{Union{Float64, Missing}}(undef, len, len)
    r_ind = Vector{Vector{Int32}}(undef, len-1)
    c_ind = Vector{Vector{Int32}}(undef, len-1)
    arg_p = Vector{Vector{Vector{Tuple{String, String}}}}(undef, len-1)

    # lower triangular matrix excluding diagonal elements
    println("Begin calculating proportions of differential signature scores...")
    for (c, r) in enumerate([2:len;])
        r_ind[c] = repeat([r], c)
        c_ind[c] = [1:c;]
        arg_p[c] = [[lb_ds[i], lb_ds[j]] for i in r for j in 1:c]
    end

    @threads for i in eachindex(arg_p) # @threads
        @threads for j in eachindex(arg_p[i])
            label, ds_id = arg_p[i][j][1]
            label_2, ds_id_2 = arg_p[i][j][2]
            m = r_ind[i][j]
            n = c_ind[i][j]

            mat, genes, cell_ann = read_pbmc(label, ds_id; min_features=min_features, min_cells=min_cells)
            mat_2, genes_2, cell_ann_2 = read_pbmc(label_2, ds_id_2; min_features=min_features, min_cells=min_cells)

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
                prop[m, n] = missing
            else
                obj, obj_2 = scSigScores_pair_from_pf(ds_id, mat, genes, cell_ann, ds_id_2, mat_2, genes_2, cell_ann_2, filename_sig)
                compute_sig_scores(obj; sig_thresh=sig_thresh)
                compute_sig_scores(obj_2; sig_thresh=sig_thresh)
                compute_pvalues(obj, obj_2, "Mann_Whitney_U"; tail=tail)
                diff_scores = diff_sig_scores(obj, obj_2; cut_off=fc_thresh, adj_p_thresh=adj_p_thresh)
                prop[m, n] = prop_diff_scores(diff_scores)
            end
        end
    end
    
    # complete symmetric elements
    for r in 1:(len-1)
        for c in (r+1):len
            prop[r, c] = prop[c, r]
        end
    end

    # complete diagonal elements
    for i in 1:len
        prop[i, i] = 0.0
    end

    df_p = hcat(DataFrame(row_name=r_c_nm), DataFrame(prop, r_c_nm))
end

"""
    Workflow to export signature scores, p-values, and fold changes over combinations of cell labels and datasets.
"""
function write_csv_sig_scores_p_fc_wf(labels::Vector{String}, dataset_ids::Vector{String},
    filename_sig::String, res_dir::String, gene_set_id::T;
    min_features::Int64=2000, min_cells::Int64=2, sig_thresh::Float64=.6, tail::String="both",
    fc_thresh::Float64=1.15, adj_p_thresh::Float64=.01, n_c_thresh::Int64=100) where {T<:AbstractString}
    println("Cell labels: ", labels)
    println("Dataset IDs: ", dataset_ids)
    n_lb = length(labels)
    n_ds = length(dataset_ids)
    len = n_lb * n_ds
    lb_ds = [(i, j) for i in labels for j in dataset_ids]
    r_c_nm = [join([lb_ds[i][1], lb_ds[i][2]], "-") for i in eachindex(lb_ds)]

    arg_p = Vector{Vector{Vector{Tuple{String, String}}}}(undef, len-1)
    # lower triangular matrix excluding diagonal elements
    println("Begin calculating proportions of differential signature scores...")
    for (c, r) in enumerate([2:len;])
        arg_p[c] = [[lb_ds[i], lb_ds[j]] for i in r for j in 1:c]
    end

    @threads for i in eachindex(arg_p) # @threads
        @threads for j in eachindex(arg_p[i])
            label, ds_id = arg_p[i][j][1]
            label_2, ds_id_2 = arg_p[i][j][2]

            mat, genes, cell_ann = read_pbmc(label, ds_id; min_features=min_features, min_cells=min_cells)
            mat_2, genes_2, cell_ann_2 = read_pbmc(label_2, ds_id_2; min_features=min_features, min_cells=min_cells)

            obj, obj_2 = scSigScores_pair_from_pf(ds_id, mat, genes, cell_ann, ds_id_2, mat_2, genes_2, cell_ann_2, filename_sig)
            compute_sig_scores(obj; sig_thresh=sig_thresh)
            compute_sig_scores(obj_2; sig_thresh=sig_thresh)
            compute_pvalues(obj, obj_2, "Mann_Whitney_U"; tail=tail)
            calc_fold_changes(obj, obj_2; modify_obj=true)
            write_csv_sig_scores_p(obj, ds_id, label, obj_2, ds_id_2, label_2, res_dir, gene_set_id)
            write_csv_p_fc(obj, ds_id, label, ds_id_2, label_2, res_dir, gene_set_id)
        end
    end

    println("Workflow done.")
end
