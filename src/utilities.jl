#!/usr/bin/env julia
#= 
utility functions to support the single-cell signature scoring workflow
Version 1.4
=#

"""
    Initialise a scSigScores type.
"""
function init_scSigScores(filename_dt::String, filename_sig::String, dataset_id::String,
    path_exp_dataset::String, path_gene_names::String, path_cell_annot::String,
    filters::Dict{String, Vector{String}};
    kw...)
    (mat, genes, cell_ann) = get_subset_h5ad(filename_dt, path_exp_dataset, path_gene_names,
    path_cell_annot, filters; kw...)
    (mat_f, genes_f, cell_ann_f) = prefilter(mat, genes, cell_ann)
    signatures = read_gmt(filename_sig)
    obj = scSigScores(dataset_id, cell_ann_f, mat_f, genes_f, signatures)
    return obj
end

"""
    Create a pair of scSigScores objects with common genes from H5AD files.

Parameters
----------
filename_dt
    Name of the H5 file
filename_sig
    Name of the signature (signature pathway) file
dataset_id
    Dataset ID
fn_cell_annot
    Filename of the cell annotation CSV table (without header)
labels
    A vector of values of second column of cell annotation to select cells with, e.g. ["monocyte", "b_cell"]
path_exp_dataset
    Path to the gene expression dataset in the file, e.g. "matrix/data"
path_gene_names
    Path to the gene name dataset in the file, e.g. "matrix/features/name"
path_cell_annot
    Path to the cell type or cell annotation dataset in the file
filters
    A dictionary that specifies the filtering rule,
    mapping the path to the dataset and a vector of group values
convert_data_Int
    Whether to convert gene expression data to Int32 or not
drop_zeros
    Whether to drop zeros from the sparse matrix or not
remove_dup
    Whether to remove duplicated genes in the sparse matrix and genes
convert_data_Int
    Whether to convert data to integer values or not (true by default)
drop_zeros
    Whether to drop zeros in data (true by default)
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple of two scSigScores objects
"""
function scSigScores_pair_h5ad(filename_dt::String, filename_sig::String, dataset_id::String,
    path_exp_dataset::String, path_gene_names::String, path_cell_annot::String, filters::Dict{String, Vector{String}},
    filename_dt_2::String, filename_sig_2::String, dataset_id_2::String, # file_2
    path_exp_dataset_2::String, path_gene_names_2::String, path_cell_annot_2::String, filters_2::Dict{String, Vector{String}};
    kw...)
    # file_1
    kw_gsh = filter_kw(get_subset_h5ad; kw...)
    (mat, genes, cell_ann) = get_subset_h5ad(filename_dt, path_exp_dataset, path_gene_names,
    path_cell_annot, filters; kw_gsh...)
    kw_pf = filter_kw(prefilter; kw...)
    (mat, genes, cell_ann) = prefilter(mat, genes, cell_ann; kw_pf...)
    # file_2
    (mat_2, genes_2, cell_ann_2) = get_subset_h5ad(filename_dt_2, path_exp_dataset_2, path_gene_names_2,
    path_cell_annot_2, filters_2; kw_gsh...)
    (mat_2, genes_2, cell_ann_2) = prefilter(mat_2, genes_2, cell_ann_2; kw_pf...)

    (mat, genes, mat_2, genes_2) = match_genes(mat, genes, mat_2, genes_2)
    signatures = read_gmt(filename_sig)
    obj = scSigScores(dataset_id, cell_ann, mat, genes, signatures)
    obj_2 = scSigScores(dataset_id_2, cell_ann_2, mat_2, genes_2, signatures)
    return (obj, obj_2)
end

"""
    Create a pair of scSigScores objects with common genes from H5 files.

Parameters
----------
filename_dt
    Name of the H5 file
filename_sig
    Name of the signature (signature pathway) file
dataset_id
    Dataset ID
fn_cell_annot
    Filename of the cell annotation CSV table (without header)
labels
    A vector of values of second column of cell annotation to select cells with, e.g. ["monocyte", "b_cell"]
path_gene_names
    Path to the gene name dataset in the file, e.g. "matrix/features/name"
path_exp_dataset
    Path to the gene expression dataset in the file, e.g. "matrix/data"
path_cell_barcodes
    Path to the cell barcode dataset in the file, e.g. "matrix/barcodes"
convert_data_Int
    Whether to convert gene expression data to Int32 or not
drop_zeros
    Whether to drop zeros from the sparse matrix or not
remove_dup
    Whether to remove duplicated genes in the sparse matrix and genes
convert_data_Int
    Whether to convert data to integer values or not (true by default)
drop_zeros
    Whether to drop zeros in data (true by default)
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple of two scSigScores objects
"""
function scSigScores_pair_h5(filename_dt::String, filename_sig::String, dataset_id::String,
    fn_cell_annot::AbstractString, labels::Vector{<:AbstractString}, path_gene_names::AbstractString,
    path_exp_dataset::AbstractString, path_cell_barcodes::AbstractString,
    filename_dt_2::String, filename_sig_2::String, dataset_id_2::String, # file_2
    fn_cell_annot_2::AbstractString, labels_2::Vector{<:AbstractString}, path_gene_names_2::AbstractString,
    path_exp_dataset_2::AbstractString, path_cell_barcodes_2::AbstractString;
    kw...)
    # file_1
    kw_gsh = filter_kw(get_subset_h5; kw...)
    (mat, genes, cell_ann) = get_subset_h5(filename_dt, fn_cell_annot, labels, path_gene_names, 
    path_exp_dataset, path_cell_barcodes; kw_gsh...)
    kw_pf = filter_kw(prefilter; kw...)
    (mat, genes, cell_ann) = prefilter(mat, genes, cell_ann; kw_pf...)
    # file_2
    (mat_2, genes_2, cell_ann_2) = get_subset_h5(filename_dt_2, fn_cell_annot_2, labels_2, path_gene_names_2,
    path_exp_dataset_2, path_cell_barcodes_2; kw_gsh...)
    (mat_2, genes_2, cell_ann_2) = prefilter(mat_2, genes_2, cell_ann_2; kw_pf...)

    (mat, genes, mat_2, genes_2) = match_genes(mat, genes, mat_2, genes_2)
    signatures = read_gmt(filename_sig)
    obj = scSigScores(dataset_id, cell_ann, mat, genes, signatures)
    obj_2 = scSigScores(dataset_id_2, cell_ann_2, mat_2, genes_2, signatures)
    return (obj, obj_2)
end

"""
    Create scSigScores object pair from prefiltered data
"""
function scSigScores_pair_from_pf(dataset_id::String,
    data::SparseMatrixCSC{<:Real, <:Integer},
    genes::Vector{<:AbstractString},
    cell_annot::Vector{<:AbstractString},
    dataset_id_2::String,
    data_2::SparseMatrixCSC{<:Real, <:Integer}, # dataset 2
    genes_2::Vector{<:AbstractString},
    cell_annot_2::Vector{<:AbstractString},
    filename_sig::String;
    start_pos::Int64=3)
    if size(data)[2] < 100
        println("# cells less than 100 in $(dataset_id).")
    end
    if size(data_2)[2] < 100
        println("# cells less than 100 in $(dataset_id_2).")
    end
    (data, genes, data_2, genes_2) = match_genes(data, genes, data_2, genes_2)
    signatures = read_gmt(filename_sig; start_pos=start_pos)
    obj = scSigScores(dataset_id, cell_annot, data, genes, signatures)
    obj_2 = scSigScores(dataset_id_2, cell_annot_2, data_2, genes_2, signatures)
    return (obj, obj_2)
end

"""
    Create triple scSigScores objects from prefiltered data
"""
function scSigScores_tri_from_pf(dataset_id::String, # dataset 1
    data::SparseMatrixCSC{<:Real, <:Integer},
    genes::Vector{<:AbstractString},
    cell_annot::Vector{<:AbstractString},
    dataset_id_2::String, # dataset 2
    data_2::SparseMatrixCSC{<:Real, <:Integer},
    genes_2::Vector{<:AbstractString},
    cell_annot_2::Vector{<:AbstractString},
    dataset_id_3::String, # dataset 3
    data_3::SparseMatrixCSC{<:Real, <:Integer},
    genes_3::Vector{<:AbstractString},
    cell_annot_3::Vector{<:AbstractString},
    filename_sig::String;
    start_pos::Int64=3)
    if size(data)[2] < 100
        println("# cells less than 100 in $(dataset_id).")
    end
    if size(data_2)[2] < 100
        println("# cells less than 100 in $(dataset_id_2).")
    end
    if size(data_3)[2] < 100
        println("# cells less than 100 in $(dataset_id_2).")
    end
    (data, genes, data_2, genes_2) = match_genes(data, genes, data_2, genes_2)
    (data, genes, data_3, genes_3) = match_genes(data, genes, data_3, genes_3)
    signatures = read_gmt(filename_sig; start_pos=start_pos)
    obj = scSigScores(dataset_id, cell_annot, data, genes, signatures)
    obj_2 = scSigScores(dataset_id_2, cell_annot_2, data_2, genes_2, signatures)
    obj_3 = scSigScores(dataset_id_3, cell_annot_3, data_3, genes_3, signatures)
    return (obj, obj_2, obj_3)
end

"""
    Pre-filter genes and cells with zero gene expression along each dimension respectively.

Parameters
----------
data
    An m-by-n CSC sparse matrix of gene expression values (m genes, n cells)
genes
    A vector of genes
convert_data_Int
    Whether to convert data to integer values or not (true by default)
drop_zeros
    Whether to drop zeros in data (true by default)
min_features
    Minimum number of features (genes) to filter cells
min_cells
    Minimum number of cells to filter features

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
and a vector of gene symbols
"""
function prefilter(data::SparseMatrixCSC{<:Real, <:Integer}, genes::Vector{<:AbstractString},
    cell_annot::Vector{<:AbstractString};
    convert_data_Int::Bool=false, drop_zeros::Bool=false, min_cells::Integer=0, min_features::Integer=0)
    if drop_zeros
        dropzeros!(data) # remove stored zeros
    end

    if min_features >= 1
        data_f, cell_ann_f = filter_cells_n_f(data, cell_annot; min_features=min_features)
    else
        data_f = data
        cell_ann_f = cell_annot
    end

    if min_cells >= 1
        data_f, genes_f = filter_genes(data_f, genes; min_cells=min_cells)
    else
        genes_f = genes
    end

    if convert_data_Int
        try
            data_f = SparseMatrixCSC{Int32, Int32}(data_f)
        catch
            print("
            Cannot convert a matrix of floating-point numbers to a sparse matrix of integers. Skip type conversion.
            ")
            data_f = SparseMatrixCSC{eltype(data_f), Int32}(data_f)
        end
    else
        data_f = SparseMatrixCSC{eltype(data_f), Int32}(data_f)        
    end

    println("Number of filtered genes: $(length(genes_f)).")
    println("Number of filtered cells: $(size(data_f)[2]).")
    return (data_f, genes_f, cell_ann_f)
end

"""
    prefilter function for data of type Matrix.
"""
prefilter(data::Matrix{<:Real}, genes::Vector{<:AbstractString}, cell_annot::Vector{<:AbstractString}) = prefilter(sparse(data), genes, cell_annot; convert_data_Int=true)

"""
    prefilter function for data of type DataFrame.
"""
prefilter(data::DataFrame, cell_annot::Vector{<:AbstractString}) = prefilter(sparse(Matrix(data[:, 2:end])), data[:, 1], cell_annot; convert_data_Int=true)

"""
    Keep only common genes in two scSigScores objects.
"""
function match_genes(data::SparseMatrixCSC{<:Real, <:Integer}, genes::Vector{<:AbstractString},
    data_2::SparseMatrixCSC{<:Real, <:Integer}, genes_2::Vector{<:AbstractString})
    c_genes = intersect(genes, genes_2)
    keep = findall(x -> in(x, c_genes), genes)
    keep_2 = findall(x -> in(x, c_genes), genes_2)
    data_f = data[keep, :]
    genes_f = genes[keep]
    data_2_f = data_2[keep_2, :]
    genes_2_f = genes_2[keep_2]
    println("Number of genes in common: $(length(c_genes)).")
    return (data_f, genes_f, data_2_f, genes_2_f)
end

"""
    Filter cells by number of features.
"""
function filter_cells_n_f(data::SparseMatrixCSC{<:Real, <:Integer}, cell_annot::Vector{<:AbstractString};
    min_features::Integer=2000)
    n_f = sum(data .> 0, dims=1)
    f = vec(n_f .>= min_features)
    data_f = data[:, f]
    cell_ann_f = cell_annot[f]
    println("Number of filtered cells: $(size(data_f)[2]).")
    return (data_f, cell_ann_f)
end

"""
    Filter cells by library size.
"""
function filter_cells_lib_s(data::SparseMatrixCSC{<:Real, <:Integer}, cell_annot::Vector{<:AbstractString};
    min_lib_size::Integer=5000)
    n_c = sum(data, dims=1)
    f = vec(n_c .>= min_lib_size)
    data_f = data[:, f]
    cell_ann_f = cell_annot[f]
    println("Number of filtered cells: $(size(data_f)[2]).")
    return (data_f, cell_ann_f)
end

"""
    Filter genes by number of cells expressing the genes.
"""
function filter_genes(data::SparseMatrixCSC{<:Real, <:Integer}, genes::Vector{<:AbstractString};
    min_cells::Integer=3)
    f = vec(sum(data .> 0, dims=2) .>= min_cells)
    data_f = data[f, :]
    genes_f::Vector{String} = genes[f]
    println("Number of filtered genes: $(length(genes_f)).")
    return (data_f, genes_f)
end
