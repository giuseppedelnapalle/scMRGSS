#!/usr/bin/env julia
#= 
read and write to files
nikola
Dec 2022
version 1.3
=#

# using HDF5, SparseArrays, CSV, DataFrames

"""
    Read a gmt file and output a dictionary of gene sets.
"""
function read_gmt(filename::AbstractString; start_pos::Int64=3, delim::String="\t") # line order not preserved
    sig = Dict{String, Vector{String}}()
    open(filename, "r") do f
        for ln in eachline(f)
            s = split(ln, delim)
            sig[string(s[1])] = s[start_pos:end]
        end
    end
    return sig
end

"""
    Access an H5AD file and output a subset of gene expression data.

Gene expression data in the H5AD file are saved in the Compressed Sparse Row (CSR)
sparse matrix format specifying three NumPy arrays, i.e. data, indices, and indptr with 0-based indexing.
A transpose of the matrix subset in the format of Compressed Sparse Column (CSR) will be returned.

Parameters
----------
filename
    Name of the H5AD file
path_exp_dataset
    Path to the gene expression dataset in the file, e.g. "raw/X/data"
path_gene_names
    Path to the gene name dataset in the file
path_cell_annot
    Path to the cell type or cell annotation dataset in the file
filters
    A dictionary that specifies the filtering rule, mapping the path to the dataset and a vector of group values
path_obs_categories
    Path to the obs (row) categories group that stores coding information for each category
path_var_categories
    Path to the var (column) categories group that stores coding information for each category
full_path_gene_cat
    Path to gene name categories (e.g. var/gene_symbol/categories). Overrides path_var_categories
path_obs_categories
    Path to the obs (row) categories group that stores coding information for each category
full_path_cell_cat
    Path to cell annotation categories (e.g. obs/cell_ontology_class/categories). Overrides path_obs_categories
convert_data_Int
    Whether to convert gene expression data to Int32 or not
drop_zeros
    Whether to drop zeros from the sparse matrix or not
remove_dup
    Whether to remove duplicated genes in the sparse matrix and genes

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell labels for cells
"""
function get_subset_h5ad(filename::AbstractString, path_exp_dataset::AbstractString,
    path_gene_names::AbstractString, path_cell_annot::AbstractString, filters::Dict{Tk, Vector{Tv}};
    path_var_categories::AbstractString="var/__categories", full_path_gene_cat::AbstractString="",
    path_obs_categories::AbstractString="obs/__categories", full_path_cell_cat::AbstractString="",
    convert_data_Int::Bool=true, drop_zeros::Bool=true,
    remove_dup::Bool=true) where {Tk<:AbstractString, Tv<:AbstractString}
    genes = gene_names_h5ad(filename, path_gene_names; path_var_categories=path_var_categories, full_path_gene_cat=full_path_gene_cat)
    (mat, cell_ann) = main_subset_h5ad(filename, genes, path_exp_dataset, path_cell_annot, filters; 
    path_obs_categories=path_obs_categories, full_path_cell_cat=full_path_cell_cat, convert_data_Int=convert_data_Int, drop_zeros=drop_zeros)
    if remove_dup
        (mat, genes) = remove_dup_genes(mat, genes)
    end
    return (mat, genes, cell_ann)
end

"""
    Access an H5AD file and output a subset of gene expression data.
"""
function main_subset_h5ad(filename::AbstractString, genes::Vector{<:AbstractString}, 
    path_exp_dataset::AbstractString, path_cell_annot::AbstractString, filters::Dict{Tk, Vector{Tv}}; 
    path_obs_categories::AbstractString="obs/__categories", full_path_cell_cat::AbstractString="",
    convert_data_Int::Bool=true, drop_zeros::Bool=true) where {Tk<:AbstractString, Tv<:AbstractString}
    path_indices = replace(path_exp_dataset, "data" => "indices")
    path_indptr = replace(path_exp_dataset, "data" => "indptr")
    ind_ipt = indices_indptr(filename, filters; path_obs_categories=path_obs_categories, full_path_cell_cat=full_path_cell_cat)
    cell_ann = get_cell_annot(filename, path_cell_annot, ind_ipt; path_obs_categories=path_obs_categories, full_path_cell_cat=full_path_cell_cat)
    aux_i = ind_ipt .+ 1
    indptr = h5read(filename, path_indptr)[ind_ipt]
    indptr_aux = h5read(filename, path_indptr)[aux_i]
    r_indices::Vector{Int32} = row_indices(indptr, indptr_aux) # J
    ipt_filter = indptr_filter(indptr, indptr_aux)
    indices::Vector{Int32} = h5read(filename, path_indices)[ipt_filter] .+ 1 # I
    if convert_data_Int
        try
            data = Vector{Int32}(h5read(filename, path_exp_dataset)[ipt_filter]) # V
        catch
            data = h5read(filename, path_exp_dataset)[ipt_filter]
            print("
            Cannot convert a vector of floating-point numbers to a vector of integers. Skip type conversion.
            ")
        end
    else
        data = h5read(filename, path_exp_dataset)[ipt_filter]
    end
    m = length(genes)
    n = maximum(r_indices) # number of cells
    mat = sparse(indices, r_indices, data, m, n) # I, J, V
    if drop_zeros
        dropzeros!(mat) # remove stored zeros
    end
    return (mat, cell_ann) # SparseMatrixCSC{eltype(data), Int32}
end

"""
    Read an H5 file.

Parameters
----------
filename
    Name of the Loom file
path_gene_names
    Path to the gene name dataset in the file, e.g. "matrix/features/name"
path_exp_dataset
    Path to the gene expression dataset in the file, e.g. "matrix/data"
convert_data_Int
    Whether to convert gene expression data to Int32 or not
drop_zeros
    Whether to drop zeros from the sparse matrix or not
remove_dup
    Whether to remove duplicated genes in the sparse matrix and genes

Output
----------
    A tuple consisting of an m-by-n matrix of gene expression values (m genes, n cells)
    and a vector of gene symbols
"""
function read_h5(filename::AbstractString, path_gene_names::AbstractString="matrix/features/name",
    path_exp_dataset::AbstractString="matrix/data";
    convert_data_Int::Bool=true, drop_zeros::Bool=true, remove_dup::Bool=true)
    path_indptr = replace(path_exp_dataset, "data" => "indptr")
    path_indices = replace(path_exp_dataset, "data" => "indices")
    indptr = h5read(filename, path_indptr)
    r_indices::Vector{Int32} = get_row_ind(indptr) # J
    indices::Vector{Int32} = h5read(filename, path_indices) .+ 1 # I
    if convert_data_Int
        try
            data = Vector{Int32}(h5read(filename, path_exp_dataset)) # V
        catch
            data = h5read(filename, path_exp_dataset)
            print("
            Cannot convert a vector of floating-point numbers to a vector of integers. Skip type conversion.
            ")
        end
    else
        data = h5read(filename, path_exp_dataset)
    end
    genes = h5read(filename, path_gene_names)
    m = length(genes)
    n = maximum(r_indices) # number of cells
    mat = sparse(indices, r_indices, data, m, n) # I, J, V
    if drop_zeros
        dropzeros!(mat) # remove stored zeros
    end
    if remove_dup
        (mat, genes) = remove_dup_genes(mat, genes)
    end
    return (mat, genes)
end

"""
Create row indices from indptr.
"""
function get_row_ind(indptr::Vector{<:Integer}) # 1-based indexing
    r_indices = Vector{Vector{Int32}}(undef, length(indptr)-1)
    for i in eachindex(r_indices)
        r_indices[i] = repeat([i], indptr[i+1]-indptr[i])
    end
    r_ind_vc = vcat(r_indices...)
    return r_ind_vc
end

"""
    Access an H5 file and output a subset of gene expression data.

Gene expression data in the H5AD file are saved in the Compressed Sparse Row (CSR)
sparse matrix format specifying three NumPy arrays, i.e. data, indices, and indptr with 0-based indexing.
A transpose of the matrix subset in the format of Compressed Sparse Column (CSR) will be returned.

Parameters
----------
filename
    Name of the H5 file
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

Output
----------
A tuple consisting of an m-by-n CSC sparse matrix of gene expression values (m genes, n cells),
a vector of gene symbols, and a vector of cell labels for the cells
"""
function get_subset_h5(filename::AbstractString, fn_cell_annot::AbstractString, labels::Vector{<:AbstractString},
    path_gene_names::AbstractString="matrix/features/name",
    path_exp_dataset::AbstractString="matrix/data", path_cell_barcodes::AbstractString="matrix/barcodes";
    convert_data_Int::Bool=true, drop_zeros::Bool=true, remove_dup::Bool=true)
    path_indptr = replace(path_exp_dataset, "data" => "indptr")
    path_indices = replace(path_exp_dataset, "data" => "indices")

    c_ann_df = CSV.read(fn_cell_annot, DataFrame; header=0)
    df_sub = subset(c_ann_df, :Column2 => ByRow(x -> x in labels))
    cells_f::Vector{String} = df_sub[:, 1] # filtered cells
    cell_ann::Vector{String} = df_sub[:, 2]
    
    c_barcodes = h5read(filename, path_cell_barcodes)
    ind_ipt = findall(x -> x in cells_f, c_barcodes) # row indices of select cells
    aux_i = ind_ipt .+ 1
    indptr = h5read(filename, path_indptr)[ind_ipt]
    indptr_aux = h5read(filename, path_indptr)[aux_i]
    r_indices::Vector{Int32} = row_indices(indptr, indptr_aux) # J
    ipt_filter = indptr_filter(indptr, indptr_aux)
    indices::Vector{Int32} = h5read(filename, path_indices)[ipt_filter] .+ 1 # I
    genes = h5read(filename, path_gene_names)
    if convert_data_Int
        try
            data = Vector{Int32}(h5read(filename, path_exp_dataset)[ipt_filter]) # V
        catch
            data = h5read(filename, path_exp_dataset)[ipt_filter]
            print("
            Cannot convert a vector of floating-point numbers to a vector of integers. Skip type conversion.
            ")
        end
    else
        data = h5read(filename, path_exp_dataset)[ipt_filter]
    end
    m = length(genes)
    n = maximum(r_indices) # number of cells
    mat = sparse(indices, r_indices, data, m, n) # I, J, V
    if drop_zeros
        dropzeros!(mat) # remove stored zeros
    end
    if remove_dup
        (mat, genes) = remove_dup_genes(mat, genes)
    end
    return (mat, genes, cell_ann)
end

"""
    Access a Loom file and output a subset of gene expression data.

Parameters
----------
filename
    Name of the Loom file
path_gene_names
    Path to the gene name dataset in the file
path_cell_names
    Path to the cell name dataset in the file
filtered_cells
    A vector of cell names
path_exp_dataset
    Path to the gene expression dataset in the file, e.g. "matrix"
output_csc_matrix
    Whether to output a CSC sparse matrix or not
convert_data_Int
    Whether to convert gene expression data to Int32 or not
drop_zeros
    Whether to drop zeros from the sparse matrix or not
remove_dup
    Whether to remove duplicated genes in the sparse matrix and genes

Output
----------
A tuple consisting of an m-by-n matrix of gene expression values (m genes, n cells)
and a vector of gene symbols
"""
function get_subset_loom(filename::AbstractString, path_gene_names::AbstractString, path_cell_names::AbstractString,
    filtered_cells::Vector{<:AbstractString}, path_exp_dataset::AbstractString="matrix";
    output_csc_matrix::Bool=true, convert_data_Int::Bool=true, drop_zeros::Bool=true, remove_dup::Bool=false)
    genes = h5read(filename, path_gene_names)
    cells = h5read(filename, path_cell_names)
    col_filter = indexin(filtered_cells, cells)
    if in(nothing, col_filter)
        println("Some filtered cells not found in the matrix.")
    end
    filter!(x -> x !== nothing, col_filter)
    if length(col_filter) != length(unique(col_filter))
        println("Duplicated cells found and removed.")
        unique!(col_filter)
    end
    if length(col_filter) == 0
        throw(error("No filtered cell found in the matrix. Please perform another filtering."))
    end
    if convert_data_Int
        try
            data = Matrix{Int32}(h5read(filename, path_exp_dataset)[:, col_filter])
        catch
            data = h5read(filename, path_exp_dataset)[:, col_filter]
            print("
            Cannot convert a matrix of floating-point numbers to a matrix of integers. Skip type conversion.
            ")
        end
    else
        data = h5read(filename, path_exp_dataset)[:, col_filter]
    end
    if output_csc_matrix
        mat = SparseMatrixCSC{eltype(data), Int32}(data)
        if drop_zeros
            dropzeros!(mat) # remove stored zeros
        end
    else
        mat = data
    end
    if remove_dup
        (mat, genes) = remove_dup_genes(mat, genes)
    end
    return (mat, genes)
end

"""
    Create a dictionary for the filters argument.

Parameters
----------
args
    multiple Pairs or Tuples specifying the filtering rule
    each Pair or Tuple stores the path to the dataset and a vector of values, 
    e.g. "obs/cell_type/codes" => ["t cell"], or ("obs/__categories/cell_type", ["t cell"]).
"""
function create_filters(args::Union{
    Pair{Tk, Vector{Tv}}, Tuple{Tk, Vector{Tv}}
    }...) where {Tk<:AbstractString, Tv<:AbstractString}
    filters = Dict{String, Vector{String}}()
    for arg in args
        if !(typeof(arg[2]) <: Vector{<:AbstractString})
            throw(ArgumentError("The second field of the argument should be a a vector of strings."))
        end
        filters[String(arg[1])] = Vector{String}(arg[2])
    end
    return filters
end

"""
    Write the result to a CSV file.

Parameters
----------
object
    A scSigScores object
path_file
    Path to the file
key
    A key of the results field to specify the DataFrame of interest.
    Choose from "sig_scores", "two_sample_test", "fold_change", and "is_diff".

Output
A CSV file
"""
function write_csv(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    path_to_file::AbstractString, key::AbstractString; kw...)
    k = unique(keys(object.results[key]))[1]
    df = object.results[key][k]
    CSV.write(path_to_file, df; kw...)
end

"""
    Access a H5AD file and output a vector of gene names.
"""
function gene_names_h5ad(filename::AbstractString, path_gene_names::AbstractString;
    path_var_categories::AbstractString="var/__categories", full_path_gene_cat::AbstractString="")
    gene_names = h5read(filename, path_gene_names)
    if eltype(gene_names) <: AbstractString # situation in which gene names saved in the object
        genes = gene_names
    else # situation in which code numbers corresponding to gene names saved in the object
        if !isempty(full_path_gene_cat) # full_path_gene_cat overrides path_var_categories
            coding = h5read(filename, full_path_gene_cat)
        else
            p = joinpath(path_var_categories, split(path_gene_names, "/")[end]) # compatible with old version of TS dataset
            coding = h5read(filename, p)
        end

        codes = Dict{Int, String}(zip([0:(length(coding) - 1);], coding))
        genes = Vector{String}(undef, length(gene_names))

        for i in eachindex(gene_names)
            genes[i] = codes[gene_names[i]]
        end
    end

    return genes
end

"""
    Obtain indices of indptr corresponding to selected cells.

Parameters
----------
full_path_cell_cat
    Full path to the cell category (coding) dataset (values of type string) in the file.
"""
function indices_indptr(filename::AbstractString, filters::Dict{Tk, Vector{Tv}};
    path_obs_categories::AbstractString="obs/__categories", full_path_cell_cat::AbstractString="") where {Tk<:AbstractString, Tv<:AbstractString}
    i_grp = Vector{Vector{Int32}}()
    for (key, value) in filters
        cat = h5read(filename, key) # Int vector encoding cell categories

        if !isempty(full_path_cell_cat)
            if occursin("codes", key)
                coding = h5read(filename, replace(key, "codes" => "categories"))
            else
                coding = h5read(filename, joinpath(path_obs_categories, split(key, "/")[end]))
            end
        else
            coding = h5read(filename, full_path_cell_cat)
        end

        coding = [rstrip(str) for str in coding] # remove trailing space

        i_cat = Int32[]
        for v in value
            c = findfirst(x -> x == v, coding) - 1 # integer code for the value in the categories
            if c===nothing
                throw(ArgumentError("Current filter: $(key) => $(v). No cell matching the filter found."))
            end
            r = findall(x -> x == c, cat) # row indices for selected cells
            append!(i_cat, r)
        end
        push!(i_grp, i_cat)
    end
    ind_ipt = intersect(i_grp...) # 1-based indexing
    return ind_ipt # Vector{Int32}
end

"""
    Find cell type labels of filtered cells from a H5AD file.
"""
function get_cell_annot(filename::AbstractString, path_cell_annot::AbstractString, ind_ipt::Vector{Int32};
    path_obs_categories::AbstractString="obs/__categories", full_path_cell_cat::AbstractString="")
    c_t = h5read(filename, path_cell_annot) # 0-based indexing

    if !isempty(full_path_cell_cat)
        coding = h5read(filename, full_path_cell_cat)
    else
        coding = h5read(filename, joinpath(path_obs_categories, split(path_cell_annot, "/")[end]))
    end

    c_t_f = c_t[ind_ipt] # int code of filtered cells
    cell_ann = Vector{String}(undef, length(c_t_f))
    for i in eachindex(c_t_f)
        cell_ann[i] = coding[c_t_f[i]+1]
    end
    return cell_ann
end

"""
    Obtain row indices of sparse matrix.
"""
function row_indices(indptr::Vector{T}, indptr_aux::Vector{T}) where {T<:Integer}
    r_indices = Vector{Vector{Int32}}(undef, length(indptr))
    for i in eachindex(indptr)
        r_indices[i] = repeat([i], indptr_aux[i]-indptr[i])
    end
    r_ind_vc = vcat(r_indices...)
    return r_ind_vc
end

"""
    Remove duplicated genes by keeping genes of first appearance.
"""
function remove_dup_genes(data::SparseMatrixCSC{<:Real, <:Integer}, genes::Vector{<:AbstractString})
    u_genes = unique(genes)
    filter = zeros(Int32, length(u_genes))
    for i in eachindex(filter)
        filter[i] = findfirst(x -> x == u_genes[i], genes)
    end
    data_f = data[filter, :]
    genes_f = genes[filter]
    return (data_f, genes_f)
end

"""
    Obtain indices of elements from data and indices for all selected cells.
"""
function indptr_filter(indptr::Vector{T}, indptr_aux::Vector{T}) where {T<:Integer}
    ipt_filter = Vector{Vector{Int32}}(undef, length(indptr))
    for i in eachindex(indptr)
        ipt_filter[i] = [(indptr[i]+1):indptr_aux[i];]
    end
    ipt_f_vc = vcat(ipt_filter...)
    return ipt_f_vc
end

"""
    Select cells from meta data.
"""
function filter_cells(filename::AbstractString, filters::Dict{Tk, Vector{Tv}},
    column::Union{Integer, AbstractString}=1, delim::String=",") where {Tk<:AbstractString, Tv<:AbstractString}
    meta = CSV.read(filename, DataFrame; delim=delim)
    i_all = Vector{Vector{Int32}}()
    for (key, value) in filters
        cat = meta[!, Symbol(key)]
        i_cat = Int32[]
        for v in value
            rows = findall(x -> x == v, cat)
            if length(rows)==0
                throw(ArgumentError("Current filter: $(key) => $(v). No cell matching the filter found."))
            end
            append!(i_cat, rows)
        end
        push!(i_all, i_cat)
    end
    filter = intersect(i_all...)
    if isa(column, Integer) # row name column set as first column by default
        c = column
    else
        c = Symbol(column)
    end
    cells::Vector{String} = meta[filter, c]
    return cells
end

"""
    Write data frames of signature scores for two objects to separate CSV file.
"""
function write_csv_sig_scores_p(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    ds_id::String, label::String,
    object_ref::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    ds_id_ref::String, label_ref::String,
    res_dir::String, gene_set_id::T; kw...) where {T<:AbstractString}
    pfx = "sig_scores_$(gene_set_id)_pair_$(ds_id)_$(label)_ref_$(ds_id_ref)_$(label_ref)"
    fn = joinpath(res_dir, "$(pfx)-$(ds_id)_$(label).csv")
    write_csv(object, fn, "sig_scores"; kw...)

    fn_2 = joinpath(res_dir, "$(pfx)-$(ds_id_ref)_$(label_ref).csv")
    write_csv(object_ref, fn_2, "sig_scores"; kw...)
end

"""
    Write p-values, adj. p-values, fold changes for the target object to a CSV file.
"""
function write_csv_p_fc(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    ds_id::String, label::String, ds_id_2::String, label_2::String, res_dir::String, gene_set_id::T;
    kw...) where {T<:AbstractString}
    fn = joinpath(res_dir, "p_fc_$(gene_set_id)_pair_$(ds_id)_$(label)_ref_$(ds_id_2)_$(label_2).csv")
    p = access_result(object, "two_sample_test")
    fc = access_result(object, "fold_change")
    df = innerjoin(p, fc; on=:signature, makeunique=true)
    CSV.write(fn, df; kw...)
end
