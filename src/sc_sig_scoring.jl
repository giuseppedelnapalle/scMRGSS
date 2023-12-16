#!/usr/bin/env julia
#= 
Single-cell signature scoring method based on gene ranks using a reference cell type
for significance evaluation
Input
    1. A single-cell RNA-seq data matrix (m-by-n) or data frame (m-by-(n+1)) for cells of interest, 
    where rows represent genes and columns represent cells
    2. A single-cell RNA-seq data matrix or data frame for the reference cell type
    3. A vector of genes of the signature (a priori gene set) for scoring
Output
    1. A data frame of signature scores for cells of interest
    2. A data frame of p-values and adjusted p-values for the associated signature scores
nikola
Nov 2023
Version 1.6.7
=#

# using DataFrames, Statistics, SparseArrays, HypothesisTests, MultipleTesting, Loess

"""
scSigScores{Ts<:AbstractString, Tk<:AbstractString, Tv<:AbstractString}

A composite type for storing the input and output data to implement
the rank-based signature scoring method.
"""
struct scSigScores{Ts<:AbstractString, Tk<:AbstractString, Tv<:AbstractString}
dataset_id::Ts
cell_annot::Vector{<:AbstractString}
data::SparseMatrixCSC{<:Real, <:Integer}
genes::Vector{<:AbstractString}
signatures::Dict{Tk, Vector{Tv}}
results::Dict{String, Dict{String, DataFrame}}
end

"""
    Set the results field to an empty dictionary by default.
"""
scSigScores(dataset_id::Ts, cell_annot::Vector{<:AbstractString},
data::SparseMatrixCSC{<:Real, <:Integer},
genes::Vector{<:AbstractString},
signatures::Dict{Tk, Vector{Tv}}) where {Ts<:AbstractString, Tk<:AbstractString, Tv<:AbstractString} =
scSigScores(dataset_id, cell_annot, data, genes, signatures, Dict{String, Dict{String, DataFrame}}())

"""
    Constructor for data of type Matrix.
"""
scSigScores(dataset_id::Ts, cell_annot::Vector{<:AbstractString},
data::Matrix, genes::Vector{<:AbstractString},
signatures::Dict{Tk, Vector{Tv}}) where {Ts<:AbstractString, Tk<:AbstractString, Tv<:AbstractString} =
scSigScores(dataset_id, cell_annot, sparse(data), genes, signatures)

"""
    Constructor for data of type DataFrame.
"""
scSigScores(dataset_id::Ts, cell_annot::Vector{<:AbstractString},
data::DataFrame, signatures::Dict{Tk, Vector{Tv}}) where {Ts<:AbstractString, Tk<:AbstractString, Tv<:AbstractString} =
scSigScores(dataset_id, cell_annot, df_to_csc_sparse(data), genes_from_df(data), signatures)

"""
    Convert a DataFrame to a CSC sparse matrix.
"""
function df_to_csc_sparse(df::DataFrame; convert_data_Int::Bool=true)
    if convert_data_Int
        try
            mat = SparseMatrixCSC{Int32, Int32}(Matrix(df[:, 2:end]))
        catch
            print("
            Cannot convert a data frame of floating-point numbers to a sparse matrix of integers. Skip type conversion.
            ")
            mat = Matrix(df[:, 2:end])
            mat = SparseMatrixCSC{eltype(mat), Int32}(mat)
        end
    else
        mat = Matrix(df[:, 2:end])
        mat = SparseMatrixCSC{eltype(mat), Int32}(mat)
    end
    return mat
end

"""
    Obtain genes from a DataFrame of gene expression values.
"""
genes_from_df(df::DataFrame)::Vector{String} = df[:, 1]

"""
    Create scSigScores object with signatures field derived from a gmt file.
"""
function scSigScores_sig_gmt(dataset_id::String,
    data::SparseMatrixCSC{<:Real, <:Integer},
    genes::Vector{<:AbstractString},
    cell_annot::Vector{<:AbstractString},
    filename_sig::String)
    if size(data)[2] < 100
        println("# cells less than 100 in $(dataset_id).")
    end
    signatures = read_gmt(filename_sig)
    obj = scSigScores(dataset_id, cell_annot, data, genes, signatures)
end

"""
    Compute signature scores and adjust for a covariate for a pair of scSigScores objects.
"""
function compute_adj_sig_scores_p(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    object_2::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString}, covariate::String="n_counts";
    kw...)
    compute_sig_scores(object; kw...)
    compute_sig_scores(object_2; kw...)
    adjust_mean_ranks(object, covariate)
    adjust_mean_ranks(object_2, covariate)
end

"""
    Compute signature scores based on gene ranks across cells for each signature (gene set).

Parameters
----------
object
    A scSigScores object
sig_thresh
    Signature threshold. Calculation of signature scores will be aborted if the proportion
    of signature genes in the data is below the threshold.

Output
----------
An m-by-n DataFrame of signature scores across cells for each signature (m gene sets, n cells)
"""
function compute_sig_scores(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString};
    kw...)
    n_c = size(object.data)[2]
    n_s = length(object.signatures)
    println("Number of genes: $(size(object.data)[1]).")
    println("Number of cells: $(n_c).")
    println("Number of signatures: $(n_s).")
    if n_c==0
        error("No cell fouond. Please adjust the threshod of filtering.")
    end
    ranks = mat_to_ranks(object.data, object.genes)
    keys = Vector{String}(undef, n_s)
    scores = Matrix{Union{Float64, Missing}}(undef, n_s, n_c)
    for (i, item) in enumerate(object.signatures)
        keys[i] = first(item)
        scores[i, :] = rank_sig_scores(ranks, object.genes, last(item); kw...)
    end
    df = hcat(DataFrame(signature=keys), DataFrame(scores, :auto))
    object.results["sig_scores"] = Dict{String, DataFrame}(object.dataset_id => df)
end

"""
    Convert a sparse matrix to a vector of pairs for ranks.
"""
function mat_to_ranks(data::SparseMatrixCSC{<:Real, <:Integer}, genes::Vector{<:AbstractString})
    n_c = size(data)[2] # num of cells
    ranks = Vector{Pair{Vector{String}, Vector{Int32}}}(undef, n_c)
    for i in 1:n_c
        v = data[:, i] # SparseVector
        i_exp_g = v .!= 0
        v_s = v[i_exp_g] # only expressed genes are considered in mean rank calculation
        p = sortperm(v_s) # indices of sorted array in ascending order
        rk_v = zeros(Int32, sum(i_exp_g))
        rk_v[p] = 1:length(v_s)
        ranks[i] = genes[i_exp_g] => rk_v
    end
    return ranks # ranks[i] > 0 for all genes
end

"""
    Compute rank based signature scores across cells for a signature.
"""
function rank_sig_scores(ranks::Vector{Pair{Vector{String}, Vector{Int32}}}, genes::Vector{<:AbstractString},
    signature::Vector{<:AbstractString}; sig_thresh::T=.6) where {T<:AbstractFloat}
    test = test_signature(genes, signature; sig_thresh=sig_thresh)
    if !test
        print("
        Proportion of signature genes in the data is below the threshold $(sig_thresh). Procedure aborted.
        ")
        return repeat([missing], length(ranks)) # Vector{Missing}
    end
    mean_ranks = zeros(Float64, length(ranks))
    for i in eachindex(ranks)
        rk_sig = last(ranks[i])[findall(x -> x in signature, first(ranks[i]))]
        mean_rk = length(rk_sig) > 0 ? mean(rk_sig) : 1e-6
        mean_ranks[i] = mean_rk / length(last(ranks[i])) # normalise by number of expressed genes
    end
    return mean_ranks
end

""""
    Compute p-values for the signature scores over signatures.

The p-value is defined as the probability of obtaining more extreme OR by chance, 
under the assumption of the null hypothesis, i.e. OR of cells of intest = OR of reference cell type.

Parameters
----------
object
    A scSigScores object
ref_object
    A scSigScores object storing the reference data
proc
    Procedure to perform between-group test on signature scores.
    Choose from "Mann_Whitney_U" (default) and "approx_permutation".
tail
    Type of the hypothesis test. Choose from "right" , "left" or "both (right)".
centroid
    Measure of central tendency. Either median or mean. Only applicable to the "approx_permutation" procedure.
n
    Sampling size of permutations. Only applicable to the "approx_permutation" procedure.
    Defaults to 10_000.
adjust_p
    Procedure to perform p-value adjustment for multiple comparisons.
    Choose from "Bonferroni" (default), "BenjaminiHochberg", "BenjaminiYekutieli", "Hochberg","Holm",
    "Hommel", "Sidak", or "forward_stop".

Output
----------
A DataFrame of p-values and adjusted p-values with each row corresponding to a signature
"""
function compute_pvalues(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    ref_object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString}, proc::String="Mann_Whitney_U";
    tail::String="both", centroid::String="median", n::Int=10000, adjust_p::String="Bonferroni", filter_scores::Bool=true)
    if !in(tail, ["right", "left", "both"])
        throw(ArgumentError("""
        Incorrect argument for tail. Please choose from "right", "left", or "both".
            """))
    end
    t = ifelse(tail=="both", "two", tail)
    println("Perform $(t)-sided $(proc) procedure.")

    scores = Matrix{Union{Float64, Missing}}(access_result(object, "sig_scores")[:,2:end])
    ref_scores = Matrix{Union{Float64, Missing}}(access_result(ref_object, "sig_scores")[:,2:end])
    p_values = two_sample_test(scores, ref_scores, proc; tail=tail, centroid=centroid, n=n, filter_scores=filter_scores)
    adj_p = adjust_p_values(p_values; adjust_p=adjust_p)

    if any(ismissing, p_values)
        if !all(ismissing, p_values)
            adj_p = insert_missing(adj_p, p_values)
        end
    end

    sig_names = DataFrame(signature=object.results["sig_scores"][object.dataset_id][:, 1])
    if !isempty(p_values)
        df = hcat(sig_names, DataFrame(p_value=p_values), DataFrame(adj_p=adj_p))
    else
        v_missing = repeat([missing], size(sig_names)[1])
        df = hcat(sig_names, DataFrame(p_value=v_missing), DataFrame(adj_p=v_missing))
    end

    if "two_sample_test" in keys(object.results)
        object.results["two_sample_test"][ref_object.dataset_id] = df
    else
        object.results["two_sample_test"] = Dict{String, DataFrame}(ref_object.dataset_id => df)
    end
end

"""
    Perform a two-sample test on signature scores from two groups.
"""
function two_sample_test(scores::Matrix{Union{Float64, Missing}}, ref_scores::Matrix{Union{Float64, Missing}},
    proc::String="Mann_Whitney_U";
    tail::String="both", centroid::String="median", n::Int=10000, filter_scores::Bool=true)
    if proc=="Mann_Whitney_U"
        p_values = proc_dispatch(scores, ref_scores, mwu_test; tail=tail, filter_scores=filter_scores)
    elseif proc=="approx_permutation"
        p_values = proc_dispatch(scores, ref_scores, approx_perm_test; tail=tail, centroid=centroid, n=n, filter_scores=filter_scores)
    else
        throw(ArgumentError(("""
        Incorrect argument for proc. Please choose from "Mann_Whitney_U", or "approx_permutation".
            """)))
    end
    return p_values
end

"""
    A wrapper function to compute signature scores, associated p-values and adjusted p-values.
"""
function compute_scores_pvalues(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    ref_object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    proc::String="Mann_Whitney_U"; kw...)
    kw_ss = filter_kw(compute_sig_scores_p; kw...)
    compute_sig_scores_p(object, ref_object; kw_ss...)
    kw_pv = filter_kw(compute_pvalues; kw...)
    compute_pvalues(object, ref_object, proc; kw_pv...)
end

"""
    Identify differential signature scores based on fold changes and adjusted p-values.
"""
function diff_sig_scores(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    object_2::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString};
    cut_off::Float64=1.2, adj_p_thresh::Float64=.01, filter_scores::Bool=true, modify_obj::Bool=false)
    scores = access_result(object, "sig_scores")[:,2:end]
    scores_2 = access_result(object_2, "sig_scores")[:,2:end]
    n_sig = length(object.signatures)
    adj_p = access_result(object, "two_sample_test")[:,3]
    i_n_m = findall(!ismissing, adj_p)
    if !isempty(i_n_m)
        is_diff = Vector{Union{Bool, Missing}}(undef, n_sig)
        for i in i_n_m
            s = Vector{Float64}(scores[i, :])
            s_2 = Vector{Float64}(scores_2[i, :])

            if filter_scores
                s = s[s .> 1e-6] # s = Float64[] if all values of s < 1e-6 
                s_2 = s_2[s_2 .> 1e-6]
            end
            m_s = mean(s) # m_s = NaN if s = Float64[]
            m_s_2 = mean(s_2)
            f_c = maximum([m_s / m_s_2, m_s_2 / m_s])
            is_diff[i] = (adj_p[i] < adj_p_thresh) && (f_c >= cut_off)
        end
    else
        println("Adjusted p-values missing. `is_diff` coerced into a vector of `missing`.")
        is_diff = repeat([missing], n_sig) # Vector{Missing}
    end
    df = DataFrame(signature=access_result(object, "sig_scores")[:,1], is_diff=is_diff)

    if modify_obj
        if "is_diff" in keys(object.results)
            object.results["is_diff"][object_2.dataset_id] = df
        else
            object.results["is_diff"] = Dict{String, DataFrame}(object_2.dataset_id => df)
        end
    end
    
    return df
end

"""
    Calculate fold changes of signature scores between two objects.
"""
function calc_fold_changes(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    object_2::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString};
    filter_scores::Bool=true, modify_obj::Bool=false)
    scores = access_result(object, "sig_scores")[:,2:end]
    scores_2 = access_result(object_2, "sig_scores")[:,2:end]
    n_sig = length(object.signatures)

    f_c = Vector{Union{Float64, Missing}}(undef, n_sig)
    for i in 1:n_sig
        s = Vector{Union{Float64, Missing}}(scores[i, :])
        s_2 = Vector{Union{Float64, Missing}}(scores_2[i, :])
        if test_missing(s, s_2)
            f_c[i] = missing
        else
            if filter_scores
                s = s[s .> 1e-6]
                s_2 = s_2[s_2 .> 1e-6]
            end
            m_s = mean(s)
            m_s_2 = mean(s_2)
            f_c[i] = m_s / m_s_2
        end
    end

    if modify_obj
        df = DataFrame(signature=access_result(object, "sig_scores")[:,1], fold_change=f_c)
        if "fold_change" in keys(object.results)
            object.results["fold_change"][object_2.dataset_id] = df
        else
            object.results["fold_change"] = Dict{String, DataFrame}(object_2.dataset_id => df)
        end
    end

    return f_c
end

"""
    Transform fold chages to values greater than 1, i.e.,
    replace the fold change with its reciprocal if it is less than 1.
"""
function transform_f_c(f_c::Vector{Union{Missing, Float64}})
    f_c_new = copy(f_c)
    for i in eachindex(f_c_new)
        if f_c_new[i] < 1
            f_c_new[i] = 1.0 / f_c_new[i]
        end
    end
    return f_c_new
end

"""
    Calculate maximum fold change of signature scores.
"""
function max_fold_change(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    object_2::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString})
    res = access_result(obj, "sig_scores")
    res_2 = access_result(obj_2, "sig_scores")
    adj_p = access_result(obj, "two_sample_test")[:,3]
    i_n_m = findall(!ismissing, adj_p)
    res_f = res[i_n_m, :]
    res_2_f = res_2[i_n_m, :]
    f_c = zeros(Float64, length(i_n_m))
    for i in eachindex(f_c)
        s = Vector{Float64}(res_f[i,2:end])
        s_2 = Vector{Float64}(res_2_f[i,2:end])
        f_c[i] = maximum([mean(s)/mean(s_2), mean(s_2)/mean(s)])
    end
    return maximum(f_c)
end

"""
    Calculate proportion of differential signature scores.
"""
function prop_diff_scores(diff_scores::DataFrame)
    is_diff = diff_scores[:,2]
    if !all(ismissing, is_diff)
        p::Union{Float64, Missing} = sum(skipmissing(is_diff)) / (length(is_diff) - count(ismissing, is_diff))
    else
        p = missing
    end
    return p
end

"""
    Access the result using the key of the results field.

Parameters
----------
object
    A scSigScores object
key
    A key of the results field to specify the DataFrame of interest.
    Choose from "sig_scores", "two_sample_test", "fold_change", and "is_diff".

Output
A DataFrame
"""
function access_result(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    key::T) where {T<:AbstractString}
    k = unique(keys(object.results[key]))[1]
    df = object.results[key][k]
    return df
end

"""
    Access the result using two keys.

Parameters
----------
object
    A scSigScores object
key
    A key of the results field to specify the DataFrame of interest.
    Choose from "sig_scores" or "two_sample_test".
key_2
    A dataset ID as the key for the nested dictionary object.

Output
A DataFrame
"""
function access_result_k(object::scSigScores{<:AbstractString, <:AbstractString, <:AbstractString},
    key::Tu, key_2::Tv) where {Tu<:AbstractString, Tv<:AbstractString}
    df = object.results[key][key_2]
    return df
end

"""
    Select required keyword arguments for the function.
"""
function filter_kw(func; kw...)
    kw_f = Base.kwarg_decl.(methods(func))
    kw_sub = filter(x -> x[1] in kw_f[1], kw)
    return kw_sub
end

"""
    Test if the proportion of signature genes in the data is above the threshold. Genes with all
    zero expression values across cells should be eliminated before using the function.
"""
function test_signature(genes::Vector{<:AbstractString}, signature::Vector{<:AbstractString};
    sig_thresh::T=.75) where {T<:AbstractFloat}
    test = length(intersect(genes, signature)) / length(signature) >= sig_thresh
    return test
end

"""
    Create a vector of boolean values to select genes in the signature.
"""
function select_genes(genes::Vector{<:AbstractString}, signature::Vector{<:AbstractString})
    select = zeros(Bool, length(genes))
    for i in eachindex(genes)
        select[i] = in(genes[i], signature)
    end
    return select
end

"""
    Count number of genes with expression values above a threshold for each cell.
"""
function count_exp_genes(data::SparseMatrixCSC{<:Real, <:Integer};
    exp_thresh::Union{<:AbstractFloat, <:Integer}=0)
    n_c = size(data)[2]
    counts = zeros(Int32, n_c)
    for i in 1:n_c
        col = data[:, i]
        counts[i] = length(findall(x -> x > exp_thresh, col))
    end
    return counts
end

"""
    Perform hypothesis tests over signatures using the specified procedure.
"""
function proc_dispatch(scores::Matrix{Union{Float64, Missing}}, ref_scores::Matrix{Union{Float64, Missing}},
    func::Function; filter_scores::Bool=true, kw...)
    p_values = Vector{Union{Float64, Missing}}(undef, size(scores)[1])
    for i in eachindex(p_values)
        s = scores[i, :]
        r_s = ref_scores[i, :]
        if test_missing(s, r_s)
            p_values[i] = missing
            println("Scores missing for the signature in either group. Skip the hypothesis test.")
        else
            s = Vector{Float64}(s)
            r_s = Vector{Float64}(r_s)
            if filter_scores
                s = s[s .> 1e-6]
                r_s = r_s[r_s .> 1e-6]
            end

            if isempty(s) || isempty(r_s)
                p_values[i] = 1
            else
                p_values[i] = func(s, r_s; kw...)
            end
        end
    end
    return p_values
end

"""
    Test if missing is in the vector.
"""
function test_missing(scores::Vector{Union{Float64, Missing}}, ref_scores::Vector{Union{Float64, Missing}})
    test = any(ismissing, scores) || any(ismissing, ref_scores)
    return test
end

"""
    Insert missing to adj_p with same indices as in p_values.
"""
function insert_missing(adj_p::Vector{Float64}, p_values::Vector{Union{Float64, Missing}})
    i_m = findall(ismissing, p_values)
    a_p_new = Vector{Union{Float64, Missing}}(undef, length(adj_p))
    for i in eachindex(a_p_new)
        a_p_new[i] = adj_p[i]
    end
    for j in i_m
        insert!(a_p_new, j, missing)
    end
    return a_p_new
end

"""
    Perform the Mann-Whitney U test.
"""
function mwu_test(scores::Vector{Float64}, ref_scores::Vector{Float64};
    tail::String="both")
    p_value = pvalue(MannWhitneyUTest(scores, ref_scores); tail=Symbol(tail))
    return p_value
end

"""
    Perform the approximate permutation test.
"""
function approx_perm_test(scores::Vector{Float64}, ref_scores::Vector{Float64};
    tail::String="both", centroid::String="median", n::T=10000) where {T<:Int}
    if centroid=="median"
        f = median
    elseif centroid=="mean"
        f = mean
    else
        throw(ArgumentError(("""
        Incorrect argument for centroid. Please choose from "median" or "mean".
        """)))
    end
    p_value = pvalue(ApproximatePermutationTest(scores, ref_scores, f, n); tail=Symbol(tail))
    return p_value
end

"""
    Adjust p-values for multiple comparisions.
"""
function adjust_p_values(p_values::Vector{Union{Float64, Missing}}; adjust_p::String="Bonferroni")
    p_no_missing = collect(skipmissing(p_values)) # Vector{Float64}
    if !isempty(p_no_missing)
        if adjust_p=="Bonferroni"
            adj_p = adjust(p_no_missing, Bonferroni())
        elseif adjust_p=="BenjaminiHochberg"
            adj_p = adjust(p_no_missing, BenjaminiHochberg())
        elseif adjust_p=="BenjaminiYekutieli"
            adj_p = adjust(p_no_missing, BenjaminiYekutieli())
        elseif adjust_p=="Hochberg"
            adj_p = adjust(p_no_missing, Hochberg())
        elseif adjust_p=="Holm"
            adj_p = adjust(p_no_missing, Holm())
        elseif adjust_p=="Hommel"
            adj_p = adjust(p_no_missing, Hommel())
        elseif adjust_p=="Sidak"
            adj_p = adjust(p_no_missing, Sidak())
        elseif adjust_p=="forward_stop"
            adj_p = adjust(p_no_missing, ForwardStop())
        else
            throw(ArgumentError(("""
            Incorrect argument for adjust_p.
            Please choose from "Bonferroni", "BenjaminiHochberg", "BenjaminiYekutieli", "Hochberg",
            "Holm", "Hommel", "Sidak", or "forward_stop".
                """)))
        end
    else
        adj_p = p_values # vector of all missing values
    end
    return adj_p
end
