#!/usr/bin/env julia
#= 
compare signature scores between cell types of PBMC from same or distinct datasets
nikola
Jan 2023
=#

# check # threads
Threads.nthreads()

import CSV

# set working directory
const wd = dirname(@__FILE__)
cd(wd)

# result directory
const res_dir = joinpath(wd, "results")
if !isdir(res_dir)
    mkdir(res_dir)
end

# figure directory
const fig_dir = joinpath(wd, "figures")
if !isdir(fig_dir)
    mkdir(fig_dir)
end

# import functions
include("func_comp_sig_scores_pbmc.jl")

# set constants
const study_id = "pbmc"

const fn_sig = "/home/nikola/Documents/Rsch/resources/gene_sets/MSigDB/C2_curated_gene_sets/c2.cp.biocarta.v2023.1.Hs.symbols.gmt"

# cell labels
const labels = ["cd14_mono", "naive_cd4_t", "b"]

# datasets
const dataset_ids = ["pbmc_10k_3p", "pbmc_10k_5p", "TabulaSapiens"]

# const min_features = 2000
const min_features = 2000
const min_cells = 2
const sig_thresh = .6
const fc_thresh = 1.2
const adj_p_thresh = .01

# if signature file not "h.all.v2022.1.Hs.symbols.gmt"
g_s_id = split(split(fn_sig, "/")[10], ".")[3]

# include("func_comp_sig_scores_pbmc.jl")

# calculate proportions of differential signature scores
@time prop = prop_diff_scores_export_res_wf(labels, dataset_ids, fn_sig, res_dir, g_s_id;
min_features = min_features, min_cells = min_cells, sig_thresh = sig_thresh,
fc_thresh = fc_thresh, adj_p_thresh = adj_p_thresh);

println(prop)
# adj_scores
# println(size(adj_scores))

# export the result
lb_p = join([labels]..., "-")
ds_p = join([dataset_ids]..., "-")
CSV.write(joinpath(res_dir, "prop_diff_scores_$(study_id)_$(g_s_id)_$(lb_p)_$(ds_p)_$(sig_thresh)_$(fc_thresh)_$(min_features)_upd_v1.6.7.csv"), prop)

# heatmaps
# signature scores

# proportions of significant differential scores
fn_p_pdf = "heatmap_prop_diff_scores_$(study_id)_$(g_s_id)_$(lb_p)_$(ds_p)_$(sig_thresh)_$(fc_thresh)_$(min_features)_$(min_cells)_upd_v1.6.7.pdf"
heatmap_df(prop, joinpath(fig_dir, fn_p_pdf); label="proportion", text=true)
# no text
fn_p_pdf = "heatmap_prop_diff_scores_$(study_id)_$(g_s_id)_$(lb_p)_$(ds_p)_$(sig_thresh)_$(fc_thresh)_$(min_features)_$(min_cells)_no_text_upd_v1.6.7.pdf"
heatmap_df(prop, joinpath(fig_dir, fn_p_pdf); label="proportion", text=false)

println("task done")
