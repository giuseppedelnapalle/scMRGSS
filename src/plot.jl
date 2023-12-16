#!/usr/bin/env julia
#= 
functions to plot scSigScores output
=#

# using CairoMakie
# using DataFrames: DataFrame
# using Distances: pairwise, Euclidean
# using Clustering: hclust
# using Random: seed!, randperm

"""
    Plot the numeric values in a DataFrame on a heatmap.
"""
function heatmap_df(df::DataFrame, path_to_file::String, seed::Int64=6358;
    resolution::Tuple{Int64, Int64}=(1200, 1000), fontsize::Int64=20, label::String="value", title::String="",
    hclust_x=false, hclust_y=false, x_lab=true, y_lab=true, linkage::String="average", text::Bool=false,
    colormap::Symbol=:plasma, pt_per_unit::Int64=2, width_colorbar::Int64=15, ticksize::Int64=15)
    fig = Figure(resolution = resolution, fontsize = fontsize)
    df_r_o = reorder_df_hclust(df, seed; hclust_x=hclust_x, hclust_y=hclust_y, linkage=linkage)
    mat = rotate_matrix_cw(Matrix{Float64}(df_r_o[:, 2:end]))
    m, n = size(mat)
    xticks = names(df_r_o)[2:end]
    yticks = reverse(df_r_o[:, 1])

    ax = get_axis(fig, m, n, xticks, yticks; x_lab=x_lab, y_lab=y_lab)
    hmap = heatmap!(ax, mat, colormap = colormap)
    if text
        thresh = (maximum(mat) - minimum(mat)) / 6
        for i in 1:m, j in 1:n
            txtcolor = mat[i, j] < thresh ? :white : :black
            text!(ax, "$(round(mat[i,j], digits = 2))", position = (i, j),
                color = txtcolor, align = (:center, :center))
        end
    end
    Colorbar(fig[1, 2], hmap; label = label, width = width_colorbar, ticksize = ticksize)
    ax.xticklabelrotation = π / 3
    ax.xticklabelalign = (:right, :center)
    if !isempty(title)
        ax.title = title
    end
    
    save(path_to_file, fig; pt_per_unit = pt_per_unit)
end

"""
    Reorder rows and/or columns of a DataFrame by hierarchical clustering.
"""
function reorder_df_hclust(df::DataFrame, seed::Int64=8369;
    hclust_x=false, hclust_y=false, linkage::String="average")
    seed!(seed)
    s_1, s_2 = rand(1:10_000, 2)
    lkg = Symbol(linkage)
    r_nm = df[:, 1]
    c_nm = names(df)[2:end]
    mat = Matrix{Float64}(df[:, 2:end])
    m, n = size(mat)

    if hclust_x || hclust_y
        if hclust_x
            seed!(s_1)
            perm = randperm(n)
            mat = mat[:, perm]
            c_nm = c_nm[perm]
            dis_mat = pairwise(Euclidean(), mat, dims=2) # col
            hcl = hclust(dis_mat, linkage=lkg)
            mat = mat[:, hcl.order]
            c_nm = c_nm[hcl.order]
        end
        if hclust_y
            seed!(s_2)
            perm = randperm(m)
            mat = mat[perm, :]
            r_nm = r_nm[perm]
            dis_mat = pairwise(Euclidean(), mat, dims=1) # row
            hcl = hclust(dis_mat, linkage=lkg)
            mat = mat[hcl.order, :]
            r_nm = r_nm[hcl.order]
        end
        df_r_o = hcat(DataFrame(names(df)[1] => r_nm), DataFrame(mat, c_nm))
    else
        df_r_o = df
        println("hierarchical clustering along both axes not performed. Original df returned.")
    end

    return df_r_o
end

"""
    Rotate a matrix by 90 degrees clockwise.
"""
function rotate_matrix_cw(mat::Matrix{Float64})
    mat_r = transpose(mat)
    for r in eachrow(mat_r)
        reverse!(r)
    end
    return mat_r
end

"""
    Obtain the 2D axis for the heatmap.
"""
function get_axis(fig::Figure, m::Int64, n::Int64, xticks::Vector{String}, yticks::Vector{String};
    x_lab::Bool=true, y_lab::Bool=false)
    if x_lab
        if y_lab
            ax = Axis(fig[1, 1], xticks = (1:m, xticks), yticks = (1:n, yticks))
        else
            ax = Axis(fig[1, 1], xticks = (1:m, xticks))
        end
    else
        if y_lab
            ax = Axis(fig[1, 1], yticks = (1:n, yticks))
        else
            ax = Axis(fig[1, 1])
        end
    end
    return ax
end
