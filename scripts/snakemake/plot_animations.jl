# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
using Graphs
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "plotting.jl"))
include(srcdir("julia", "utils.jl"))

# Run code
dataDict = wload(datadir("raw", "timeseries", savename(wildcards,"jld2")))
graph = SimpleDiGraph(get_adj_matrices(;adj_matrix_source=wildcards.adj_matrix_source)[1])
result = plot_graph_evolution(dataDict, graph)

# Write out data
mkpath(plotsdir("animations"))
save(plotsdir("animations",savename(wildcards,"mp4")), result)
