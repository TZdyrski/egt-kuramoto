# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
using Graphs
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "utils.jl"))
include(srcdir("julia", "postprocess.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))

# Get data
graph = SimpleDiGraph(get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1])

# Run code
vertex_coords, edge_list = generate_nodes_edges(graph)

# Write out data
mkpath(datadir("processed", "graph_structure"))
CSV.write(datadir("processed","graph_structure", savename("vertices",wildcards,"csv")), vertex_coords)
CSV.write(datadir("processed","graph_structure", savename("edges",wildcards,"csv")), edge_list)
