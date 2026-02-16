# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "utils.jl"))
include(srcdir("julia", "postprocess.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
time_step = pop!(loadDict, :time_step)

# Get data
graph = SimpleDiGraph(get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1])

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "timeseries", savename(loadDict, "jld2")))
processing_fun = dataDict -> generate_nodes_edges(graph; data=dataDict, time_step)[1]
data_generation_fun = processing_fun ∘ loading_fun

# Run code
result = loadDict |> data_generation_fun

# Write out data
mkpath(datadir("processed", "graph_structure_snapshot"))
CSV.write(datadir("processed","graph_structure_snapshot",savename("vertices", wildcards,"csv")), result)
