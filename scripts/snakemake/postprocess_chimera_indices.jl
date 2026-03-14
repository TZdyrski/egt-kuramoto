# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
using SimpleWeightedGraphs
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "postprocess.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
postprocessKeys = [:community_algorithm, :covariance_cutoff_fraction, :covariance_data,
                   :walktrap_steps, :community_resolution, :community_beta, :community_n_iter]
postprocessDict = Dict(k => pop!(loadDict,k) for k in keys(loadDict) if k in postprocessKeys)
numSeeds = pop!(loadDict, :num_seeds)

# Get data
graph = SimpleWeightedDiGraph(get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1])
communities = generate_communities(graph, pop!(postprocessDict, :community_algorithm);
                                   postprocessDict...)

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "timeseries", savename(loadDict, "jld2")))
processing_fun = dataDict -> get_chimera_indices(dataDict["initial_actions"], dataDict["deltas"], communities, loadDict[:nb_phases])
data_generation_fun = DataFrame ∘ processing_fun ∘ loading_fun

# Run code
result = loadDict |> apply_all_asymm(apply_seed_avg(data_generation_fun, numSeeds))

# Write out data
mkpath(datadir("processed", "chimeraindex"))
# Shorten community_ keywords to avoid file name length limit
wildcards_dict = Dict(pairs(wildcards))
wildcards_dict = Dict(replace(string(k), "community_" => "comm_") => v for (k,v) in wildcards_dict)
CSV.write(datadir("processed","chimeraindex",savename(wildcards_dict,"csv")), result)
