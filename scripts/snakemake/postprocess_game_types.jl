# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
using StatsBase
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "postprocess.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
postprocessKeys = [:only_mixed_games]
postprocessDict = Dict(k => pop!(loadDict,k) for k in keys(loadDict) if k in postprocessKeys)
numSeeds = pop!(loadDict, :num_seeds)

# Get data
adj_matrix = get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1]
# Make list of games following Bruns 2015 from top-right to bottom-left
game_ordering = ("asymmetry", string.(instances(GameType))...)

# Define processing functions
loading_fun = loadDict -> merge(wload(datadir("raw", "timeseries", savename(loadDict, "jld2"))),
		Dict("symmetry_breaking" => loadDict[:symmetry_breaking]))
processing_fun = dataDict -> calc_timeseries_statistics(dataDict["initial_actions"],
	dataDict["deltas"],
	loadDict[:nb_phases], dataDict["symmetry_breaking"],
	loadDict[:B_to_c], loadDict[:beta_to_B], loadDict[:cost], adj_matrix;
	postprocessDict...)
game_proportion_fun = (dict -> Dict(string(k) => v for (k,v) in pairs(dict))) ∘ proportionmap ∘ (df -> df["most_common_game_types"])
data_generation_fun = DataFrame ∘ game_proportion_fun ∘ processing_fun ∘ loading_fun
sort_game_types_fun = dataFrame -> select(dataFrame, intersect(("asymmetry",game_ordering.*"_mean"...,game_ordering.*"_std"...), names(dataFrame)))

# Run code
result = loadDict |> apply_all_asymm(apply_seed_avg(data_generation_fun, numSeeds)) |> sort_game_types_fun

# Write out data
mkpath(datadir("processed", "gametype"))
CSV.write(datadir("processed","gametype",savename(wildcards,"csv")), result)
