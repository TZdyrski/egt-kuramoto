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
time_snapshot = Integer(pop!(loadDict, :time_snapshot))

# Get data
interaction_adj_matrix = get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1]

# Define processing functions
loading_fun = loadDict -> merge(wload(datadir("raw", "timeseries", savename(loadDict, "jld2"))),
		Dict("symmetry_breaking" => loadDict[:symmetry_breaking]))
decode_data_fun = dataDict -> merge(dataDict, Dict("state" => decode_delta_encoded(
		dataDict["initial_actions"], dataDict["deltas"], time_snapshot)))
gen_game_type_fun = dataDict -> merge(dataDict, Dict("game_types" => game_types_per_strategy_pair(
		loadDict[:B_to_c]*loadDict[:cost], loadDict[:beta_to_B]*loadDict[:B_to_c]*loadDict[:cost], loadDict[:cost],
		dataDict["symmetry_breaking"], loadDict[:nb_phases]; postprocessDict...)))
processing_fun = dataDict -> count_games(dataDict["state"]; game_types=dataDict["game_types"], interaction_adj_matrix)[2]
remove_missing_fun = filter(p -> !ismissing(p.first))
remove_zeros_fun = filter(p -> p.second != 0)
game_proportion_fun = dataDict -> Dict(k => v/sum(collect(values(dataDict))) for (k,v) in dataDict)
stringify_game_names = dict -> Dict(string(k) => v for (k,v) in dict)
data_generation_fun = DataFrame ∘ stringify_game_names ∘ game_proportion_fun ∘ remove_zeros_fun ∘ remove_missing_fun ∘ processing_fun ∘ gen_game_type_fun ∘ decode_data_fun ∘ loading_fun

# Run code
result = loadDict |> apply_all_asymm(apply_seed_avg(data_generation_fun, numSeeds))

# Write out data
mkpath(datadir("processed", "game_snapshot"))
CSV.write(datadir("processed","game_snapshot",savename(wildcards,"csv")), result)
