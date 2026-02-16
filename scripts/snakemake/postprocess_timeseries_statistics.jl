# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "postprocess.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
postprocessKeys = [:only_mixed_games]
postprocessDict = Dict(k => pop!(loadDict,k) for k in keys(loadDict) if k in postprocessKeys)
numSeeds = pop!(loadDict, :num_seeds)
early_cutoff_fraction = pop!(loadDict, :early_cutoff_fraction, nothing)
num_samples = pop!(loadDict, :num_samples, nothing)

# Get data
adj_matrix = get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1]

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "timeseries", savename(loadDict, "jld2")))
processing_fun = dataDict -> calc_timeseries_statistics(dataDict["initial_actions"],
	dataDict["deltas"],
	loadDict[:nb_phases], loadDict[:symmetry_breaking],
	loadDict[:B_to_c], loadDict[:beta_to_B], loadDict[:cost], adj_matrix;
	postprocessDict...)
insert_time_fun = df -> insertcols(df, 1, :time => 0:loadDict[:time_steps])
data_generation_fun = insert_time_fun ∘ DataFrame ∘ processing_fun ∘ loading_fun
rename_col = df -> rename(df, :most_common_game_types => :game_type, :order_parameters => :order_parameter, :fraction_communicative => :communicative_fraction)
average_fun_games = df -> combine(groupby(df,:time), :game_type => mode ∘ sort => :game_type, Not([:game_type,:seed,:time]) .=> mean, Not([:game_type,:seed,:time]) .=> std)
apply_seed_avg = (generation_function,numSeeds) -> (average_fun_games ∘ rename_col ∘ apply_over_param(:seed => 1:numSeeds)(generation_function))

# Run code
result = loadDict |> apply_seed_avg(data_generation_fun, numSeeds)

if !isnothing(early_cutoff_fraction)
	# Note: the populations include the initial statistics, so we need one more than time-steps
	new_end_time = Int(round(loadDict[:time_steps]*early_cutoff_fraction))+1
	result = result[1:new_end_time,:]
end

if isnothing(num_samples)
	result = downsample_data(result)
else
	result = downsample_data(result, num_samples)
end

# Write out data
mkpath(datadir("processed", "timeseries_statistics"))
CSV.write(datadir("processed","timeseries_statistics",savename(wildcards,"csv")), result)
