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
early_cutoff_fraction = pop!(loadDict, :early_cutoff_fraction, nothing)
num_samples = pop!(loadDict, :num_samples, nothing)

# Get data
graph = SimpleWeightedDiGraph(get_adj_matrices(;adj_matrix_source=loadDict[:adj_matrix_source])[1])
communities = generate_communities(graph, pop!(postprocessDict, :community_algorithm);
                                   postprocessDict...)

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "timeseries", savename(loadDict, "jld2")))
processing_fun = dataDict -> order_parameters_by_community(dataDict["initial_actions"], dataDict["deltas"], communities, loadDict[:nb_phases])
insert_time_fun = df -> insertcols(df, 1, :time => 0:loadDict[:time_steps])
array_cols_to_dict_fun = array -> Dict("orderParamGroup"*string(n) => v for (n,v) in enumerate(eachcol(array)))
data_generation_fun = insert_time_fun ∘ DataFrame ∘ array_cols_to_dict_fun ∘ processing_fun ∘ loading_fun

# Run code
result = loadDict |> data_generation_fun

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
mkpath(datadir("processed", "order_parameter"))
# Shorten community_ keywords to avoid file name length limit
wildcards_dict = Dict(pairs(wildcards))
wildcards_dict = Dict(replace(string(k), "community_" => "comm_") => v for (k,v) in wildcards_dict)
CSV.write(datadir("processed","order_parameter",savename(wildcards_dict,"csv")), result)
