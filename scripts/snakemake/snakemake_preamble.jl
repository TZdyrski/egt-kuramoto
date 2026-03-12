# SPDX-License-Identifier: GPL-3.0-or-later

# Load SparseArrays since all scripts either save or load
# a dataDict containing a SparseArray
using SparseArrays
using CSV
using DataFrames
using StatsBase

function parse_string(string::String)::Union{String,Int64,Float64}
  value = tryparse(Float64,string)
  if isnothing(value)
    return string
  end
  return value
end

# Note: snakemake object is prepopulated by snakemake executable
wildcards_dict = deepcopy(snakemake.wildcards)

# snakemake.wildcards doubles the wildcards: each wildcard has a key
# with the (string) name of the wildcard as well as a key with the (integer) position
# Only keep the (string) name key
filter!(x -> isa(x.first, String), wildcards_dict)

# Remove empty optional arguments (denoted by _flag)
optional_flags = filter(endswith("_flag"), keys(wildcards_dict))
for flag in optional_flags
	if pop!(wildcards_dict, flag) == ""
		# Remove last 5 character ("_flag") to get parameter name
		param_name = chop(flag, tail=5)
		delete!(wildcards_dict, param_name)
	end
end

wildcards_dict = Dict((Symbol(k),parse_string(v)) for (k,v) in wildcards_dict)

integer_params = [:time_steps, :time_snapshot, :nb_phases, :nb_players, :time_step, :decimation_factor, :community_n_iter, :num_seeds, :seed, :walktrap_steps]
wildcards_dict = Dict((k, k in integer_params ? Integer(v) : v) for (k,v) in wildcards_dict)
bool_params = [:only_mixed_games]
wildcards_dict = Dict((k, k in bool_params ? parse(Bool,v) : v) for (k,v) in wildcards_dict)

wildcards = NamedTuple(wildcards_dict)

# Define function that applies processing_function over loadDict for
# each value of paramPair.first in paramPair.second
function apply_over_param(paramPair::Pair, loadDict::Dict, processing_function::Function)
	results = DataFrame()
	for val = paramPair.second
		additionalConfig = Dict(paramPair.first => val)
		result = processing_function(merge(loadDict, additionalConfig))
		insertcols!(result, 1, additionalConfig...)
		append!(results, result, cols=:union)
	end
	# Replace any Missing with 0.0
	results = coalesce.(results, 0.0)
	return results
end

# Curry apply_over_param
function apply_over_param(paramPair::Pair)
	return processing_function -> (loadDict -> apply_over_param(paramPair, loadDict, processing_function))
end

symmetry_breaking_range = [0.0,0.25,0.5,0.75,1.0]

rename_asymmetry = df -> rename(df, :symmetry_breaking => :asymmetry)
average_fun = df -> combine(df, Not(:seed) .=> mean, Not(:seed) .=> std)

# Define convenience functions
apply_all_asymm = generation_function -> (rename_asymmetry ∘ apply_over_param(:symmetry_breaking => symmetry_breaking_range)(generation_function))
apply_seed_avg = (generation_function,numSeeds) -> (average_fun ∘ apply_over_param(:seed => 1:numSeeds)(generation_function))
