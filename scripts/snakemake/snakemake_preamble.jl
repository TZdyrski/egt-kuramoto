# SPDX-License-Identifier: GPL-3.0-or-later

function parse_string(string::String)::Union{String,Int64,Float64}
  value = tryparse(Float64,string)
  if isnothing(value)
    return string
  end
  return value
end

# Note: snakemake object is prepopulated by snakemake executable
wildcards_dict = copy(snakemake.wildcards)
# Check for and remove optional wildcards
if haskey(wildcards_dict, "nb_players_flag") && pop!(wildcards_dict, "nb_players_flag") == ""
  delete!(wildcards_dict, "nb_players")
end
if haskey(wildcards_dict, "covariance_cutoff_fraction_flag") && pop!(wildcards_dict, "covariance_cutoff_fraction_flag") == ""
  delete!(wildcards_dict, "covariance_cutoff_fraction")
end
if haskey(wildcards_dict, "community_resolution_flag") && pop!(wildcards_dict, "community_resolution_flag") == ""
  delete!(wildcards_dict, "community_resolution")
end
if haskey(wildcards_dict, "community_beta_flag") && pop!(wildcards_dict, "community_beta_flag") == ""
  delete!(wildcards_dict, "community_beta")
end
if haskey(wildcards_dict, "community_n_iter_flag") && pop!(wildcards_dict, "community_n_iter_flag") == ""
  delete!(wildcards_dict, "community_n_iter")
end
if haskey(wildcards_dict, "walktrap_steps_flag") && pop!(wildcards_dict, "walktrap_steps_flag") == ""
  delete!(wildcards_dict, "walktrap_steps")
end
if haskey(wildcards_dict, "early_cutoff_fraction_flag") && pop!(wildcards_dict, "early_cutoff_fraction_flag") == ""
  delete!(wildcards_dict, "early_cutoff_fraction")
end
if haskey(wildcards_dict, "seed_flag") && pop!(wildcards_dict, "seed_flag") == ""
  delete!(wildcards_dict, "seed")
end

# snakemake.wildcards doubles the wildcards: each wildcard has a key
# with the (string) name of the wildcard as well as a key with the (integer) position
# Only keep the (string) name key
filter!(x -> isa(x.first, String), wildcards_dict)

wildcards_dict = Dict((Symbol(k),parse_string(v)) for (k,v) in wildcards_dict)

integer_params = [:time_steps, :nb_phases, :nb_players, :time_step, :decimation_factor, :community_n_iter, :seed, :walktrap_steps]
wildcards_dict = Dict((k, k in integer_params ? Integer(v) : v) for (k,v) in wildcards_dict)
bool_params = [:only_mixed_games]
wildcards_dict = Dict((k, k in bool_params ? parse(Bool,v) : v) for (k,v) in wildcards_dict)

wildcards = NamedTuple(wildcards_dict)
