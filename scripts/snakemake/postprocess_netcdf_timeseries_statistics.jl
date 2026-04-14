# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
using DimensionalData
using DimensionalData.Lookups
using FillArrays
using YAXArrays
using NetCDF
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

include(srcdir("julia", "postprocess.jl"))

# Get list of input data filenames
input_filenames = filter(endswith(".jld2"), collect(values(snakemake.input)))

# If the first parameter in the filename is multiple words joined by an
# underscore, parse_savename strips the first word;
# Therefore, add an extra underscore to prevent this
strip_prefix_and_add_underscore = filename -> "_"*splitpath(filename)[end]

load_dir = allequal(dirname.(input_filenames)) ? dirname(input_filenames[1]) : throw(ArgumentError("All input .jld2 files must share same directory"))

# Convert filenames to parameter dictionaries
param_dicts = input_filenames .|> strip_prefix_and_add_underscore .|> filename -> parse_savename(filename)[2]

# Get list of parameters that are the same between all inputs
properties = Dict(intersect(param_dicts...))

# Remove properties from each dict
unique_params = setdiff.(param_dicts, Ref(properties))

# Merge unique parameters into single vectors
param_dict = mergewith(vcat, Dict.(unique_params)...)

# Remove duplicate param values
param_dict = Dict(k => unique(v) for (k,v) in param_dict)

# Remove seed
numSeeds = length(pop!(param_dict, "seed"))

# Get data
adj_matrix = get_adj_matrices(;adj_matrix_source=properties["adj_matrix_source"])[1]

# Define processing functions
loading_fun = (loadDict,properties) -> merge(wload(datadir("raw", "timeseries", savename(loadDict, "jld2"))),
	loadDict,properties)
	processing_fun = dataDict -> merge(calc_timeseries_statistics(dataDict["initial_actions"],
		dataDict["deltas"],
		dataDict["nb_phases"], dataDict[:symmetry_breaking],
		dataDict[:B_to_c], dataDict["beta_to_B"], dataDict["cost"], adj_matrix;
		only_mixed_games=false),
	Dict("most_common_game_types_only_mixed_games" => calc_timeseries_statistics(dataDict["initial_actions"],
		dataDict["deltas"],
		dataDict["nb_phases"], dataDict[:symmetry_breaking],
		dataDict[:B_to_c], dataDict["beta_to_B"], dataDict["cost"], adj_matrix;
		only_mixed_games=true)["most_common_game_types"]),
	Dict("time" => 0:dataDict["time_steps"]),
	)
convert_game_types_fun = df -> transform(df,
	:most_common_game_types => ByRow(Integer) => :most_common_game_types,
	:most_common_game_types_only_mixed_games => ByRow(Integer) => :most_common_game_types_only_mixed_games,
	)
data_generation_fun = convert_game_types_fun ∘ DataFrame ∘ processing_fun ∘ loading_fun
data_generation_fun_curried = properties -> (loadDict -> data_generation_fun(loadDict, properties))
average_fun_games = df -> combine(groupby(df,:time),
  :most_common_game_types => mode ∘ sort => :most_common_game_types,
  :most_common_game_types_only_mixed_games => mode ∘ sort => :most_common_game_types_only_mixed_games,
	)
apply_seed_avg = (generation_function,numSeeds) -> (average_fun_games ∘ apply_over_param(:seed => 1:numSeeds)(generation_function))

loading_fun = (loadDict,properties) -> merge(wload(datadir("raw", "timeseries", savename(loadDict, "jld2"))),
	loadDict,properties)
	processing_fun = dataDict -> merge(calc_timeseries_statistics(dataDict["initial_actions"],
		dataDict["deltas"],
		dataDict["nb_phases"], dataDict[:symmetry_breaking],
		dataDict[:B_to_c], dataDict["beta_to_B"], dataDict["cost"], adj_matrix),
	)
filter_games_fun = df -> select(df, Not(:most_common_game_types))
data_generation_fun_single = filter_games_fun ∘ DataFrame ∘ processing_fun ∘ loading_fun
data_generation_fun_curried_single = properties -> (loadDict -> data_generation_fun_single(loadDict, properties))

# Convert to spans
function convert_to_span(vec::Vector{<:Number})
	deltas = unique(diff(vec))
	if length(deltas) == 1
		return Sampled(vec; span=Regular(deltas[1]))
	end
	return Sampled(sort(vec))
end
convert_to_span(vec::Vector) = Sampled(sort(vec))
ds_axes_base = Tuple(Dim{Symbol(k)}(convert_to_span(v)) for (k,v) in param_dict)
times_downsampled = downsample_data(DataFrame(:time => 0:properties["time_steps"])).time
time_axis = Dim{:time_step}(times_downsampled)
ds_axes = (time_axis, ds_axes_base...)
player_axis = Dim{:player}(1:size(adj_matrix,1))
ds_axes_player = (time_axis, player_axis, ds_axes_base...)

# Add lookup dictionary
game_lookup_dict = Dict("enum_lookup:" .* string.(Integer.(instances(GameType)))
                     .=> string.(Symbol.(instances(GameType))))

# Generate data skeleton
skeleton = YAXArray(ds_axes, Fill{Union{Missing,Float64}}(missing, length.(ds_axes)))
# Note: it seems using a properties dict and an array is Missing data
# types triggers a bug, throwing
# "MethodError: Cannot `convert` an object of type Missing to an object of type Int64"
# Instead, use -1 as a missing value
skeleton_games = YAXArray(ds_axes, Fill{Union{Missing,Int64}}(missing, length.(ds_axes)))
dataset_names = ["order_parameters","fraction_communicative"]
dataset_names_games =  ["most_common_game_types","most_common_game_types_only_mixed_games"]
skeleton_ds = Dataset(; properties=merge(properties,game_lookup_dict),
	Dict(Symbol(name) => skeleton for name in dataset_names)...,
	Dict(Symbol(name) => skeleton_games for name in dataset_names_games)...,
	:timeseries => YAXArray(ds_axes_player, Fill{Union{Missing,Int64}}(missing, length.(ds_axes_player)))
	)
# Write out skeleton data
mkpath(datadir("processed", "netcdf"))
file_savename = datadir("processed","netcdf",savename("timeseries_statistics",wildcards,"nc"))
savedataset(skeleton_ds; path=file_savename, driver=:netcdf, skeleton=true, overwrite=true, compress=9)

# Lazily load skeleton data
ds_array = open_dataset(file_savename)

# Populate skeleton data
param_names = Symbol.(keys(param_dict))
param_dims = dims(ds_array.fraction_communicative, param_names...)
for vals in Iterators.product(values(param_dict)...)
	param_set = Dict(param_names .=> vals)
	loadDict = merge(param_set,properties)

	# Run code
	results_avg = loadDict |> apply_seed_avg(data_generation_fun_curried(properties), numSeeds)
	results_avg = downsample_data(results_avg)

	results_single = merge(loadDict,Dict("seed"=>1)) |> data_generation_fun_curried_single(properties)
	results_single = downsample_data(results_single)

	results = hcat(results_avg, results_single)

  # Use `At` to specify exact value match
	params_at = (rebuild(k,At(v)) for (k,v) in zip(param_dims,vals))

	for variable in setdiff(ds_array.cubes.keys, [:timeseries])
	  ds_array[variable][params_at...] .= results[!,variable]
	end

	dataDict = loading_fun(merge(loadDict,Dict("seed"=>1)), properties)
	timeseries_data = transpose(decode_delta_encoded_all(dataDict["initial_actions"],dataDict["deltas"]))
	timeseries_data_downsampled = Matrix(downsample_data(DataFrame(timeseries_data, :auto)))
	ds_array[:timeseries][params_at...] .= timeseries_data_downsampled
end
