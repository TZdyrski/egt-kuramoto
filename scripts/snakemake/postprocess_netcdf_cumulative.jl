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

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "cumulative", savename(loadDict, "jld2")))
processing_fun = dataDict -> Dict("communicative_fraction" => dataDict["fraction_communicative"],
	"B0" => dataDict["Bs"])
data_generation_fun = DataFrame ∘ processing_fun ∘ loading_fun
average_fun = df -> combine(groupby(df,:B0), Not(:seed, :B0) .=> mean, Not(:seed,:B0) .=> std)
apply_seed_avg = (generation_function,numSeeds) -> (average_fun ∘ apply_over_param(:seed => 1:numSeeds)(generation_function))

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

# Convert to spans
function convert_to_span(vec::Vector{<:Number})
	deltas = unique(diff(vec))
	if length(deltas) == 1
		return Sampled(vec; span=Regular(deltas[1]))
	end
	return Sampled(sort(vec))
end
convert_to_span(vec::Vector) = Sampled(sort(vec))
ds_axes = Tuple(Dim{Symbol(k)}(convert_to_span(v)) for (k,v) in param_dict)

# Add axis for maximum_joint_benefit
default_nb_players = 20
nb_players = get(properties, "nb_players", default_nb_players)
B_crit = 2 * properties["cost"] * (nb_players - 1) / (nb_players - 2)
nb_Bs = 11
step_size_Bs = 0.04
Bs = B_crit .+ ((0:(nb_Bs - 1)) .- (nb_Bs - 1) / 2) .* step_size_Bs
ds_axes = (Dim{:maximum_joint_benefit}(Bs), ds_axes...)

# Generate data skeleton
skeleton = YAXArray(ds_axes, Fill{Union{Missing,Float64}}(missing, length.(ds_axes)))
dataset_names = [name*postfix for postfix in ["_mean","_std"] for name in
		["communicative_fraction"]]
skeleton_ds = Dataset(; properties, Dict(Symbol(name) => skeleton for name in dataset_names)...)

# Write out skeleton data
mkpath(datadir("processed", "netcdf"))
file_savename = datadir("processed","netcdf",savename("cumulative",wildcards,"nc"))
savedataset(skeleton_ds; path=file_savename, driver=:netcdf, skeleton=true, overwrite=true, compress=9)

# Lazily load skeleton data
ds_array = open_dataset(file_savename)

# Populate skeleton data
param_names = Symbol.(keys(param_dict))
param_dims = dims(ds_array.communicative_fraction_mean, param_names...)
for vals in Iterators.product(values(param_dict)...)
	param_set = Dict(param_names .=> vals)
	loadDict = merge(param_set,properties)

	# Run code
	result = loadDict |> apply_seed_avg(data_generation_fun, numSeeds)

  # Use `At` to specify exact value match
	params_at = (rebuild(k,At(v)) for (k,v) in zip(param_dims,vals))

	for variable in ds_array.cubes.keys
	  ds_array[variable][params_at...] .= result[!,variable]
	end
end
