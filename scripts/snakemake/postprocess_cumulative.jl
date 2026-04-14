# SPDX-License-Identifier: GPL-3.0-or-later

# Workaround to force snakemake to use Project.toml
# Source: https://github.com/snakemake/snakemake/issues/2215#issuecomment-1747136802
dirname(Base.active_project()) != pwd() && exit(run(`julia --project=@. $(@__FILE__)`).exitcode)

# Creates wildcards NamedTuple with snakemake wildcards
using DrWatson
@quickactivate "Chimera_EGT_Kuramoto"
include(scriptsdir("snakemake","snakemake_preamble.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
numSeeds = pop!(loadDict, :num_seeds)

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "cumulative", savename(loadDict, "jld2")))
processing_fun = dataDict -> Dict("communicative_fraction" => dataDict["fraction_communicative"],
	"B0" => dataDict["Bs"])
data_generation_fun = DataFrame ∘ processing_fun ∘ loading_fun
average_fun = df -> combine(groupby(df,:B0), Not(:seed, :B0) .=> mean, Not(:seed,:B0) .=> std)
apply_seed_avg = (generation_function,numSeeds) -> (average_fun ∘ apply_over_param(:seed => 1:numSeeds)(generation_function))

# Run code
result = loadDict |> apply_seed_avg(data_generation_fun, numSeeds)

# Write out data
mkpath(datadir("processed", "cumulative"))
CSV.write(datadir("processed","cumulative",savename(wildcards,"csv")), result)
