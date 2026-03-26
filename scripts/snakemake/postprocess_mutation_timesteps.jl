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

# Define processing functions
loading_fun = loadDict -> wload(datadir("raw", "timeseries", savename(loadDict, "jld2")))
processing_fun = dataDict -> Dict("mutation_timesteps" => dataDict["steps_following_mutation"])
data_generation_fun = DataFrame ∘ processing_fun ∘ loading_fun

# Run code
result = loadDict |> data_generation_fun

# Write out data
mkpath(datadir("processed", "mutation_timesteps"))
CSV.write(datadir("processed","mutation_timesteps",savename(wildcards,"csv")), result)
