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
include(srcdir("julia", "utils.jl"))

# Get parameter sets
loadDict = Dict(pairs(wildcards))
adj_matrix_source = pop!(loadDict,:adj_matrix_source)

# Get data
graph = SimpleWeightedDiGraph(get_adj_matrices(;adj_matrix_source)[1])

# Define processing functions
processing_fun = dataDict -> communicative_fraction_theory(graph; dataDict...)
rename_col = df -> rename(df, :Bs => :B0)
data_generation_fun = rename_col ∘ DataFrame ∘ processing_fun

# Run code
result = loadDict |> data_generation_fun

# Write out data
mkpath(datadir("processed", "cumulative_theory"))
CSV.write(datadir("processed","cumulative_theory",savename(wildcards,"csv")), result)
