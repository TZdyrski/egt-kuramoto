using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Load source files
include(srcdir("julia", "postprocess.jl"))

# Print number of unidirectional and bidirectional edges
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                          "c-elegans-undirected"]
    calc_number_unidirection_bidirectional_edges(; adj_matrix_source=adj_matrix_source)
end

# Save well-mixed graph
export_graph_nodes_edges(; adj_matrix_source="well-mixed")

# Save specific c-elegans graph timestep
time_step = 560000
for symmetry_breaking in [0.0,0.75,1.0]
        export_graph_nodes_edges(;
            time_step=time_step, B_to_c=1.5, selection_strength=0.2,
            symmetry_breaking=symmetry_breaking, adj_matrix_source="c-elegans",
            time_steps=8000000)
end

# Save cumulative
for type in ["theory","approx","simulation"]
    for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                              "c-elegans-undirected"]
        for symmetry_breaking in [0.0, 0.25, 0.5, 0.75, 1.0]
                extract_cumulative(;
                    type=type, selection_strength=0.2,
                    symmetry_breaking=symmetry_breaking,
                    adj_matrix_source=adj_matrix_source,
                    time_steps=200000000)
        end
    end
end

# Save timeseries_statistics
for adj_matrix_source in ["well-mixed", "c-elegans"]
        for symmetry_breaking in [0.0, 0.25, 0.5, 0.75, 1.0]
            extract_timeseries_statistics(;
                B_to_c=1.5, selection_strength=0.2,
                symmetry_breaking=symmetry_breaking, adj_matrix_source=adj_matrix_source,
                early_cutoff_fraction=0.1, time_steps=800000)
        end
end

# Save mutation timesteps
for adj_matrix_source in ["well-mixed", "c-elegans"]
        for symmetry_breaking in [0.0, 0.25, 0.5, 0.75, 1.0]
            extract_mutation_timesteps(;
                B_to_c=1.5, selection_strength=0.2,
                symmetry_breaking=symmetry_breaking, adj_matrix_source=adj_matrix_source,
                time_steps=800000)
        end
end

# Save chimera indices
covariance_cutoff=1500
for community_algorithm in ["covariance", "infomap"]
    extract_chimera_indices(;
        community_algorithm=community_algorithm,
        B_to_c=1.5, selection_strength=0.2,
        adj_matrix_source="c-elegans",
        time_steps=8000000,
        covariance_cutoff=covariance_cutoff)
end

# Save game types
for adj_matrix_source in ["well-mixed", "c-elegans"]
    extract_game_types(
        B_to_c=1.5, selection_strength=0.2,
        adj_matrix_source=adj_matrix_source,
        time_steps=8000000)
end

# Save results in NetCDF format
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                          "c-elegans-undirected"]
    for (data_type, decimation_factor) in zip(["cumulative", "timeseries-statistics", "timeseries-statistics"], [nothing, nothing, 10000])
        create_netcdf(; data_type=data_type, adj_matrix_source=adj_matrix_source, cumulative_time_steps=2000, timeseries_time_steps=8000, decimation_factor=decimation_factor)
    end
end
