using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Here you may include files from the source directory
include(srcdir("julia", "local_moran_interaction.jl"))

# Save specific graph timestep
time_step = 560000
local_moran_interaction.export_graph_nodes_edges(
    time_step;B_factor=1.5, selection_strength=0.2,
    symmetry_breaking=0.75, adj_matrix_source="c-elegans",
    time_steps=8000000)

# Save cumulative
for type in ["theory","approx","simulation"]
    for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                              "c-elegans-undirected"]
        for symmetry_breaking in [0.0, 0.25, 0.5, 0.75, 1.0]
                local_moran_interaction.extract_cumulative(
                    type; selection_strength=0.2,
                    symmetry_breaking=symmetry_breaking,
                    adj_matrix_source=adj_matrix_source,
                    time_steps=200000000)
        end
    end
end

# Save timeseries_statistics
for adj_matrix_source in ["well-mixed", "c-elegans"]
    local_moran_interaction.extract_timeseries_statistics(
        B_factor=1.5, selection_strength=0.2,
        symmetry_breaking=0.75, adj_matrix_source=adj_matrix_source,
        time_steps=800000)
end

# Save chimera indices
covariance_cutoff=1500
for community_algorithm in ["covariance", "infomap"]
    local_moran_interaction.extract_chimera_indices(
        community_algorithm;
        B_factor=1.5, selection_strength=0.2,
        adj_matrix_source="c-elegans",
        time_steps=8000000,
        covariance_cutoff=covariance_cutoff)
end

# Save game types
for adj_matrix_source in ["well-mixed", "c-elegans"]
    local_moran_interaction.extract_game_types(
        B_factor=1.5, selection_strength=0.2,
        adj_matrix_source=adj_matrix_source,
        time_steps=8000000)
end
