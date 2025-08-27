using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Here you may include files from the source directory
include(srcdir("julia", "run_simulations.jl"))
include(srcdir("julia", "plotting.jl"))

update_method = "single-update"

# Plot cumulative
for num_time_steps in Int.([2E3, 2E6, 2E8])
    for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                              "c-elegans-undirected", "c-elegans-undirected-unweighted",
                              "random-regular-graph", "random-regular-digraph"]
        for selection_strength in [0.005, 0.2, 5]
            for symmetry_breaking in [0, 1 / 4, 1 / 2, 3 / 4, 1]
                plot_cumulative(selection_strength,
                                                        symmetry_breaking,
                                                        adj_matrix_source, update_method,
                                                        num_time_steps)
            end
        end
    end
end

# Plot time-series
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                          "c-elegans-undirected", "c-elegans-undirected-unweighted",
                          "random-regular-graph", "random-regular-digraph"]
    for (B_to_c, selection_strength) in Iterators.product([0.005, 1.5, 2.5], [0.2, 5])
        for num_time_steps in Int.([8E4, 8E5, 8E6])
            for symmetry_breaking in [0, 1 / 4, 1 / 2, 3 / 4, 1]
                plot_timeseries(B_to_c, selection_strength,
                                                        symmetry_breaking,
                                                        adj_matrix_source, update_method,
                                                        num_time_steps)
            end
        end
    end
end

# Plot heatmap
save_heatmap()

# Plot game-type regions
plot_payoff_regions()

# Plot colored graph
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                          "c-elegans-undirected", "c-elegans-undirected-unweighted",
                          "random-regular-graph", "random-regular-digraph"]
    for (B_to_c, selection_strength) in Iterators.product([1.5, 2.5], [0.2, 5])
        for num_time_steps in Int.([8E4, 8E5, 8E6])
            for symmetry_breaking in [0, 1 / 4, 1 / 2, 3 / 4, 1]
                plot_graph_evolution(B_to_c, selection_strength,
                                                             symmetry_breaking,
                                                             adj_matrix_source,
                                                             update_method, num_time_steps)
            end
        end
    end
end
