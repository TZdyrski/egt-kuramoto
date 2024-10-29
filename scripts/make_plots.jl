using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Here you may include files from the source directory
include(srcdir("julia", "local_moran_interaction.jl"))

update_method="single-update"

# Plot cumulative
for num_time_steps in Int.([2E3,2E6,2E8])
    for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
			      "c-elegans-undirected", "c-elegans-undirected-unweighted"]
	for selection_strength in [0.005, 0.2]
	    local_moran_interaction.plot_cumulative(selection_strength, adj_matrix_source, update_method, num_time_steps)
	end
    end
end

# Plot time-series
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
			      "c-elegans-undirected", "c-elegans-undirected-unweighted"]
    for (B_factor, selection_strength) in zip([1.5,2.5],[0.2,5])
	for num_time_steps in Int.([8E4, 8E5, 8E6])
            local_moran_interaction.plot_timeseries(B_factor, selection_strength,
						    adj_matrix_source, update_method, num_time_steps)
	end
    end
end

# Plot heatmap
local_moran_interaction.save_heatmap()

# Plot game-type regions
local_moran_interaction.plot_payoff_regions()

# Plot colored graph
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
			      "c-elegans-undirected", "c-elegans-undirected-unweighted"]
    for (B_factor, selection_strength) in zip([1.5,2.5],[0.2,5])
	for num_time_steps in Int.([8E4, 8E5, 8E6])
            local_moran_interaction.plot_graph_evolution(B_factor, selection_strength,
						    adj_matrix_source, update_method, num_time_steps)
	end
    end
end
