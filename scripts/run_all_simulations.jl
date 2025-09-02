using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Here you may include files from the source directory
include(srcdir("julia", "run_simulations.jl"))

update_method = "single-update"

# Calc cumulative
num_time_steps_cumulative = Int(2E8)
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
		      "c-elegans-undirected", "c-elegans-undirected-unweighted",
		      "random-regular-graph", "random-regular-digraph"]
    for selection_strength in [0.005, 0.2, 5]
        for symmetry_breaking in [0, 1 / 4, 1 / 2, 3 / 4, 1]
		calc_and_save_cumulative(selection_strength,
							symmetry_breaking,
							adj_matrix_source, update_method,
							num_time_steps_cumulative)
	    end
	end
end

# Calc time-series
num_time_steps_timeseries = Int(8E6)
for adj_matrix_source in ["well-mixed", "c-elegans", "c-elegans-unweighted",
                          "c-elegans-undirected", "c-elegans-undirected-unweighted",
                          "random-regular-graph", "random-regular-digraph"]
    for (B_to_c, selection_strength) in Iterators.product([1.5, 2.5], [0.005, 0.2, 5])
        for symmetry_breaking in [0, 1 / 4, 1 / 2, 3 / 4, 1]
            calc_and_save_timeseries(B_to_c, selection_strength,
                                                    symmetry_breaking,
                                                    adj_matrix_source, update_method,
                                                    num_time_steps_timeseries)
        end
    end
end
