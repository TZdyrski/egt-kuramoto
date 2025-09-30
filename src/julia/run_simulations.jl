# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson

include("moran.jl")
include("utils.jl")

function calc_and_save_cumulative(;selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         time_steps::Integer=2_000_000,
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 beta_to_B::Real=0.95,
			 mutation_rate::Real=0.0001,
			 nb_players::Integer=20,
			 )

    # Load results
    config = @strdict(adj_matrix_source, time_steps,
                      selection_strength, symmetry_breaking, nb_phases, cost, beta_to_B, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end
    data, _ = produce_or_load(calc_cumulative, config, datadir("raw","cumulative"))
end

function calc_cumulative(config::Dict)
    # Unpack values
    @unpack selection_strength, symmetry_breaking, adj_matrix_source, time_steps, nb_phases, cost, beta_to_B, mutation_rate = config

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(; adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Define Bs on which to run
    B_crit = 2 * cost * (nb_players - 1) / (nb_players - 2)
    nb_Bs = 11
    step_size_Bs = 0.04
    Bs = B_crit .+ ((0:(nb_Bs - 1)) .- (nb_Bs - 1) / 2) .* step_size_Bs

    # Run the model for weak selection strength
    @time cumulative_populations = [cumulative(Moran(NormalFormGame(payoff_matrix(nb_phases,
                                                                                                  B,
                                                                                                  beta_to_B *
                                                                                                  B,
                                                                                                  cost;
                                                                                                  symmetry_breaking)),
                                                                     interaction_adj_matrix,
                                                                     reproduction_adj_matrix,
                                                                     selection_strength,
                                                                     mutation_rate),
                                               time_steps)
                                    for B in Bs]

    # Plot fraction communcative
    nb_communicative = [extract_num_communicative(final_population)
                        for final_population in cumulative_populations]
    fraction_communicative = nb_communicative ./ (time_steps * nb_players)

    # Package results
    return @strdict(Bs, nb_players, fraction_communicative)
end

function calc_and_save_timeseries(;B_to_c::Real, selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         time_steps::Integer=80_000,
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 beta_to_B::Real=0.95,
			 mutation_rate::Real=0.0001,
			 nb_players::Integer=20,
			 )
    # Load results
    config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                      symmetry_breaking, selection_strength, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end
    data, _ = produce_or_load(calc_timeseries, config, datadir("raw","timeseries"))
end


function calc_timeseries(config::Dict)
    # Unpack variables
    @unpack B_to_c, selection_strength, symmetry_breaking, adj_matrix_source, time_steps, nb_phases, cost, beta_to_B, mutation_rate = config

    # Define system
    B = cost * B_to_c

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(; adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Run the model for weak selection strength
    all_populations, steps_following_mutation = time_series(Moran(NormalFormGame(payoff_matrix(nb_phases,
                                                                                     B,
                                                                                     beta_to_B *
                                                                                     B,
                                                                                     cost;
                                                                                     symmetry_breaking)),
                                                        interaction_adj_matrix,
                                                        reproduction_adj_matrix,
                                                        selection_strength, mutation_rate),
                                  time_steps)

    # Package results
    return @strdict(all_populations, steps_following_mutation, nb_phases, nb_players, interaction_adj_matrix)
end
