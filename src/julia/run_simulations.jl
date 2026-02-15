# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson

include("moran.jl")
include("utils.jl")

function calc_cumulative(;selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         time_steps::Integer=2_000_000,
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 beta_to_B::Real=0.95,
			 mutation_rate::Real=0.0001,
			 nb_players::Integer=20,
			 seed::Integer=12345,
			 fixed_aspect::Union{String,Nothing}=nothing,
			 )

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(; adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Define Bs on which to run
    B_crit = 2 * cost * (nb_players - 1) / (nb_players - 2)
    nb_Bs = 11
    step_size_Bs = 0.04
    Bs = B_crit .+ ((0:(nb_Bs - 1)) .- (nb_Bs - 1) / 2) .* step_size_Bs

		# Generate model
		nfgames = [NormalFormGame(payoff_matrix(nb_phases, B, beta_to_B * B, cost; symmetry_breaking)) for B in Bs]
		models = [Moran(nfgame, interaction_adj_matrix, reproduction_adj_matrix, selection_strength, mutation_rate, fixed_aspect) for nfgame in nfgames]

    # Run the model for weak selection strength
		cumulative_populations = [cumulative(model, time_steps, seed) for model in models]

    # Plot fraction communicative
    nb_communicative = [extract_num_communicative(final_population)
                        for final_population in cumulative_populations]
    fraction_communicative = nb_communicative ./ (time_steps * nb_players)

    # Package results
    return @strdict(Bs, nb_players, fraction_communicative)
end

function calc_timeseries(;B_to_c::Real, selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         time_steps::Integer=80_000,
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 beta_to_B::Real=0.95,
			 mutation_rate::Real=0.0001,
			 nb_players::Integer=20,
			 seed::Integer=12345,
			 fixed_aspect::Union{String,Nothing}=nothing,
			 )

    # Define system
    B = cost * B_to_c

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(; adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

		# Generate model
		nfgame = NormalFormGame(payoff_matrix(nb_phases, B, beta_to_B * B, cost; symmetry_breaking))
		model = Moran(nfgame, interaction_adj_matrix, reproduction_adj_matrix, selection_strength, mutation_rate, fixed_aspect)

    # Run the model for weak selection strength
    initial_actions, deltas, steps_following_mutation = time_series(model, time_steps, seed)

    # Package results
    return @strdict(initial_actions, deltas, steps_following_mutation, nb_phases, nb_players, interaction_adj_matrix)
end
