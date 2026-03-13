# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson, Test
using DimensionalData
quickactivate("..", "Chimera_EGT_Kuramoto")

include(srcdir("julia", "postprocess.jl"))

@testset "Order parameter" begin
	# All players have the same strategy
	@test extract_order_parameters([0,5,0,0],2) == 1
	# All players have different strategies
	@test extract_order_parameters([1,1,1,1],2) ≈ 0 atol=1e-16
	# Players split into 2 equally sized groups at 90 degrees
	@test extract_order_parameters([0,3,3,0,0,0,0,0],4) ≈ abs(1+1im)/2
	# Communicative half and non-communicative half correctly
	# identified with each other
	@test extract_order_parameters([0,1,0,2],2) == 1
end

@testset "Same strategy check" begin
	nb_phases = 2
	@test check_all_same_strategy([1,1,2,1,2,1],nb_phases) == allCommunicative
	@test check_all_same_strategy([3,3,3,3,4,3],nb_phases) == allNoncommunicative
	@test isnothing(check_all_same_strategy([1,2,3,3,4,3],nb_phases))
end

@testset "Disconnected synchronized check" begin
	@test check_disconnected_synchronized(Dict{Union{Missing,GameType},Integer}(missing => 4)) == disconnectedSynchronizedPopulations
	@test isnothing(check_disconnected_synchronized(Dict{Union{Missing,GameType},Integer}(missing => 4, chicken => 1)))
end

@testset "Game types" begin
	B0 = 4
	beta0 = 2
	cost = 1
	alpha = 1/2
	nb_phases = 2

	game_types = game_types_per_strategy_pair(B0, beta0, cost, alpha, nb_phases)
	# Payoff matrix:
	#        C(1) C(2) N(3) N(4)
	# C(1) [ 4    0    2    0 ]
	# C(2) [ 0    4    0    2 ]
	# N(3) [ 3    1    1    1 ]
	# N(4) [ 1    3    1    1 ]
	@test payoff_matrix(nb_phases, B0, beta0, cost; symmetry_breaking=alpha) ==
	  [4 0 2 0;0 4 0 2;3 1 1 1;1 3 1 1]

	nb_players = 6
	interaction_adj_matrix = ones(Int64,nb_players,nb_players) - I

	# Of 15 (double-sided) edges, 10 are neutral
	expected_games = [missing neutral neutral neutral neutral concord;
									  neutral missing neutral neutral neutral concord;
									  neutral neutral missing neutral neutral concord;
									  neutral neutral neutral missing neutral concord;
									  neutral neutral neutral neutral missing concord;
										concord concord concord concord concord missing];
	games, game_counts = count_games([3,3,3,3,3,1];game_types,interaction_adj_matrix)
	@test isequal(games, expected_games)
	@test extract_most_common_game_types(game_counts) == neutral

	# All 5 mixed-type games are concord
	game_types_only_mixed = game_types_per_strategy_pair(B0, beta0, cost, alpha, nb_phases, only_mixed_games=true)
	expected_games = [missing missing missing missing missing concord;
									  missing missing missing missing missing concord;
									  missing missing missing missing missing concord;
									  missing missing missing missing missing concord;
									  missing missing missing missing missing concord;
										concord concord concord concord concord missing];
	games, game_counts = count_games([3,3,3,3,3,1];game_types=game_types_only_mixed,interaction_adj_matrix)
	@test isequal(games, expected_games)
	@test extract_most_common_game_types(game_counts) == concord

	# Of 15 (double-sided) edges, 9 are 1-3 interactions:
	#        C(1) N(3)
	# C(1) [ 4    2 ]
	# N(3) [ 3    1 ] (left-up convention)
	# which is a concord-type game
	expected_games = [missing neutral neutral concord concord concord;
									  neutral missing neutral concord concord concord;
									  neutral neutral missing concord concord concord;
									  concord concord concord missing neutral neutral;
									  concord concord concord neutral missing neutral;
										concord concord concord neutral neutral missing];
	games, game_counts = count_games([3,3,3,1,1,1];game_types,interaction_adj_matrix)
	@test isequal(games, expected_games)
	@test extract_most_common_game_types(game_counts) == concord

	# Ensure that repeated evaluations (with optimizations) gives the same result as one repetition
        initial_strategies = [3,2,3,1,4,1]
        games, game_counts = count_games(initial_strategies;game_types,interaction_adj_matrix)
        new_strategies = [3,2,3,2,4,3]
        changes = new_strategies .!= initial_strategies
        games_partial, game_counts_partial = count_games(new_strategies,games,game_counts,changes;game_types,interaction_adj_matrix)
				game_type_partial = extract_most_common_game_types(game_counts)

	# Compare to just calculating the new game type
        games_full, game_counts_full = count_games(new_strategies;game_types,interaction_adj_matrix)
				game_type_full = extract_most_common_game_types(game_counts_full)

	@test (games_partial, game_counts_partial, game_type_partial) == (games_full, game_counts_full, game_type_full)
end

@testset "Extract counts" begin
	@test extract_counts([1,1,1,1],2) == [4,0,0,0]
	@test extract_counts([2,2,2,2],2) == [0,4,0,0]
	@test extract_counts([4,4,4,4],2) == [0,0,0,4]
	@test extract_counts([1,2,1,3],2) == [2,1,1,0]
end

@testset "Metastability index" begin
	communities = ones(Int64,3)
	nb_phases = 3

	# Always synchronized
	data = DimArray([1 1 1;1 1 1;2 2 2;1 1 1;2 2 2;2 2 2], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["metastability_index"] ≈ 0 atol=1e-16

	# Always disordered
	data = DimArray([1 2 3;2 3 1;3 2 1;3 1 2;2 1 3;1 2 3], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["metastability_index"] == 0

	# Equal times synchronized and disordered (order parameter = 1 or 0)
	time_steps = 6
	corrected = time_steps/(time_steps-1)
	data = DimArray([1 2 3;1 1 1;2 2 2;3 1 2;1 1 1;1 2 3], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["metastability_index"] ≈ 1/4*corrected
end

@testset "Chimera index" begin
	communities = [1,1,1,2,2,2]
	nb_phases = 3

	# Both communities fully synchronized
	data = DimArray([1 1 1 1 1 1;2 2 2 1 1 1], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["chimera_index"] ≈ 0 atol=1e-16

	# Both communities fully disordered
	data = DimArray([1 2 3 3 2 1;3 1 2 1 2 3], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["chimera_index"] == 0

	# Half of the communities synchronized and half disoredered
	data = DimArray([1 2 3 3 3 3;2 2 2 3 1 2], (:time_step, :player_index))
	starting_data, deltas = encode_delta_encoded(transpose(data))
	@test get_chimera_indices(order_parameters_by_community(starting_data, deltas, communities, nb_phases))["chimera_index"] == 1/2
end
