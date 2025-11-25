# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson, Test
using LinearAlgebra
using Random
using Graphs
using SimpleWeightedGraphs
quickactivate("..", "Chimera_EGT_Kuramoto")

include(srcdir("julia", "moran.jl"))

@testset "Deterministic, undirected, prisoner's dilemma, well-mixed graph" begin
    payoff_matrix = [3 1;4 2]
    N = 4
    interaction_graph = ones(Int,N,N) - I
    reproduction_graph = interaction_graph + I
    selection_strength = 100
    mutation_rate = 0
    time_steps = 10

    lmi = Moran(NormalFormGame(payoff_matrix),
      interaction_graph, reproduction_graph,
      selection_strength, mutation_rate)
    actions = [1,1,2,1]
    aux = WorkParams(lmi, actions)
    rng = Xoshiro(1)
    for i = 1:time_steps
      play!(actions,rng,lmi,aux)
    end

    # After enough turns, the single defect (strategy==2)
    # invader should take over
    @test actions == [2,2,2,2]
end

@testset "Directed, prisoner's dilemma, spreader graph" begin
    payoff_matrix = [3 1;4 2]
    N = 5
    # Note: star digraph gives a graph with edges radiating *outward*
    interaction_graph = adjacency_matrix(star_digraph(N))
    reproduction_graph = interaction_graph + I
    # Neutral selection if we want cooperative strategy to win
    selection_strength = 0
    mutation_rate = 0
    time_steps = 100

    lmi = Moran(NormalFormGame(payoff_matrix),
      interaction_graph, reproduction_graph,
      selection_strength, mutation_rate)
    actions = [1,2,2,2,2]
    aux = WorkParams(lmi, actions)

    # Check the payoff (payoffs_colplayer_per_strategy[i,alpha])
    # that player i gives to a row player wish strategy alpha
    # If row player cooperates (strategy==1), player 1 gives reward (3)
    # and everyone else gives sucker (1); if row player defects (strategy==2),
    # then player 1 gives temptation (4) and everyone else give punishment (2)
    @test aux.payoffs_colplayer_per_strategy == [[3,1,1,1,1] [4,2,2,2,2]]
    # Check the payoff (payoffs_rowplayer_per_strategy[i,alpha])
    # that player i receives from a row player wish strategy alpha
    # If row player cooperates (strategy==1), player 1 gets reward (3)
    # and everyone else gets temptation (4); if row player defects (strategy==2),
    # then player 1 gets sucker (1) and everyone else gets punishment (2)
    @test aux.payoffs_rowplayer_per_strategy == [[3,4,4,4,4] [1,2,2,2,2]]
    # Check actual payoffs between each player (respecting interaction
    # graph); the only actual payoffs are from player 1 to the other
    # players (temptation, 4)
    @test transpose(aux.payoffs_player_pairwise_transpose) ==
    [0 0 0 0 0;4 0 0 0 0;4 0 0 0 0;4 0 0 0 0;4 0 0 0 0]
    # Check each player's total payoff; player 1 has no in-neighbors,
    # and every other player gets the temptation (4) reward from the
    # player 1
    calc_payoffs!(aux,actions,lmi)
    @test aux.payoffs == [0,4,4,4,4]

    # After enough turns, the single cooperator (strategy==1)
    # in the center node should spread to all the other nodes
    # because of the directionality
    actions = [1,2,2,2,2]
    aux = WorkParams(lmi, actions)
    rng = Xoshiro(1)
    for i = 1:time_steps
      play!(actions,rng,lmi,aux)
    end
    @test actions == [1,1,1,1,1]
end

@testset "Directed, prisoner's dilemma, random graph" begin
    payoff_matrix = [1 2 3;4 5 6;7 8 9]
    N = 3
    graph = cycle_digraph(N)
    add_edge!(graph, 2, 1)

    # Selection strength will not be used, so set it arbitrarily
    selection_strength = 0
    mutation_rate = 0

    lmi = Moran(NormalFormGame(payoff_matrix),
      adjacency_matrix(graph), adjacency_matrix(graph),
      selection_strength, mutation_rate)
    actions = [1,2,3]
    aux = WorkParams(lmi, actions)

    # Check initial payoffs
    @test transpose(aux.payoffs_player_pairwise_transpose) ==
    [0 2 3;4 0 0;0 8 0]
    @test aux.payoffs_colplayer_per_strategy == [1 4 7;2 5 8;3 6 9]
    @test aux.payoffs_rowplayer_per_strategy == [1 2 3;4 5 6;7 8 9]

    # Update node 2 to strategy 3
    update_aux!(aux, 3, 2, lmi)
    @test transpose(aux.payoffs_player_pairwise_transpose) ==
    [0 3 3;7 0 0;0 9 0]
    @test aux.payoffs_colplayer_per_strategy == [1 4 7;3 6 9;3 6 9]
    @test aux.payoffs_rowplayer_per_strategy == [1 2 3;7 8 9;7 8 9]
end

@testset "Weighted, undirected prisoner's dilemma" begin
  payoff_matrix = [1 2 3;4 5 6;7 8 9]
    N = 3
    graph = SimpleWeightedGraph{Int64,Int64}(N)
    add_edge!(graph, 1, 2, 10)
    add_edge!(graph, 2, 3, 1)
    add_edge!(graph, 1, 3, 7)

    # Selection strength will not be used, so set it arbitrarily
    selection_strength = 0
    mutation_rate = 0

    lmi = Moran(NormalFormGame(payoff_matrix),
      adjacency_matrix(graph), adjacency_matrix(graph),
      selection_strength, mutation_rate)
    actions = [1,2,3]
    aux = WorkParams(lmi, actions)

    # Check initial payoffs
    @test transpose(aux.payoffs_player_pairwise_transpose) ==
    [0 20 21;40 0 6;49 8 0]

    # Update node 1 to strategy 2
    update_aux!(aux, 2, 1, lmi)
    @test transpose(aux.payoffs_player_pairwise_transpose) ==
    [0 50 42;50 0 6;56 8 0]
end
