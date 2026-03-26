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
      play!(rng,lmi,aux)
    end

    # After enough turns, the single defect (strategy==2)
    # invader should take over
    @test aux.interaction_graph.ndata.strategy == [2,2,2,2]
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

    # Check each edge's payoff; every edge connects player 1 (cooperate)
    # to a defector, giving each edge the temptation (4) reward
    @test aux.interaction_graph.edata.payoffs == [4,4,4,4]
    # Check each player's total payoff; player 1 has no in-neighbors,
    # and every other player gets the temptation (4) reward from the
    # player 1
    @test aux.interaction_graph.ndata.payoffs == [0,4,4,4,4]

    # After enough turns, the single cooperator (strategy==1)
    # in the center node should spread to all the other nodes
    # because of the directionality
    actions = [1,2,2,2,2]
    aux = WorkParams(lmi, actions)
    rng = Xoshiro(1)
    for i = 1:time_steps
      play!(rng,lmi,aux)
    end
    @test aux.interaction_graph.ndata.strategy == [1,1,1,1,1]
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

    # Check edge order
    @test [src(e) => dst(e) for e in edges(aux.interaction_graph)] ==
	    [2 => 1, 3 => 1, 1 => 2, 2 => 3]

    # Check initial payoffs
    @test aux.interaction_graph.edata.payoffs == [2,3,4,8]

    # Update node 2 to strategy 3
    aux.interaction_graph.ndata.strategy[2] = 3
    update_aux!(aux, 2, lmi)
    @test aux.interaction_graph.edata.payoffs == [3,3,7,9]
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

    # Check edge order
    @test [src(e) => dst(e) for e in edges(aux.interaction_graph)] ==
	    [2 => 1, 3=> 1, 1 => 2, 3 => 2, 1 => 3, 2 => 3]

    # Check initial payoffs
    @test aux.interaction_graph.edata.payoffs == [20,21,40,6,49,8]

    # Update node 1 to strategy 2
    aux.interaction_graph.ndata.strategy[1] = 2
    update_aux!(aux, 1, lmi)
    @test aux.interaction_graph.edata.payoffs == [50,42,50,6,56,8]
end
