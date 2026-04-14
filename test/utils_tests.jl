# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson, Test
using LinearAlgebra
using Graphs
quickactivate("..", "Chimera_EGT_Kuramoto")

include(srcdir("julia", "utils.jl"))

@testset "Extract number communicative" begin
    @test extract_num_communicative([5,0]) == 5
    @test extract_num_communicative([0,3]) == 0
    @test extract_num_communicative([2,6]) == 2
    @test extract_num_communicative([2,8,6,5]) == 10
    @test_throws AssertionError extract_num_communicative([2,6,4])
end

@testset "Payoff matrix" begin
    @test size(payoff_matrix(3,0,0,0)) == (6,6)
    @test payoff_matrix(1,0,0,0) == zeros(2,2)
    @test payoff_matrix(1,3,0,0) == [3 0; 0 0]
    @test payoff_matrix(1,4,2,0) == [4 2; 2 0]
    @test payoff_matrix(1,5,2,1) == [5 2; 3 1]
    @test payoff_matrix(1,5,4,2;symmetry_breaking=0.25) == [5 2; 8 2]
    @test payoff_matrix(2,6,2,1) == [6 0 2 0; 0 6 0 2; 3 1 1 1; 1 3 1 1]
end

@testset "Adjacency matrices" begin
    # Test well-mixed gives all-ones except diagonal (self-loops) for interaction matrix
    # and reproduction matrix
    @test get_adj_matrices(adj_matrix_source="well-mixed") == (ones(20,20) - I, ones(20,20) - I)
    # Test well-mixed has size (nb_players, nb_players)
    @test get_adj_matrices(adj_matrix_source="well-mixed",nb_players=5) == (ones(5,5) - I, ones(5,5) - I)
    # Test that c-elegans has size (300,300) since we remove the 2
    # unconnected neurons
    @test size(get_adj_matrices(adj_matrix_source="c-elegans")[1]) == (300,300)
		# Test that the c-elegans interaction graph does have dangling nodes
		@test any(outdegree(SimpleDiGraph(get_adj_matrices(adj_matrix_source="c-elegans")[1])) .== 0)
		# Test that the c-elegans reproduction graph has no dangling nodes
		@test !any(outdegree(SimpleDiGraph(get_adj_matrices(adj_matrix_source="c-elegans")[2])) .== 0)
    # Test that c-elegans undirected matrices are symmetric
    @test issymmetric(get_adj_matrices(adj_matrix_source="c-elegans-undirected")[1])
    @test issymmetric(get_adj_matrices(adj_matrix_source="c-elegans-undirected-unweighted")[1])
    # Test that c-elegans unweighted matrices only have elements 0 or 1
    @test all(x in [0,1] for x in get_adj_matrices(adj_matrix_source="c-elegans-unweighted")[1])
    @test all(x in [0,1] for x in get_adj_matrices(adj_matrix_source="c-elegans-undirected-unweighted")[1])
    # Test that random regular graph is symmetric
    @test issymmetric(get_adj_matrices(adj_matrix_source="random-regular-graph")[1])
    # Test that random regular digraph has size (20,20)
    @test size(get_adj_matrices(adj_matrix_source="random-regular-digraph")[1]) == (20,20)
end

@testset "C-elegans dataset" begin
    connectome = get_celegans_connectome_labelled()
    # Compare to published data set
    # I2L -> I4 (good check, because I4 -> I2l does *not* exist)
    @test connectome["row_labels"][3] == "I2L"
    @test connectome["col_labels"][6] == "I4"
    @test connectome["connectome"][3,6] == 13
    # VC04 -> VC03 (good check, because it's near the end and VC03 -> VC04 does *not* exist)
    @test connectome["row_labels"][298] == "VC04"
    @test connectome["col_labels"][297] == "VC03"
    @test connectome["connectome"][298,297] == 1
    @test get_celegans_connectome()[1:3,1:6] == [0 0 10 0 3 0;0 0 0 6 1 0;2 0 0 3 0 13]
end
