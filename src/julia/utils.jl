# SPDX-License-Identifier: GPL-3.0-or-later

using DataToolkit
using LinearAlgebra
using BlockArrays
using Random
using Graphs
using SimpleWeightedGraphs
using XLSX

"Add together the elements from the first half of a list with an even number of elements."
function extract_num_communicative(players_per_strategy::AbstractVector{<:Integer})
    @assert iseven(length(players_per_strategy))
    # Define matrix that extracts number of communicative players
    N = div(length(players_per_strategy), 2)
    communicative_matrix = [ones(Int, N); zeros(Int, N)]
    # Communicative have value < nb_phases
    num_communicative = dot(players_per_strategy, communicative_matrix)
    return num_communicative
end

"""
    payoff_matrix(nb_phases, mutual_benefit_synchronous, unilateral_benefit_synchronous,
      cost; symmetry_breaking=1/2)

Create a square payoff matrix for the symmetric Evolutionary Kuramoto game.

The returned matrix has size `(2*nb_phases, 2*nb_phases)`
to represent `nb_phases` for both the communicative and non-communicative strategies.
Note: The cost is added to each cell to prevent negative entries.
"""
function payoff_matrix(nb_phases::Integer,
                       mutual_benefit_synchronous::Real,
                       unilateral_benefit_synchronous::Real, cost::Real;
                       symmetry_breaking::Real=1 / 2)
    benefit_scaling = [(1 + cospi(2 * (phi_j - phi_i))) / 2
                       for phi_i in (0:(nb_phases - 1)) / nb_phases,
                           phi_j in (0:(nb_phases - 1)) / nb_phases]
    payoff_matrix = Matrix(mortar([[mutual_benefit_synchronous * benefit_scaling .- cost,
                                    2 * (1 - symmetry_breaking) *
                                    unilateral_benefit_synchronous * benefit_scaling];;
                                   [2 * symmetry_breaking * unilateral_benefit_synchronous *
                                    benefit_scaling .- cost,
                                    zeros(eltype(benefit_scaling), size(benefit_scaling))]]) .+
                           cost)
    return payoff_matrix
end

function add_self_loops_to_dangling_nodes!(graph::AbstractGraph)
	dangling_node_idxs = findall(outdegree(graph) .== 0)

	# Add self-loops
	for node_idx in dangling_node_idxs
		add_edge!(graph, node_idx, node_idx)
	end

	return graph
end

# Note: C-elegans removes the 2 disconnected neurons
function get_adj_matrices(; adj_matrix_source::String, nb_players::Integer = 20, regular_degree::Integer = 10, rng::AbstractRNG=Xoshiro(1)) # Define interaction graph without loops
    # Define reproduction graph with loops
    if adj_matrix_source == "well-mixed"
        interaction_adj_matrix = ones(Int64, nb_players, nb_players) - I
    elseif adj_matrix_source == "random-regular-graph"
        interaction_adj_matrix = adjacency_matrix(random_regular_graph(nb_players, regular_degree; rng))
    elseif adj_matrix_source == "random-regular-digraph"
        interaction_adj_matrix = adjacency_matrix(random_regular_digraph(nb_players, regular_degree; rng))
    elseif adj_matrix_source == "c-elegans-unweighted"
        interaction_adj_matrix = collect(round.(get_celegans_connectome()) .!= 0)
    elseif adj_matrix_source == "c-elegans-undirected"
        interaction_adj_matrix = round.(get_celegans_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
    elseif adj_matrix_source == "c-elegans-undirected-unweighted"
        interaction_adj_matrix = round.(get_celegans_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
        interaction_adj_matrix = collect(interaction_adj_matrix .!= 0)
    elseif adj_matrix_source == "c-elegans"
        interaction_adj_matrix = round.(get_celegans_connectome())
    else
        throw(ArgumentError("adj_matrix_source must be a string in set [\"well-mixed\", "
                            * "\"c-elegans\", \"c-elegans-unweighted\", "
                            * "\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\", "
			    * "\"random-regular-graph\", \"random-regular-digraph\"]"))
    end
		reproduction_adj_matrix = adjacency_matrix(add_self_loops_to_dangling_nodes!(SimpleWeightedDiGraph(interaction_adj_matrix)))
    return interaction_adj_matrix, reproduction_adj_matrix
end

function get_celegans_connectome()
    connectome = get_celegans_connectome_labelled()["connectome"]
    # Replace "Missing" data with zeros
    replace!(connectome, missing => 0)
    return connectome
end

function get_celegans_connectome_labelled()
    connectome_and_muscles_with_labels = read(dataset("celegans-connectome-cook"), Matrix)
    connectome_with_labels = connectome_and_muscles_with_labels[:, [1:21..., 52:323..., 439:446...]]
    row_labels = connectome_with_labels[2:end, 1]
    col_labels = connectome_with_labels[1, 2:end]
    if row_labels != col_labels
      throw(ErrorException("Rows and columns in c-elegans adjacency table do not match up"))
    end
    # Note: Cook uses the convention row => src node
    connectome = connectome_with_labels[2:end, 2:end]
    results = Dict("connectome" => connectome, "row_labels" => row_labels,
                   "col_labels" => col_labels)
    return results
end

function encode_delta_encoded(data::AbstractMatrix)
	starting_data = data[:,1]
	deltas = diff(data; dims=2)
	return starting_data, deltas
end

function decode_delta_encoded(starting_data::AbstractVector,deltas::AbstractMatrix,time_step::Integer)
  # Note: starting_data corresponds to time_step = 0
  # starting_data should be a n-vector
  # and deltas should be an nxm matrix of differences
  # between successive time steps
  state = starting_data
  for i = 1:time_step
	  state .+= deltas[:,i]
  end
  return state
end

function decode_delta_encoded_all(starting_data::AbstractVector,deltas::AbstractMatrix,final_time_step::Integer=size(deltas,2))
  trimmed_deltas = deltas[:,1:final_time_step]
  state = Matrix{eltype(trimmed_deltas)}(undef, length(starting_data), final_time_step+1)
  cumsum!(state, hcat(starting_data, trimmed_deltas); dims=2)
  return state
end
