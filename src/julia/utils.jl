# SPDX-License-Identifier: GPL-3.0-or-later

using DataToolkit
using LinearAlgebra
using BlockArrays
using Random

function extract_num_communicative(players_per_strategy::AbstractVector{<:Integer})
    @assert iseven(length(players_per_strategy))
    # Define matrix that extracts number of communicative players
    N = div(length(players_per_strategy), 2)
    communicative_matrix = [ones(Int, N); zeros(Int, N)]
    # communicatative have value < nb_phases
    num_communicative = dot(players_per_strategy, communicative_matrix)
    return num_communicative
end

function payoff_matrix(nb_phases::Integer,
                       mutual_benefit_synchronous::Real,
                       unilateral_benefit_synchronous::Real, cost::Real;
                       symmetry_breaking::Real=1 / 2)
    benefit_scaling = [(1 + cos(2 * pi * (phi_j - phi_i))) / 2
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

function get_adj_matrices(; adj_matrix_source::String, nb_players::Integer = 20, regular_degree::Integer = 10, rng::AbstractRNG=Xoshiro(1))
    # Define interaction graph without loops
    # Define reproduction graph with loops
    if adj_matrix_source == "well-mixed"
        interaction_adj_matrix = ones(Int64, nb_players, nb_players) - I
        reproduction_adj_matrix = interaction_adj_matrix + I
    elseif adj_matrix_source == "random-regular-graph"
        interaction_adj_matrix = adjacency_matrix(random_regular_graph(nb_players, regular_degree; rng))
        reproduction_adj_matrix = interaction_adj_matrix + I
    elseif adj_matrix_source == "random-regular-digraph"
        interaction_adj_matrix = adjacency_matrix(random_regular_digraph(nb_players, regular_degree; rng))
        reproduction_adj_matrix = interaction_adj_matrix + I
    elseif adj_matrix_source == "c-elegans-unweighted"
        interaction_adj_matrix = collect(round.(get_celegans_connectome()) .!= 0)
        reproduction_adj_matrix = collect((interaction_adj_matrix + I) .!= 0)
    elseif adj_matrix_source == "c-elegans-undirected"
        interaction_adj_matrix = round.(get_celegans_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
        reproduction_adj_matrix = interaction_adj_matrix + I
    elseif adj_matrix_source == "c-elegans-undirected-unweighted"
        interaction_adj_matrix = round.(get_celegans_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
        interaction_adj_matrix = collect(interaction_adj_matrix .!= 0)
        reproduction_adj_matrix = collect((interaction_adj_matrix + I) .!= 0)
    elseif adj_matrix_source == "c-elegans"
        interaction_adj_matrix = round.(get_celegans_connectome())
        reproduction_adj_matrix = interaction_adj_matrix + I
    else
        throw(ArgumentError("adj_matrix_source must be a string in set [\"well-mixed\", "
                            * "\"c-elegans\", \"c-elegans-unweighted\", "
                            * "\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\", "
			    * "\"random-regular-graph\", \"random-regular-digraph\"]"))
    end
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
    connectome = connectome_with_labels[2:end, 2:end]
    results = Dict("connectome" => connectome, "row_labels" => row_labels,
                   "col_labels" => col_labels)
    return results
end
