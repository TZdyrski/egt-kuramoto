# SPDX-License-Identifier: GPL-3.0-or-later

using GameTheory
using Random
using StatsBase
using Graphs
using GNNGraphs
using DataFrames

struct Moran{N,T1<:Real,S<:Integer,S2<:Integer,T2<:Real,T3<:Real,
                             AT1<:AbstractMatrix{<:T1},AT2<:AbstractMatrix{<:T1},
                             AS1<:AbstractMatrix{S2},AS2<:AbstractMatrix{S2},
                             AS3<:AbstractMatrix{S2},AS4<:AbstractMatrix{S2}}
    players::NTuple{N,Player{2,T1}}
    num_actions::S
    payoff_matrix::AT1
    payoff_matrix_transpose::AT2
    interaction_adj_matrix::AS1 # Convention: a_ij is number of edges going from i to j
    interaction_adj_matrix_transpose::AS2
    reproduction_adj_matrix::AS3
    reproduction_adj_matrix_transpose::AS4
    beta::T2
    epsilon::T3
end

function Moran(g::NormalFormGame,
                               interaction_adj_matrix::AbstractMatrix,
                               reproduction_adj_matrix::AbstractMatrix,
                               selection_strength::Real,
                               mutation_rate::Real)
    if size(interaction_adj_matrix, 1) != size(interaction_adj_matrix, 2)
        throw(ArgumentError("Interaction adjacency matrix must be square"))
    end
    if size(reproduction_adj_matrix, 1) != size(reproduction_adj_matrix, 2)
        throw(ArgumentError("Reproduction adjacency matrix must be square"))
    end
    if size(interaction_adj_matrix, 1) != size(reproduction_adj_matrix, 1)
        throw(ArgumentError("Interaction adjacency matrix and reproduction adjacency matrix must be same size"))
    end
    if minimum(outdegree(SimpleDiGraph(reproduction_adj_matrix))) == 0
        throw(ArgumentError("Each node of the reproduction graph must have out-degree greater than zero"))
    end
    min_payoff_by_player = [minimum(g.players[i].payoff_array) for i in 1:length(g.players)]
    min_payoff = minimum(min_payoff_by_player)
    if min_payoff < 0
        throw(ArgumentError("Payoff matrix must be non-negative"))
    end
    N = size(interaction_adj_matrix, 1)
    players = ntuple(i -> g.players[1], N)
    num_actions = g.nums_actions[1]
    if num_actions != g.nums_actions[2]
        throw(ArgumentError("Payoff matrix must be square"))
    end
    return Moran(players, num_actions, g.players[1].payoff_array,
                                 transpose(g.players[1].payoff_array),
                                 interaction_adj_matrix, transpose(interaction_adj_matrix),
                                 reproduction_adj_matrix,
                                 transpose(reproduction_adj_matrix), selection_strength,
                                 mutation_rate)
end

function update_weights!(w::Weights{S,TA,V}, new_wts::V) where {S,TA,V}
    w.values = new_wts
    return w.sum = sum(w.values)
end

struct WorkParams{AS1<:AbstractVector,AT2<:AbstractVector,
                              AT3<:AbstractWeights,AS2<:AbstractWeights,
                              AT4<:AbstractGraph,
                              }
    neighbor_idxs::AS1
    fitnesses::AT2
    weights_float::AT3
    weights_int::AS2
    interaction_graph::AT4
end

function calculate_node_payoff(graph::GNNGraph, node::Integer)
    payoff = 0.0
    for game in graph.ndata.in_edges[node]
       payoff += graph.edata.payoffs[game]
    end
    return payoff
end

function WorkParams(lmi::Moran{N},
                                initial_actions::AbstractVector) where {N}
    neighbor_idxs = Vector{Int64}(undef, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))

    interaction_graph = GNNGraph(lmi.interaction_adj_matrix)
    interaction_graph.ndata.strategy = initial_actions
    nidx = 1:nv(interaction_graph)

    edges_graph = edges(interaction_graph)
    df = DataFrame(:src => src.(edges_graph), :dst => dst.(edges_graph),
		   :edge_num => 1:ne(interaction_graph))
    interaction_graph.ndata.in_edges = coalesce.(leftjoin(DataFrame(dst=1:nv(interaction_graph)),
		       combine(groupby(df,:dst),:edge_num => (x -> [x]) => :edges);
		       on=:dst,order=:left),[Int[]])[:,:edges]
    interaction_graph.ndata.out_edges = coalesce.(leftjoin(DataFrame(src=1:nv(interaction_graph)),
		       combine(groupby(df,:src),:edge_num => (x -> [x]) => :edges);
		       on=:src,order=:left),[Int[]])[:,:edges]

    interaction_graph.edata.payoffs = Vector{Float64}(undef,ne(interaction_graph))
    for (edge_idx, edge) in enumerate(edges(interaction_graph))
      src_strategy = interaction_graph.ndata.strategy[src(edge)]
      dst_strategy = interaction_graph.ndata.strategy[dst(edge)]
      payoff = lmi.payoff_matrix[dst_strategy, src_strategy]*get_edge_weight(interaction_graph)[edge_idx]
      interaction_graph.edata.payoffs[edge_idx] = payoff
    end

    interaction_graph.ndata.payoffs = Vector{Float64}(undef,nv(interaction_graph))
    for node in vertices(interaction_graph)
      interaction_graph.ndata.payoffs[node] = calculate_node_payoff(interaction_graph, node)
    end

    return WorkParams(neighbor_idxs, fitnesses,
                                  weights_float, weights_int,
                                  interaction_graph,
				  )
end

function update_aux!(aux::WorkParams,
                     death_idx::Integer,
                     lmi::Moran{N}) where {N}
    # Update focal player's payoffs
    in_neighbors = inneighbors(aux.interaction_graph, death_idx)
    death_strategy = aux.interaction_graph.ndata.strategy[death_idx]
    edge_list = edges(aux.interaction_graph)
    for edge_idx in aux.interaction_graph.ndata.in_edges[death_idx]
      neighbor = src(edge_list[edge_idx])
      neighbor_strategy = aux.interaction_graph.ndata.strategy[neighbor]
      payoff = lmi.payoff_matrix[death_strategy,neighbor_strategy]*get_edge_weight(aux.interaction_graph)[edge_idx]
      aux.interaction_graph.edata.payoffs[edge_idx] = payoff
    end
    aux.interaction_graph.ndata.payoffs[death_idx] = calculate_node_payoff(aux.interaction_graph, death_idx)

    # Update payoffs of focal player's neighbors
    for edge_idx in aux.interaction_graph.ndata.out_edges[death_idx]
      neighbor = dst(edge_list[edge_idx])
      neighbor_strategy = aux.interaction_graph.ndata.strategy[neighbor]
      aux.interaction_graph.edata.payoffs[edge_idx] = lmi.payoff_matrix[neighbor_strategy,death_strategy]*get_edge_weight(aux.interaction_graph)[edge_idx]
      aux.interaction_graph.ndata.payoffs[neighbor] = calculate_node_payoff(aux.interaction_graph, neighbor)
    end
end

function play!(rng::AbstractRNG,
               lmi::Moran{N},
               aux::WorkParams) where {N}
    # Payoff flows along and is weighted by interaction_adj_matrix
    # Node is selected for reproduction using exponential fitness weighting
    # One of its outedges is selected for replacement, weighted by edge reproduction_adj_matrix
    # Note: each node on reproduction graph must have out-degree greater than
    #   zero, otherwise there is no way to choose node for replacement

    # Calculate all fitnesses
    # Note: normally, the fitness is exp(beta*(payoffs - mean(payoffs)))
    # However, if some payoffs are much larger than the mean
    # exponentiation can cause float overflow to Inf, which will cause errors in StatsBase.sample
    # Therefore, change (payoffs - mean(payoffs)) to (payoffs - maximum(payoffs));
    # this will ensure all fitnesses are less than one,
    # but the relative ratio of each fitness is unchanged, so the weighted sampling is unchanged
    payoffs = aux.interaction_graph.ndata.payoffs
    aux.fitnesses .= exp.(lmi.beta .* (payoffs .- maximum(payoffs)))

    # Choose focal (birth) node
    update_weights!(aux.weights_float, aux.fitnesses)
    focal_idx = sample(rng, 1:N, aux.weights_float)

    # Get list of out neighbors
    @. aux.neighbor_idxs = @view lmi.reproduction_adj_matrix_transpose[:, focal_idx]
    # Choose death node
    update_weights!(aux.weights_int, aux.neighbor_idxs)
    death_idx = sample(rng, 1:N, aux.weights_int)

    # Check for mutation
    mutation = false
    if rand(rng) <= lmi.epsilon
        new_strategy = rand(rng, 1:(lmi.num_actions))
	mutation = true
    else
        new_strategy = aux.interaction_graph.ndata.strategy[focal_idx]
    end

    # Check if updates are necessary
    if aux.interaction_graph.ndata.strategy[death_idx] == new_strategy
	    return mutation
    end

    # Apply spatial Moran process
    aux.interaction_graph.ndata.strategy[death_idx] = new_strategy

    # Update aux variables
    update_aux!(aux, death_idx, lmi)

    return mutation
end

function count_actions!(counts::AbstractVector{N}, new_values::AbstractVector{N}) where {N}
    for (i, _) in enumerate(counts)
        counts[i] += count(==(i), new_values)
    end
    return counts
end

function time_series(lmi::Moran{N}, ts_length::Integer,
                     seed::Integer=12345) where {N}
    rng = Xoshiro(seed)
    initial_actions = rand(rng, 1:(lmi.num_actions), N)
    aux = WorkParams(lmi, initial_actions)
    out = Matrix{Int}(undef, N, ts_length + 1)
    steps_following_mutation = Integer[]
    for i in 1:N
        out[i, 1] = initial_actions[i]
    end
    for t in 1:ts_length
        mutation = play!(rng, lmi, aux)
        out[:, t + 1] = aux.interaction_graph.ndata.strategy
	if mutation
		push!(steps_following_mutation, t+1)
	end
    end
    return out, steps_following_mutation
end

function cumulative(lmi::Moran{N}, ts_length::Integer,
                    seed::Integer=12345) where {N}
    rng = Xoshiro(seed)
    initial_actions = rand(rng, 1:(lmi.num_actions), N)
    cumulative = zeros(Int, lmi.num_actions)
    count_actions!(cumulative, initial_actions)
    aux = WorkParams(lmi, initial_actions)
    for t in 2:ts_length
        play!(rng, lmi, aux)
        count_actions!(cumulative, aux.interaction_graph.ndata.strategy)
    end
    return cumulative
end
