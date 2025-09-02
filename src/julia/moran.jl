using GameTheory
using Random
using StatsBase
using Graphs

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

struct WorkParams{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,
                              AT3<:AbstractWeights,AS2<:AbstractWeights,
                              AT4<:AbstractMatrix,AT5<:AbstractMatrix,AT6<:AbstractMatrix,
                              AT7<:AbstractMatrix}
    neighbor_idxs::AS1
    payoffs::AT1
    payoffs_transpose::AT7
    fitnesses::AT2
    weights_float::AT3
    weights_int::AS2
    payoffs_rowplayer_per_strategy::AT4
    payoffs_colplayer_per_strategy::AT5
    payoffs_player_pairwise_transpose::AT6
end

function WorkParams(lmi::Moran{N},
                                initial_actions::AbstractVector) where {N}
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    payoffs_transpose = Matrix{Float64}(undef, 1, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    payoffs_rowplayer_per_strategy = lmi.payoff_matrix[initial_actions, :]
    payoffs_colplayer_per_strategy = lmi.payoff_matrix_transpose[initial_actions, :]
    payoffs_player_pairwise_transpose = payoffs_colplayer_per_strategy[:,
                                                                       initial_actions] .*
                                        lmi.interaction_adj_matrix
    return WorkParams(neighbor_idxs, payoffs, payoffs_transpose, fitnesses,
                                  weights_float, weights_int,
                                  payoffs_rowplayer_per_strategy,
                                  payoffs_colplayer_per_strategy,
                                  payoffs_player_pairwise_transpose)
end

function calc_payoffs!(aux::WorkParams,
                       actions::AbstractVector,
                       lmi::Moran{N}) where {N}
    # Calculate payoffs
    sum!(aux.payoffs_transpose, aux.payoffs_player_pairwise_transpose)
    return transpose!(aux.payoffs, aux.payoffs_transpose)
end

function update_aux!(aux::WorkParams,
                     new_strategy::Integer,
                     death_idx::Integer,
                     lmi::Moran{N}) where {N}
    # Update various matrices
    @. aux.payoffs_rowplayer_per_strategy[death_idx, :] = @view lmi.payoff_matrix_transpose[:,
                                                                                            new_strategy]
    @. aux.payoffs_colplayer_per_strategy[death_idx, :] = @view lmi.payoff_matrix[:,
                                                                                  new_strategy]
    @. aux.payoffs_player_pairwise_transpose[:, death_idx] = @view aux.payoffs_colplayer_per_strategy[:,
                                                                                                      new_strategy]
    for idx in 1:N
        aux.payoffs_player_pairwise_transpose[idx, death_idx] *= lmi.interaction_adj_matrix[idx,
                                                                                            death_idx]
    end
    @. aux.payoffs_player_pairwise_transpose[death_idx, :] = @view aux.payoffs_rowplayer_per_strategy[:,
                                                                                                      new_strategy]
    for idx in 1:N
        aux.payoffs_player_pairwise_transpose[death_idx, idx] *= lmi.interaction_adj_matrix[death_idx,
                                                                                            idx]
    end
end

function play!(actions::AbstractVector,
               rng::AbstractRNG,
               lmi::Moran{N},
               aux::WorkParams) where {N}
    # Payoff flows along and is weighted by interaction_adj_matrix
    # Node is selected for reproduction using exponential fitness weighting
    # One of its outedges is selected for replacement, weighted by edge reproduction_adj_matrix
    # Note: each node on reproduction graph must have out-degree greater than
    #   zero, otherwise there is no way to choose node for replacement

    # Calculate all payoffs
    # Note: we use auxiliary variables to save intermediate calculations that get reused
    # In particular, for n players and m strategies, we have
    # payoffs[i] {payoffs[i] is payoff of player i, nx1 dim}
    # = sum_j=1^n payoffs_pairwise[i,j] {payoffs_pairwise[i,j] is payoff from player j to i, possibly zero; called "payoffs_player_pairwise"}
    # = sum_j=1^n payoffs_pairwise[i,j] * ones[j] {ones is nx1 vector of ones}
    # = sum_j=1^n (payoffs_connected[i,j] * adj_matrix[i.j]) * ones[j] {payoffs_connected[i,j] is payoff from player j to i assuming they were connected; called "payoffs_colplayer_per_strategy"}
    # = sum_j=1^n (strat_payoffs[strategy[i],strategy[j]] * adj_matrix[i.j]) * ones[j] {adj_matrix[i,j] is weight of connection from player j to player i,
    # 	  strat_payoffs[alpha,beta] is mxm payoff matrix from strategy alpha to strategy beta, strategy[i] player i's strategy}
    # = sum_j=1^n sum_alpha=1^m sum_beta=1^m ((e_j[alpha] * strat_payoffs[alpha,beta] * e_i[beta]) .* adj_matrix[i.j]) * ones[j]
    # 	  e[i] is m x 1 basis vector with 1 at index of player i's strategy}
    # Or, converting to vector notation
    # payoffs = ((e^T * strat_payoffs * e) .* adj_matrix) * ones {e is mxn matrix with ith column the mx1 basis vector with 1 at index of player i's strategy, and .* the Hadamard product}
    #
    # However, since at most the death player's strategy changes each round, only a single column of e changes each round;
    # therefore, after initially calculating the matrix payoffs_player_pairwise,
    # updates to it can be accomplished in three pairs of steps (pairs since both e and e^T need updates):
    # First, update strat_payoff*e and e^Tstrat_payoff:
    #   one col of e^T * strat_payoffs {denoted "payoffs_rowplayer_per_strategy"} as (e^T*strat_payoffs)[alpha,death_idx] -> payoff_matrix[new_strat,alpha] for alpha=1:m
    #   one row of strat_payoffs*e {denoted "payoffs_colplayer_per_strategy"} as (strat_payoffs*e)[death_idx,alpha] -> payoff_matrix[alpha,new_strat] for alpha=1:m
    # Next, copy the new row/column into payoffs_pairwise (notice the reversed indexing of e^T*strat_payoffs and strat_payoffs*e compared to prior step)
    #   one col of payoffs_pairwise[i,death_idx] -> (e^T*strat_payoffs)[new_strat,i]
    #   one row of payoffs_pairwise[death_idx,j] -> (strat_payoffs*e)[i,new_strat]
    # Finally, multiply payoffs_pairwise by adj_matrix
    #   one col of payoffs_pairwise[i,death_idx] -> payoffs_pairwise[i,death_idx] .* adj_matrix[i,death_idx]
    #   one row of payoffs_pairwise[death_idx,i] -> payoffs_pairwise[death_idx,i] .* adj_matrix[death_idx,i]
    # Then, the * one[j] sum must be recomputed each round
    # Note: some steps are transposed in implementation to better align with Julia's column-major ordering
    calc_payoffs!(aux, actions, lmi)

    # Calculate all fitnesses
    # Note: normally, the fitness is exp(beta*(payoffs - mean(payoffs)))
    # However, if some payoffs are much larger than the mean
    # exponentiation can cause float overflow to Inf, which will cause errors in StatsBase.sample
    # Therefore, change (payoffs - mean(payoffs)) to (payoffs - maximum(payoffs));
    # this will ensure all fitnesses are less than one,
    # but the relative ratio of each fitness is unchanged, so the weighted sampling is unchanged
    aux.fitnesses .= exp.(lmi.beta .* (aux.payoffs .- maximum(aux.payoffs)))

    # Choose focal (birth) node
    update_weights!(aux.weights_float, aux.fitnesses)
    focal_idx = sample(rng, 1:N, aux.weights_float)

    # Get list of neighbors
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
        new_strategy = actions[focal_idx]
    end

    # Apply spatial Moran process
    actions[death_idx] = new_strategy

    # Update aux variables
    update_aux!(aux, new_strategy, death_idx, lmi)

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
    actions = rand(rng, 1:(lmi.num_actions), N)
    aux = WorkParams(lmi, actions)
    out = Matrix{Int}(undef, N, ts_length + 1)
    steps_following_mutation = Integer[]
    for i in 1:N
        out[i, 1] = actions[i]
    end
    for t in 1:ts_length
        mutation = play!(actions, rng, lmi, aux)
        out[:, t + 1] = actions
	if mutation
		push!(steps_following_mutation, t+1)
	end
    end
    return out, steps_following_mutation
end

function cumulative(lmi::Moran{N}, ts_length::Integer,
                    seed::Integer=12345) where {N}
    rng = Xoshiro(seed)
    actions = rand(rng, 1:(lmi.num_actions), N)
    cumulative = zeros(Int, lmi.num_actions)
    count_actions!(cumulative, actions)
    aux = WorkParams(lmi, actions)
    for t in 2:ts_length
        play!(actions, rng, lmi, aux)
        count_actions!(cumulative, actions)
    end
    return cumulative
end
