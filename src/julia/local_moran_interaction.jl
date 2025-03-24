module local_moran_interaction

using GameTheory
using LinearAlgebra
using BlockArrays
using SparseArrays
using Random
using Statistics
using StatsBase
using Symbolics
using BenchmarkTools
using InteractiveUtils
using LaTeXStrings
using DataToolkit
using DrWatson
using Graphs
using SimpleWeightedGraphs
using DataFrames
using DataFramesMeta
using CairoMakie
using GraphMakie
using NetworkLayout
using Polyhedra
using Colors
using Memoize
using DimensionalData
using SplitApplyCombine
using CondaPkg
using PythonCall
using CSV
using JLD2
#using CombinatorialMultiGrid

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

struct LocalMoranInteraction{N,T1<:Real,S<:Integer,S2<:Integer,T2<:Real,T3<:Real,
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

function LocalMoranInteraction(g::NormalFormGame,
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
    return LocalMoranInteraction(players, num_actions, g.players[1].payoff_array,
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

abstract type WorkParams end

struct WorkParamsLoop{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,
                      AT3<:AbstractWeights,AS2<:AbstractWeights} <: WorkParams
    neighbor_idxs::AS1
    payoffs::AT1
    fitnesses::AT2
    weights_float::AT3
    weights_int::AS2
end

function WorkParamsLoop(lmi::LocalMoranInteraction{N}) where {N}
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    return WorkParamsLoop(neighbor_idxs, payoffs, fitnesses, weights_float, weights_int)
end

struct WorkParamsSingleUpdate{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,
                              AT3<:AbstractWeights,AS2<:AbstractWeights,
                              AT4<:AbstractMatrix,AT5<:AbstractMatrix,AT6<:AbstractMatrix,
                              AT7<:AbstractMatrix} <: WorkParams
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

function WorkParamsSingleUpdate(lmi::LocalMoranInteraction{N},
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
    return WorkParamsSingleUpdate(neighbor_idxs, payoffs, payoffs_transpose, fitnesses,
                                  weights_float, weights_int,
                                  payoffs_rowplayer_per_strategy,
                                  payoffs_colplayer_per_strategy,
                                  payoffs_player_pairwise_transpose)
end

struct WorkParamsMatrix{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,
                        AT3<:AbstractWeights,AS2<:AbstractWeights,
                        AS3<:AbstractMatrix,AT4<:AbstractMatrix,AS4<:AbstractVector} <:
       WorkParams
    neighbor_idxs::AS1
    payoffs::AT1
    fitnesses::AT2
    weights_float::AT3
    weights_int::AS2
    actions_matrix::AS3
    actions_matrix_T::AS3
    work_array_1::AT4
    work_array_2::AT4
    work_array_3::AT4
    ones_matrix::AS4
end

function WorkParamsMatrix(lmi::LocalMoranInteraction{N}) where {N}
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    actions_matrix = Matrix{Int64}(undef, N, lmi.num_actions)
    actions_matrix_T = Matrix{Int64}(undef, lmi.num_actions, N)
    work_array_1 = Matrix{Float64}(undef, N, lmi.num_actions)
    work_array_2 = Matrix{Float64}(undef, N, N)
    work_array_3 = Matrix{Float64}(undef, N, N)
    ones_matrix = ones(Int64, N)
    return WorkParamsMatrix(neighbor_idxs, payoffs, fitnesses, weights_float, weights_int,
                            actions_matrix, actions_matrix_T, work_array_1, work_array_2,
                            work_array_3, ones_matrix)
end

function calc_payoffs!(aux::WorkParamsLoop,
                       actions::AbstractVector,
                       lmi::LocalMoranInteraction{N}) where {N}
    @inbounds for focal_idx in 1:length(actions)
        focal_strategy = actions[focal_idx]
        total = 0.0
        @inbounds for neighbor_idx in 1:length(actions)
            neighbor_strategy = actions[neighbor_idx]
            total += lmi.payoff_matrix_transpose[neighbor_strategy, focal_strategy] *
                     lmi.interaction_adj_matrix[neighbor_idx, focal_idx]
        end
        aux.payoffs[focal_idx] = total
    end
end

function calc_payoffs!(aux::WorkParamsSingleUpdate,
                       actions::AbstractVector,
                       lmi::LocalMoranInteraction{N}) where {N}
    # Calculate payoffs
    sum!(aux.payoffs_transpose, aux.payoffs_player_pairwise_transpose)
    return transpose!(aux.payoffs, aux.payoffs_transpose)
end

function update_aux!(aux::WorkParams,
                     new_strategy::Integer,
                     death_idx::Integer,
                     lmi::LocalMoranInteraction{N}) where {N}
end

function update_aux!(aux::WorkParamsSingleUpdate,
                     new_strategy::Integer,
                     death_idx::Integer,
                     lmi::LocalMoranInteraction{N}) where {N}
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

function calc_payoffs!(aux::WorkParamsMatrix,
                       actions::AbstractVector,
                       lmi::LocalMoranInteraction{N}) where {N}
    @inbounds for col in 1:(lmi.num_actions), row in 1:N
        aux.actions_matrix[row, col] = 0
    end
    @inbounds for player_idx in 1:N
        strategy_idx = actions[player_idx]
        aux.actions_matrix[player_idx, strategy_idx] = 1
    end
    transpose!(aux.actions_matrix_T, aux.actions_matrix)
    mul!(aux.work_array_1, aux.actions_matrix, lmi.players[1].payoff_array)
    mul!(aux.work_array_2, aux.work_array_1, aux.actions_matrix_T)
    aux.work_array_3 .= aux.work_array_2 .* lmi.interaction_adj_matrix_transpose
    return mul!(aux.payoffs, aux.work_array_3, aux.ones_matrix)
end

function play!(actions::AbstractVector,
               rng::AbstractRNG,
               lmi::LocalMoranInteraction{N},
               aux::WorkParams) where {N}
    # Payoff flows along and is weighted by interaction_adj_matrix
    # Node is selected for reproduction using exponential fitness weighting
    # One of its outedges is selected for replacement, weighted by edge reproduction_adj_matrix
    # Note: each node on reproduction graph must have out-degree greater than
    #   zero, otherwise there is no way to choose node for replacement

    # Calculate all payoffs
    calc_payoffs!(aux, actions, lmi)

    # Calculate all fitnesses
    aux.fitnesses .= exp.(lmi.beta .* (aux.payoffs .- mean(aux.payoffs)))

    # Choose focal (birth) node
    update_weights!(aux.weights_float, aux.fitnesses)
    focal_idx = sample(rng, 1:N, aux.weights_float)

    # Get list of neighbors
    @. aux.neighbor_idxs = @view lmi.reproduction_adj_matrix_transpose[:, focal_idx]
    # Choose death node
    update_weights!(aux.weights_int, aux.neighbor_idxs)
    death_idx = sample(rng, 1:N, aux.weights_int)

    # Check for mutation
    if rand(rng) <= lmi.epsilon
        new_strategy = rand(rng, 1:(lmi.num_actions))
    else
        new_strategy = actions[focal_idx]
    end

    # Apply spatial Moran process
    actions[death_idx] = new_strategy

    # Update aux variables
    update_aux!(aux, new_strategy, death_idx, lmi)

    return nothing
end

function count!(counts::AbstractVector{N}, new_values::AbstractVector{N}) where {N}
    for (i, _) in enumerate(counts)
        counts[i] += count(==(i), new_values)
    end
    return counts
end

function time_series(lmi::LocalMoranInteraction{N}, ts_length::Integer,
                     payoff_counting_method::String="single-update",
                     seed::Integer=12345) where {N}
    rng = Xoshiro(seed)
    actions = rand(rng, 1:(lmi.num_actions), N)
    if payoff_counting_method == "loop"
        aux = WorkParamsLoop(lmi)
    elseif payoff_counting_method == "single-update"
        aux = WorkParamsSingleUpdate(lmi, actions)
    elseif payoff_counting_method == "matrix"
        aux = WorkParamsMatrix(lmi)
    else
        throw(ArgumentError("payoff_counting_method must be a string in set [\"loop\", \"single-update\", \"matrix\"]"))
    end
    out = Matrix{Int}(undef, N, ts_length + 1)
    for i in 1:N
        out[i, 1] = actions[i]
    end
    for t in 1:ts_length
        play!(actions, rng, lmi, aux)
        out[:, t + 1] = actions
    end
    return out
end

function cumulative(lmi::LocalMoranInteraction{N}, ts_length::Integer,
                    payoff_counting_method::String="single-update",
                    seed::Integer=12345) where {N}
    rng = Xoshiro(seed)
    actions = rand(rng, 1:(lmi.num_actions), N)
    cumulative = zeros(Int, lmi.num_actions)
    count!(cumulative, actions)
    if payoff_counting_method == "loop"
        aux = WorkParamsLoop(lmi)
    elseif payoff_counting_method == "single-update"
        aux = WorkParamsSingleUpdate(lmi, actions)
    elseif payoff_counting_method == "matrix"
        aux = WorkParamsMatrix(lmi)
    else
        throw(ArgumentError("payoff_counting_method must be a string in set [\"loop\", \"single-update\", \"matrix\"]"))
    end
    for t in 2:ts_length
        play!(actions, rng, lmi, aux)
        count!(cumulative, actions)
    end
    return cumulative
end

function check_formulas(nb_phases::Integer, nb_players::Integer)
    # Define symbols
    @variables B_0 beta_0 phi delta c i j k q r d n

    # Define constants
    d = nb_phases
    n = nb_players

    # Define functions
    B(phi) = B_0 * (1 + cos(phi)) / 2
    beta(phi) = beta_0 * (1 + cos(phi)) / 2

    # Define fixation probabilities
    function rho_CC(phi)
        return 1 / (1 + sum([exp(k * (k + 1 - n) * delta / (n - 1) * (B(phi) - B_0))
                             for k in 1:(n - 1)]))
    end
    function rho_NC(phi)
        return 1 / (1 + sum([exp(delta / (n - 1) *
                                 (j^2 * (beta(phi) - 1 / 2 * B_0) +
                                  j * (1 / 2 * B_0 - (n - 1) * beta(phi) + (n - 1) * c)))
                             for j in 1:(n - 1)]))
    end
    omega = 1 / exp(delta * ((n - 1) * c - (n - 2) / 2 * B_0))
    rho_CN(phi) = rho_NC(phi) / omega

    # Define matrices
    B_1 = [i != j ?
           rho_CC(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) :
           1 - sum([rho_CC(floor(d / 2) - abs(floor(d / 2) - abs(delta_phi_prime)))
                    for delta_phi_prime in 1:(d - 1)])
           for i in 1:d, j in 1:d]
    B_2 = [rho_CN(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) for i in 1:d, j in 1:d]
    B_3 = [rho_NC(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) for i in 1:d, j in 1:d]
    B_4 = [1 / n for i in 1:d, j in 1:d]

    # Define full Markov transition matrix
    M = mortar([[B_1, B_3] [B_2, B_4]])

    # Define stationary distribution
    # Tripp (2024) shows s_1, s_2 independent of q, so f
    q = 1
    s_1 = sum([M[r, q] for r in (d + 1):(2d)]) /
          (d * sum([M[q, r] + M[r, q] for r in (d + 1):(2d)]))
    s_2 = sum([M[q, r] for r in (d + 1):(2d)]) /
          (d * sum([M[q, r] + M[r, q] for r in (d + 1):(2d)]))
    return s = hcat([s_1 * ones(d), s_2 * ones(d)])

    ## Check that the explicitly calculated stationary distribution matches the eigenvector
    #eigen = eigvals(M)[1]
    #if s != eigen
    #    throw(ErrorException(\"Explicitly calculated stationary distribution differs from eigenvector\"))

    ## Check that fraction of communicative simplifies correctly
    #if not isequal(simplify(d*s_1, threaded=true), 1/(1+1/omega))
    #    throw(ErrorException(\"Cooperative fraction does not simplify to known correct answer 1/(1+1/omega)))
end

function extract_num_communicative(players_per_strategy::AbstractVector{<:Integer})
    @assert iseven(length(players_per_strategy))
    # Define matrix that extracts number of communicative players
    N = div(length(players_per_strategy), 2)
    communicative_matrix = [ones(Int, N); zeros(Int, N)]
    # communicatative have value < nb_phases
    num_communicative = dot(players_per_strategy, communicative_matrix)
    return num_communicative
end

function fraction_communicative(cumulative_populations, time_steps, nb_players)
    nb_communicative = extract_num_communicative(cumulative_populations)
    fraction_communicative = nb_communicative ./ (time_steps * nb_players)
    return fraction_communicative
end

@enum GameType harmony chicken battle hero compromise concord staghunt dilemma deadlock assurance coordination peace

const game_type_colors = Dict(instances(GameType) .=>
                                  distinguishable_colors(length(instances(GameType))))
const game_type_full_names = Dict(chicken => "Snowdrift", # Also called chicken
                                  battle => "Battle of the Sexes",
                                  hero => "Hero",
                                  compromise => "Compromise",
                                  deadlock => "Deadlock",
                                  dilemma => "Prisoner's Dilemma",
                                  staghunt => "Mutualism", # Also called Staghunt
                                  assurance => "Assurance",
                                  coordination => "Coordination",
                                  peace => "Peace",
                                  harmony => "Harmony",
                                  concord => "Concord")

function game_type_inequalities(R::Real, S::Real, T::Real, P::Real)
    # Break ties in direction T < R < S < P
    # Note: swap columns from Brun 2015 to convert their convention of
    #   [ (C,N)   (C,C) ]
    #   [ (N,N)   (N,C) ]
    # to our convention of
    #   [ (C,C)   (C,N) ]
    #   [ (N,C)   (N,N) ]
    return Dict{GameType,Vector}(chicken => [R <= T, T <= P, P <= S],
                                 battle => [R <= P, P < T, T <= S],
                                 hero => [P < R, R <= T, T <= S],
                                 compromise => [P < T, T < R, R <= S],
                                 deadlock => [T <= P, P < R, R <= S],
                                 dilemma => [T < R, R <= P, P <= S],
                                 staghunt => [T < R, R <= S, S < P],
                                 assurance => [T <= S, S < R, R <= P],
                                 coordination => [S < T, T < R, R <= P],
                                 peace => [S < R, R <= T, T <= P],
                                 harmony => [R <= S, S < T, T <= P],
                                 concord => [R <= T, T <= S, S < P])
end

function combine_communicative_noncommunicative(players_per_strategy::AbstractVector{<:Integer},
                                                nb_phases::Integer)
    # Combine communicatative (value < nb_phases)
    # and non-communicatative (values >= nb_phases)
    combine_com_noncom_matrix = hcat(Matrix{Int}(I, nb_phases, nb_phases),
                                     Matrix{Int}(I, nb_phases, nb_phases))
    phase_idxs = combine_com_noncom_matrix * players_per_strategy
    return phase_idxs
end

function swap_strategies!(payoff_matrix::AbstractMatrix)
    return payoff_matrix[2, 2], payoff_matrix[2, 1], payoff_matrix[1, 2], payoff_matrix[1, 1] = payoff_matrix[1,
                                                                                                              1],
                                                                                                payoff_matrix[1,
                                                                                                              2],
                                                                                                payoff_matrix[2,
                                                                                                              1],
                                                                                                payoff_matrix[2,
                                                                                                              2]
end

function canonical_payoff!(payoff_matrix::AbstractMatrix)
    # Using our convention that symmetric games imply the col-player's
    # payoff matrix is the transpose of the provided (row-player's) payoff
    # matrix, apply Robinson & Goforth's convention of highest row-player payoff in
    # right-column and highest col-player payoff in lower-row (note:
    # col-player criteria is switched to match our transpose-convention
    # above) by swapping rows and columns
    # For symmetric games, this is equivalent to renaming strategies (C <->
    # N) for both players simultaneously to ensure the max payoff is in the
    # right column
    # Note: assumes there is a unique maximum payoff
    max_index_orig = argmax(payoff_matrix)
    # Ensure highest row-player payoff is in right-column
    if max_index_orig[2] == 1
        swap_strategies!(payoff_matrix)
    end
end

function game_type(payoff_matrix::AbstractMatrix)
    # Canonicalize
    canonical_payoff!(payoff_matrix)
    # Check ordering of payoffs
    R = payoff_matrix[1, 1]
    S = payoff_matrix[1, 2]
    T = payoff_matrix[2, 1]
    P = payoff_matrix[2, 2]
    for game_type in instances(GameType)
        if (&)(game_type_inequalities(R, S, T, P)[game_type]...)
            return game_type
        end
    end
    # Non-strict inequality implying one or more ties
    throw(ErrorException("Should never hit this since tie-breaking is handled with non-strict inequalities"))
end

@memoize function game_types_per_strategy_pair(mutual_benefit_synchronous::Real,
                                               unilateral_benefit_synchronous::Real,
                                               cost::Real,
                                               symmetry_breaking::Real,
                                               nb_phases::Integer)

    # Define payoff response submatrices
    payoff = payoff_matrix(nb_phases, mutual_benefit_synchronous,
                           unilateral_benefit_synchronous, cost; symmetry_breaking)

    # Define game type per strategy pair
    game_types = Matrix{GameType}(undef, nb_phases, nb_phases)
    for idx in eachindex(IndexCartesian(), game_types)
        strategy_offset = CartesianIndices((0:nb_phases:nb_phases, 0:nb_phases:nb_phases))
        indices = idx .+ strategy_offset
        game_types[idx] = game_type(@view payoff[indices])
    end
    return game_types
end

@enum StrategyParity all_communicative all_noncommunicative mixed

function check_all_same_strategy(strategies_per_player::AbstractVector{<:Integer},
                                 nb_phases::Integer)
    nb_players = length(strategies_per_player)
    players_per_strategy = extract_counts(strategies_per_player, nb_phases)
    # Check if all players were communicative/noncommunicative
    nb_communicative = extract_num_communicative(players_per_strategy)
    if nb_communicative == 0
        return all_noncommunicative
    elseif nb_communicative == nb_players
        return all_communicative
    end
    return mixed
end

function extract_most_common_game_types(strategies_per_player::AbstractVector{<:Integer},
                                        mutual_benefit_synchronous::Real,
                                        unilateral_benefit_synchronous::Real,
                                        cost::Real,
                                        symmetry_breaking::Real,
                                        nb_phases::Integer,
                                        interaction_adj_matrix::AbstractMatrix{<:Integer})
    players_per_phase = mod1.(strategies_per_player, nb_phases)

    # Get game type of each startegy interaction pair
    game_types = game_types_per_strategy_pair(mutual_benefit_synchronous,
                                              unilateral_benefit_synchronous, cost,
                                              symmetry_breaking, nb_phases)

    # Count game types
    game_counts = Dict{GameType,Integer}(instances(GameType) .=> 0)
    for (cart_idx, value) in pairs(interaction_adj_matrix)
        if value == 0
            continue
        end
        (row, col) = Tuple(cart_idx)
        row_phase = players_per_phase[row]
        col_phase = players_per_phase[col]
        game_type = game_types[row_phase, col_phase]
        game_counts[game_type] += value
    end

    # Find most common game type
    most_common_game_type = findmax(game_counts)[2]

    return most_common_game_type
end

function extract_phases(players_per_phase::AbstractVector{<:Integer}, nb_phases::Integer)
    # Convert phase idxs to phases
    phases = 2 * pi / nb_phases .* players_per_phase
    return phases
end

function extract_order_parameters(players_per_strategy::AbstractVector{<:Integer},
                                  nb_phases::Integer)
    # Combine communicative and noncommunicative
    phase_indxs = combine_communicative_noncommunicative(players_per_strategy, nb_phases)
    # Convert strategy idxs to phases
    phases = extract_phases(phase_indxs, nb_phases)
    order_parameters = abs(mean(exp.(im * phases)))
    return order_parameters
end

function extract_counts(strategies_per_player::AbstractVector{<:Integer}, nb_phases::Integer)
    nb_strategies = 2*nb_phases
    counts = zeros(Int, nb_strategies)
    count!(counts, strategies_per_player)
    return counts
end

function generate_communities(graph::AbstractSimpleWeightedGraph)
    ## Label propagation
    #communities = label_propagation(graph; rng=Xoshiro(12345))[1]

    ## Strongly connected components
    #connected_components = strongly_connected_components(graph)
    #communities = Vector{Int64}(undef, nv(graph))
    #for (community, idxs) in pairs(connected_components)
    #    for idx in idxs
    #        communities[idx] = community
    #    end
    #end

    # InfoMap
    CSV.write(datadir("processed","InfoMapOutput","c-elegans-network.txt"), edges(graph); writeheader=false, delim=" ")
    CondaPkg.add("infomap")
    infomap = pyimport("infomap")
    infomap.Infomap(infomap.Config("-d -2 --preferred-number-of-modules 2 --variable-markov-time $(datadir("processed","InfoMapOutput","c-elegans-network.txt")) $(datadir("processed","InfoMapOutput"))",true)).run()
    df = DataFrame(CSV.File(datadir("processed","InfoMapOutput","c-elegans-network.tree");comment="#",delim=" ",header=["path","flow","name","node_id"]))
    df[!,"community"] = parse.(Int64,(map(x -> x[1], split.(df[!,"path"],":"))))
    df_new = df[!,["node_id","community"]]
    sort!(df_new, "node_id")
    communities = df_new[!, "community"]

    # Covariance
    config = Dict("B_factor"=>1.5,"adj_matrix_source"=>"c-elegans","nb_phases"=>20,"payoff_update_metho
		  d"=>"single-update","selection_strength"=>0.2,"symmetry_breaking"=>1.0,"time_steps"=>8000000)
    results = wload(datadir("raw","timeseries",savename(config,"jld2")))
    data = DimArray(results["all_populations"], (:player_index, :time_step))
    covariances = cov(data, dims=:time_step)

    covariance_cutoff = 1500
    communities = 2*ones(Int64,nv(graph))
    communities[map(x -> x[2], findall(sum(covariances,dims=1) .>= covariance_cutoff))] .= 1

    return communities
end

function get_chimera_indices(data::DimArray,communities::AbstractVector{<:Integer},nb_phases::Integer)
   strategies_grouped = groupby(data, Dim{:player_index}=>(x -> communities[x]))
   phase_parameters = map(x -> extract_order_parameters.(extract_counts.(eachslice(x,dims=:time_step),nb_phases),nb_phases), strategies_grouped)
   metastability = mean(cov.(phase_parameters))
   chimera_index = mean(cov.(invert(phase_parameters)))

   results = Dict("metastability_index"=>metastability, "chimera_index"=>chimera_index)
   return results
end

function get_adj_matrices(adj_matrix_source::String; nb_players::Integer = 20, regular_degree::Integer = 10, rng::AbstractRNG=Xoshiro(1))
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
        interaction_adj_matrix = collect(round.(get_connectome()) .!= 0)
        reproduction_adj_matrix = collect((interaction_adj_matrix + I) .!= 0)
    elseif adj_matrix_source == "c-elegans-undirected"
        interaction_adj_matrix = round.(get_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
        reproduction_adj_matrix = interaction_adj_matrix + I
    elseif adj_matrix_source == "c-elegans-undirected-unweighted"
        interaction_adj_matrix = round.(get_connectome())
        interaction_adj_matrix = interaction_adj_matrix + transpose(interaction_adj_matrix)
        interaction_adj_matrix = collect(interaction_adj_matrix .!= 0)
        reproduction_adj_matrix = collect((interaction_adj_matrix + I) .!= 0)
    elseif adj_matrix_source == "c-elegans"
        interaction_adj_matrix = round.(get_connectome())
        reproduction_adj_matrix = interaction_adj_matrix + I
    else
        throw(ArgumentError("adj_matrix_source must be a string in set [\"well-mixed\", "
                            * "\"c-elegans\", \"c-elegans-unweighted\", "
                            *
                            "\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\"
			    * "\"random-regular-graph\", \"random-regular-digraph\"]"))
    end
    return interaction_adj_matrix, reproduction_adj_matrix
end

function calc_cumulative(config::Dict)
    # Unpack values
    @unpack selection_strength, symmetry_breaking, adj_matrix_source, payoff_update_method, time_steps, nb_phases = config

    # Define system
    cost = 0.1
    mutation_rate = 0.0001

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Define Bs on which to run
    B_crit = 2 * cost * (nb_players - 1) / (nb_players - 2)
    nb_Bs = 11
    step_size_Bs = 0.04
    Bs = B_crit .+ ((0:(nb_Bs - 1)) .- (nb_Bs - 1) / 2) .* step_size_Bs

    # Run the model for weak selection strength
    @time cumulative_populations = [cumulative(LocalMoranInteraction(NormalFormGame(payoff_matrix(nb_phases,
                                                                                                  B,
                                                                                                  0.95 *
                                                                                                  B,
                                                                                                  cost;
                                                                                                  symmetry_breaking)),
                                                                     interaction_adj_matrix,
                                                                     reproduction_adj_matrix,
                                                                     selection_strength,
                                                                     mutation_rate),
                                               time_steps, payoff_update_method)
                                    for B in Bs]

    # Plot fraction communcative
    nb_communicative = [extract_num_communicative(final_population)
                        for final_population in cumulative_populations]
    fraction_communicative = nb_communicative ./ (time_steps * nb_players)

    # Package results
    return @strdict(Bs, nb_players, cost, fraction_communicative)
end

function plot_cumulative(selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         payoff_update_method::String="single-update",
                         time_steps::Integer=2_000_000,
			 nb_phases::Integer=20)

    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps,
                      selection_strength, symmetry_breaking, nb_phases)
    data, _ = produce_or_load(calc_cumulative, config, datadir("raw","cumulative"))

    # Produce figure
    fig = generate_cumulative_plot(data, config)

    # Save figure
    filename = plotsdir("cumulative", savename(config, "png"))
    save(filename, fig)

    return fig
end

function analytic_frac_communicative(B0,beta0;selection_strength,cost,nb_players,symmetry_breaking,nb_phases)
	n = nb_players
	alpha = symmetry_breaking
	delta = selection_strength
	d = nb_phases

	beta(delta_phi) = beta0*(1+cos(delta_phi))/2
	denom(delta_phi) = 1 + sum([
			  exp(delta*(j^2*(beta(delta_phi)-B0/2)
					   +j*(B0/2+beta(delta_phi)*(1-2*alpha*n)+cost*(n-1))))
			 for j in 1:n-1])
	omega(delta_phi) = exp(delta*(n-1)*((n-1)*cost+n*beta(delta_phi)*(1-2*alpha)-(n-2)/2*B0))

	s1_partials = [1/denom(delta_phi) for delta_phi in pi/d*(1:d)]
	s2_partials = [omega(delta_phi) for delta_phi in pi/d*(1:d)] .* s1_partials
	nu = sum(s2_partials)/sum(s1_partials)
	frac_communicative = 1/(1+nu)
	return frac_communicative
end

function generate_cumulative_plot(data::Dict, config::Dict)
    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1];
	      title=L"Selection $\delta = %$(round(config[\"selection_strength\"],sigdigits=2))$",
              xlabel=L"Maximum benefit of mutual communication, $B(0)$",
              ylabel="Frequency of communicative strategies",
              limits=(nothing, nothing, 0, 1))
    scatter!(ax, data["Bs"], data["fraction_communicative"]; label="Simulation")
    beta0(B0) = 0.95*B0

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(config["adj_matrix_source"])
    graph = SimpleDiGraph(interaction_adj_matrix)
    nb_effective = ( mean(indegree(graph))+1) # Add one since n=degree+1 for well-mixed case

    lines!(ax, data["Bs"][begin] .. data["Bs"][end],
	   B0 -> analytic_frac_communicative(B0, beta0(B0);
			selection_strength=config["selection_strength"], cost=data["cost"], nb_players=nb_effective,
			symmetry_breaking=config["symmetry_breaking"], nb_phases=config["nb_phases"]); label="Theory",
           color=:orange)
    lines!(ax, data["Bs"][begin] .. data["Bs"][end],
           B0 -> 1 / (1 + exp(config["selection_strength"] * (nb_effective - 1) *
                              ((nb_effective - 1) * data["cost"]
			       + nb_effective * beta0(B0) * (1-2*config["symmetry_breaking"])/2
                              - (nb_effective - 2) / 2 * B0))); label="Approx. Theory",
           color=:purple, linestyle=:dash)


    # Add legend
    axislegend(ax; position=:lt)

    return fig
end

function calc_timeseries(config::Dict)
    # Unpack variables
    @unpack B_factor, selection_strength, symmetry_breaking, adj_matrix_source, payoff_update_method, time_steps, nb_phases = config

    # Define system
    cost = 0.1
    mutation_rate = 0.0001
    B = cost * B_factor

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Run the model for weak selection strength
    all_populations = time_series(LocalMoranInteraction(NormalFormGame(payoff_matrix(nb_phases,
                                                                                     B,
                                                                                     0.95 *
                                                                                     B,
                                                                                     cost;
                                                                                     symmetry_breaking)),
                                                        interaction_adj_matrix,
                                                        reproduction_adj_matrix,
                                                        selection_strength, mutation_rate),
                                  time_steps, payoff_update_method)

    # Package results
    return @strdict(all_populations, cost, nb_phases, nb_players, mutation_rate, interaction_adj_matrix)
end

function calc_timeseries_statistics(config::Dict)
    # Calculate timeseries
    data = calc_timeseries(config)

    # Unpack variables
    @unpack all_populations, cost, nb_phases, nb_players, mutation_rate, interaction_adj_matrix = data
    @unpack B_factor, symmetry_breaking, selection_strength, adj_matrix_source, payoff_update_method, time_steps, nb_phases = config

    # Extract results
    most_common_game_types = dropdims(mapslices(x -> extract_most_common_game_types(x,
                                                                                    B_factor *
                                                                                    cost,
                                                                                    0.9 *
                                                                                    B_factor *
                                                                                    cost,
                                                                                    cost,
                                                                                    symmetry_breaking,
                                                                                    nb_phases,
                                                                                    interaction_adj_matrix),
                                                all_populations; dims=1); dims=1)
    strategy_parity = dropdims(mapslices(x -> check_all_same_strategy(x, nb_phases),
                                         all_populations; dims=1); dims=1)
    counts = mapslices(x -> extract_counts(x, nb_phases), all_populations; dims=1)
    nb_communicative = map(x -> extract_num_communicative(Vector(x)),
                           eachslice(counts; dims=2))
    fraction_communicative = nb_communicative / nb_players
    order_parameters = dropdims(mapslices(x -> extract_order_parameters(x, nb_phases),
                                          counts; dims=1); dims=1)

    # Package results
    return @strdict(fraction_communicative, order_parameters, most_common_game_types,
                    strategy_parity)
end

function get_frame(data, time_step::Integer; adj_matrix_source::String="well-mixed")
    # Generate graph
    graph = Graphs.SimpleDiGraph(local_moran_interaction.get_adj_matrices(adj_matrix_source)[1])

    # Create colormap
    cooperative_colors = range(colorant"navyblue";
                               stop=colorant"paleturquoise1",
                               length=data["nb_phases"])
    noncooperative_colors = range(colorant"darkred";
                                  stop=colorant"lightsalmon1",
                                  length=data["nb_phases"])
    colormap = [cooperative_colors; noncooperative_colors]

    # Apply colormap
    colors = colormap[data["all_populations"][:,time_step]]

    # Create plot
    fig = graphplot(graph; node_color=colors, layout=Stress(),
                                    arrow_show=false, edge_color=(:black, 0.05),
                                    edge_plottype=:linesegments)

    return (fig,colors)
end

function plot_timeseries(B_factor::Real, selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         payoff_update_method::String="single-update",
                         time_steps::Integer=80_000,
			 nb_phases::Integer=20)
    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps, B_factor,
                      symmetry_breaking, selection_strength, nb_phases)
    data = produce_or_load(calc_timeseries_statistics, config,
                           datadir("raw","timeseries_statistics"))[1]

    # Produce figure
    fig = generate_timeseries_plot(data; time_steps)

    # Save figure
    filename = plotsdir("timeseries", savename(config, "png"))
    save(filename, fig)

    return fig
end

function generate_timeseries_plot(data; time_steps::Integer)

    # Only plot subset of points to prevent large file sizes
    num_samples = 1000
    downsample_ratio = Int(floor((time_steps + 1) / num_samples))

    # Create array of times
    # Note: the populations include the initial data, so we need one more than time-steps
    plot_times = 1:downsample_ratio:(time_steps + 1)

    # Plot fraction communicative
    fig = Figure()
    ax1 = Axis(fig[1, 1];
               title=L"Strong selection $\delta = 0.2$",
               xlabel="Time",
               ylabel="Frequency of communicative strategies",
               limits=(nothing, nothing, -0.05, 1.05))
    strategy_parity_colors = Dict(all_communicative => :lightgrey,
                                  all_noncommunicative => :darkgrey)

    colors_game_type = getindex.(Ref(game_type_colors),
                                 data["most_common_game_types"][plot_times])
    colors = [strategy_parity != mixed ? strategy_parity_colors[strategy_parity] :
              color_game_type
              for (strategy_parity, color_game_type) in
                  zip(data["strategy_parity"], colors_game_type)]
    li1 = lines!(ax1, plot_times, data["fraction_communicative"][plot_times];
                 # Simply decimate color_indices instead of resampling because they're discrete quantities
                 color=colors)

    # Plot order parameter
    ax2 = Axis(fig[1, 1];
               ylabel="Order parameter",
               limits=(nothing, nothing, -0.05, 1.05),
               yaxisposition=:right,
               yticklabelcolor=:orange)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    li2 = lines!(ax2, plot_times, data["order_parameters"][plot_times];
                 color=:orange)

    # Add legend
    axislegend(ax1, [li1, li2], ["frequency_communicative", "Order parameter"];
               position=:rb)

    # Plot histogram of game types
    game_parity_or_type = [strategy_parity != mixed ? strategy_parity : game_type
                           for (strategy_parity, game_type) in
                               zip(data["strategy_parity"], data["most_common_game_types"])]
    hist_data = countmap(String.(Symbol.(game_parity_or_type)))
    ax3 = Axis(fig[2, 1];
               title=L"Strong selection $\delta = 0.2$",
               limits=(nothing, nothing, 0, nothing),
               xticks=(1:length(keys(hist_data)), collect(keys(hist_data))))
    barplot!(ax3, collect(values(hist_data)))

    return fig
end

function plot_connected_components(
                              adj_matrix_source::String="well-mixed",
			      )
    # Generate graph
    graph = Graphs.SimpleDiGraph(local_moran_interaction.get_adj_matrices(adj_matrix_source)[1])
    conn_comp = strongly_connected_components(graph)

    # Calculate connected components
    conn_comp_index = zeros(Int32, nv(graph))
    for (conn_comp_indx, elem) in pairs(conn_comp)
        for indx in elem
            conn_comp_index[indx] = conn_comp_indx
        end
    end
    num_conn_comp = length(conn_comp)

    # Plot graph
    fig = graphplot(graph; layout=Stress(), edge_color=(:black, 0.05),
        node_color=distinguishable_colors(num_conn_comp)[conn_comp_index],
        edge_plottype=:linesegments)

    return fig
end

function load_all_timeseries(time_steps::Integer=80_000)
    # Load dataframe
    df_raw = collect_results(datadir("raw","timeseries"); rinclude = [Regex("time_steps=$time_steps[._]")])

    # Add path
    df = transform(df_raw, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)

    return df
end

function load_all_timeseries_statistics(time_steps::Integer=80_000)
    # Load dataframe
    df_raw = collect_results(datadir("raw","timeseries_statistics"); rinclude = [Regex("time_steps=$time_steps[._]")])

    # Check if GameTypes were reconstructed from JLD2
    if typeof(df_raw.most_common_game_types) == Vector{Union{Missing, Vector{JLD2.ReconstructedPrimitive{:GameType, UInt32}}}}
        # Convert JLD2 Reconstruct variables back to GameType enums
        transform!(df_raw, :most_common_game_types => (x -> map(y -> map(z -> GameType(z.val), y), x)) => :most_common_game_types)
        transform!(df_raw, :strategy_parity => (x -> map(y -> map(z -> StrategyParity(z.val), y), x)) => :strategy_parity)
    end

    # Add path
    df = transform(df_raw, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)

    return df
end

function generate_game_type_distribution_vs_asymmetry_plot(df;
                                                  B_factor::Real,
                                                  selection_strength::Real,
                                                  adj_matrix_source::String="well-mixed",
                                                  )
    # Select subset of dataframe
    df_all_asymm = @rsubset(df, :selection_strength == selection_strength, :adj_matrix_source == adj_matrix_source, :factor == B_factor)
    game_types = proportionmap.(df_all_asymm.most_common_game_types)

    # Generate plot
    fig = Figure()
    CairoMakie.Axis(fig[1, 1], xlabel = L"Asymmetry $\alpha$", title = "Proportion of Game Types by Asymmetry")
    for (indx, row) in enumerate(game_types)
        asymm = df_all_asymm.symmetry_breaking[indx]
        barplot!(repeat([asymm], length(row)), collect(values(row)), stack = repeat([indx], length(row)), color = getindex.(Ref(local_moran_interaction.game_type_colors), collect(keys(row))), width=0.2)
    end

    return fig
end

function plot_graph_evolution(B_factor::Real, selection_strength::Real,
                              symmetry_breaking::Real,
                              adj_matrix_source::String="well-mixed",
                              payoff_update_method::String="single-update",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20)
    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps, B_factor,
                      selection_strength, symmetry_breaking, nb_phases)
    data, _ = produce_or_load(calc_timeseries, config, datadir("raw","timeseries"))

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)

    # Create animation
    animation = generate_graph_evolution(data, graph)

    # Save animation
    filename = plotsdir("animations", savename(config, "gif"))
    save(filename, fig)

end
function generate_graph_evolution(data, graph::Graphs.SimpleGraphs.AbstractSimpleGraph)
    ## Remove self edges for display
    #self_loops = Iterators.flatten(simplecycles_limited_length(graph,1))
    #rem_edge!.(Ref(graph), self_loops, self_loops)

    # Create colormap
    cooperative_colors = range(colorant"navyblue";
                               stop=colorant"paleturquoise1",
                               length=data["nb_phases"])
    noncooperative_colors = range(colorant"darkred";
                                  stop=colorant"lightsalmon1",
                                  length=data["nb_phases"])
    colormap = [cooperative_colors; noncooperative_colors]

    # Apply colormap
    colors = colormap[data["all_populations"]]

    # Set animation parameters
    num_times = size(data["all_populations"], 2)
    total_time_s = 10
    framerate_Hz = 30
    time_stride = round(num_times / (framerate_Hz * total_time_s))

    # Generate animation
    time = Observable(1)

    # Calculate correlation
    correlation = cor(transpose(data["all_populations"]))

    ## Create custom layout based on correlation
    #quantile_cutoff = 0.75
    #correlation_cutoff = quantile!(filter(!isnan, correlation), quantile_cutoff)
    #correlation_mask = correlation .> quantile_cutoff
    #layout = spring(correlation_mask)
    layout = Stress()

    # Create plot
    color_observable = @lift(colors[:, $time])
    fig, ax, graph_plot = graphplot(graph; node_color=color_observable, layout=layout,
                                    arrow_show=false, edge_color=(:black, 0.05),
                                    edge_plottype=:linesegments)

    recording = CairoMakie.Makie.Record(fig, 1:time_stride:num_times; framerate=framerate_Hz) do t
        return time[] = t
    end
    return recording
end

function save_heatmap()
    fig, ax, hm = heatmap(payoff_matrix(20, 1, 0.5, 0.1))

    # Add colorbar
    Colorbar(fig[:, end + 1], hm)

    # Save figure
    filename = plotsdir("heatmap.png")
    save(filename, fig)

    return fig
end

function plot_payoff_regions()
    @variables B_on_c beta_on_c
    payoff_mat = payoff_matrix(1, B_on_c, beta_on_c, 1)

    function inequality_to_hrep(inequality::SymbolicUtils.BasicSymbolic{Bool})
        args = arguments(inequality)
        if operation(inequality) in [<, <=]
            less_than = args[1] - args[2]
        elseif operation(inequality) in [>, >=]
            less_than = args[2] - args[1]
        else
            throw(ArgumentError("Not an inequality"))
        end
        xcoeff = Symbolics.coeff(less_than, B_on_c)
        ycoeff = Symbolics.coeff(less_than, beta_on_c)
        # constant = Symbolics.coeff(less_than) # Note: this seems give the
        # wrong answer if the less_than expression is a negative monomial
        # in one of the variables (e.g. -beta_on_c); specifically, while it
        # should return 0 as the constant term, it returns the monomal
        # (e.g. -beta_on_c)
        # As a workaround, subtract the xcoeff and ycoeff to get the
        # remainder which is equal to the constant as long as the less_than
        # expression is of the form m*B_on_c + n*beta_on_c + p
        constant = (less_than - xcoeff * B_on_c - ycoeff * beta_on_c).val
        return HalfSpace([xcoeff, ycoeff], -constant)
    end

    # Define image borders otherwise the unbounded regions will be truncated to the smallest possible subset
    image_borders = [beta_on_c > 0, B_on_c > 0, beta_on_c < 4, B_on_c < 4]

    # Plot regions
    fig = Figure()
    ax = Axis(fig[1, 1];
              xlabel=L"Benefit of mutual communication $B(\delta \phi)/c$",
              ylabel=L"Benefit of unilateral communication $\beta(\delta \phi)/c$",
              limits=(0, 4, 0, 4))

    for i in 1:2
        R = payoff_mat[1, 1] # Reward
        S = payoff_mat[1, 2] # Sucker's payoff
        T = payoff_mat[2, 1] # Temptation
        P = payoff_mat[2, 2] # Punishment

        for game_type in instances(GameType)
            poly = polyhedron(intersect(inequality_to_hrep.(map(x -> x.val,
                                                                [game_type_inequalities(R,
                                                                                        S,
                                                                                        T,
                                                                                        P)[game_type]...,
                                                                 image_borders...]))...))
            if Polyhedra.volume(poly) == 0
                continue
            end
            mesh!(ax, Polyhedra.Mesh{2}(poly); color=game_type_colors[game_type])
        end

        # Swap C <-> N
        swap_strategies!(payoff_mat)
    end

    # Plot guide lines
    ablines!(ax, 0, 1; color=:grey)
    ablines!(ax, -1, 1; color=:grey)
    hlines!(ax, [1]; color=:grey)
    vlines!(ax, [1]; color=:grey)

    # Add legend
    axislegend(ax,
               [PolyElement(; color=game_type_colors[game_type], strokewidth=1,
                            strokecolor=:grey) for game_type in instances(GameType)],
               [game_type_full_names[game_type] for game_type in instances(GameType)];
               position=:lt)

    # Save figure
    filename = plotsdir("payoff_regions.png")
    save(filename, fig)

    return fig
end

function calc_coalescence_times(adj_matrix_source::String="well-mixed")
    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)
    return calc_coalescence_times(graph)
end

function calc_coalescence_times(graph::Graphs.SimpleGraphs.AbstractSimpleGraph)
    # Take cartesian product
    cart_prod = cartesian_product(graph, graph)

    # Calculate Degree vector and Laplacian matrix
    degree_cartProd = degree(cart_prod)
    laplacian_cartProd = laplacian_matrix(cart_prod)

    # Get the off-diagonal indices
    indices_offDiag = [i
                       for i in CartesianIndices((1:size(graph, 1), 1:size(graph, 1)))
                       if i[1] != i[2]]
    indicesLinear_offDiag = LinearIndices((1:size(graph, 1), 1:size(graph, 1)))[indices_offDiag]

    # Get the off-diagonal parts of Degree vector and Laplacian matrix
    degree_cartProd_offDiag = degree_cartProd[indicesLinear_offDiag]
    laplacian_cartProd_offDiag = laplacian_cartProd[indicesLinear_offDiag,
                                                    indicesLinear_offDiag]

    # Solve for coalescence time matrix
    #coalescence_matrix_offDiag = laplacian_cartProd_offDiag \ degree_cartProd_offDiag
    (pfunc, h) = cmg_preconditioner_lap(Float64.(laplacian_cartProd_offDiag))
    coalescence_matrix_offDiag = pfunc(Float64.(degree_cartProd_offDiag))
    coalescence_matrix = sparse([v[1] for v in indices_offDiag],
                                [v[2] for v in indices_offDiag],
                                Vector(coalescence_matrix_offDiag))

    return coalescence_matrix
end

function get_connectome()
    connectome_and_muscles = get_connectome_labelled()["connectome"]
    # Remove connections to muscles
    connectome = connectome_and_muscles[:, [1:20..., 51:322..., 438:445...]]
    # Replace "Missing" data with zeros
    replace!(connectome, missing => 0)
    return connectome
end

function get_connectome_labelled()
    connectome_and_muscles_with_labels = read(dataset("connectome-cook"), Matrix)
    row_labels = connectome_and_muscles_with_labels[:, 1]
    col_labels = connectome_and_muscles_with_labels[1, :]
    connectome_and_muscles = connectome_and_muscles_with_labels[2:end, 2:end]
    results = Dict("connectome" => connectome_and_muscles, "row_labels" => row_labels,
                   "col_labels" => col_labels)
    return results
end

function export_graph_nodes_edges(time_step::Integer;
                              B_factor::Real, selection_strength::Real,
                              symmetry_breaking::Real,
                              adj_matrix_source::String="well-mixed",
                              payoff_update_method::String="single-update",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20)

    # Generate configuration
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps, B_factor,
                      selection_strength, symmetry_breaking, nb_phases)

    # Get data
    data = wload(datadir("raw","timeseries",savename(config,"jld2")))

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)

    # Calculate nodes and edges
    nodes_df, edges = generate_nodes_edges(graph, data, time_step)

    # Add index number
    insertcols!(nodes_df, 1, :index => 1:nrow(nodes_df))

    # Add time_step to config dictionary
    config["time_step"] = time_step

    # Write out vertices and edges
    CSV.write(datadir("processed",savename("vertices",config,"csv")), nodes_df)
    CSV.write(datadir("processed","edges.csv"), edges)

    return edges
end

function generate_nodes_edges(graph::AbstractGraph, data::Dict, time_step::Integer)
    # Get vertex labels
    vertex_labels = data["all_populations"][:,time_step]

    # Get vertex coordinates
    vertex_coordinates = stress(graph)

    # Split points into coordinates
    xs = (p -> p[1]).(vertex_coordinates)
    ys = (p -> p[2]).(vertex_coordinates)

    # Put into DataFrame
    df = DataFrame(x=xs, y=ys, strategyIndex=vertex_labels)

    # Get edges
    graph_edges = collect(edges(graph))

    return df, graph_edges
end

end

using .local_moran_interaction
