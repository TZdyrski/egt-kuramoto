module local_moran_interaction

using GameTheory
using LinearAlgebra
using BlockArrays
using SparseArrays
using Plots
using Random
using Statistics
using StatsBase
using Symbolics
using BenchmarkTools
using InteractiveUtils
using LaTeXStrings
using DataToolkit
using DrWatson
using ImplicitEquations
using Graphs

function payoff_matrix_template(benefit_scaling::AbstractMatrix,
		mutual_benefit_synchronous::Real,
		unilateral_benefit_synchronous::Real, cost::Real)
  payoff_matrix = Matrix(mortar([[mutual_benefit_synchronous*benefit_scaling .- cost,
				  unilateral_benefit_synchronous*benefit_scaling];;
				 [unilateral_benefit_synchronous*benefit_scaling .- cost,
				 zeros(eltype(benefit_scaling), size(benefit_scaling))]]) .+ cost);
  return payoff_matrix
end

function payoff_matrix(mutual_benefit_synchronous::Real,
		unilateral_benefit_synchronous::Real, cost::Real,
		nb_phases::Integer)
    phase_dependence = [(1 + cos(2*pi*(phi_j - phi_i)))/2 for phi_i in (0:(nb_phases-1))/nb_phases, phi_j in (0:(nb_phases-1))/nb_phases]
    payoff_matrix = payoff_matrix_template(phase_dependence,
					   mutual_benefit_synchronous,
					   unilateral_benefit_synchronous,
					   cost)
    return payoff_matrix
end;

struct LocalMoranInteraction{N,T1<:Real,S<:Integer,S2<:Integer,T2<:Real,T3<:Real,AT1<:AbstractMatrix{<:T1},AT2<:AbstractMatrix{<:T1},AS1<:AbstractMatrix{S2},AS2<:AbstractMatrix{S2},AS3<:AbstractMatrix{S2},AS4<:AbstractMatrix{S2}}
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
    return LocalMoranInteraction(players, num_actions, g.players[1].payoff_array, transpose(g.players[1].payoff_array),
				 interaction_adj_matrix, transpose(interaction_adj_matrix), reproduction_adj_matrix, transpose(reproduction_adj_matrix), selection_strength, mutation_rate)
end

function update_weights!(w::Weights{S, TA, V}, new_wts::V) where {S, TA, V}
	w.values = new_wts
	w.sum = sum(w.values)
end

abstract type WorkParams end

struct WorkParamsLoop{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,AT3<:AbstractWeights,AS2<:AbstractWeights} <: WorkParams
    neighbor_idxs::AS1
    payoffs::AT1
    fitnesses::AT2
    weights_float::AT3
    weights_int::AS2
end

function WorkParamsLoop(lmi::LocalMoranInteraction{N}) where N
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    return WorkParamsLoop(neighbor_idxs, payoffs, fitnesses, weights_float, weights_int)
end

struct WorkParamsSingleUpdate{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,AT3<:AbstractWeights,AS2<:AbstractWeights,
			      AT4<:AbstractMatrix,AT5<:AbstractMatrix,AT6<:AbstractMatrix,AT7<:AbstractMatrix} <: WorkParams
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

function WorkParamsSingleUpdate(lmi::LocalMoranInteraction{N}, initial_actions::AbstractVector) where N
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    payoffs_transpose = Matrix{Float64}(undef, 1, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    payoffs_rowplayer_per_strategy = lmi.payoff_matrix[initial_actions, :]
    payoffs_colplayer_per_strategy = lmi.payoff_matrix_transpose[initial_actions, :]
    payoffs_player_pairwise_transpose = payoffs_colplayer_per_strategy[:, initial_actions] .* lmi.interaction_adj_matrix
    return WorkParamsSingleUpdate(neighbor_idxs, payoffs, payoffs_transpose, fitnesses, weights_float, weights_int,payoffs_rowplayer_per_strategy,payoffs_colplayer_per_strategy,payoffs_player_pairwise_transpose)
end

struct WorkParamsMatrix{AS1<:AbstractVector,AT1<:AbstractVector,AT2<:AbstractVector,AT3<:AbstractWeights,AS2<:AbstractWeights,
			AS3<:AbstractMatrix,AT4<:AbstractMatrix,AS4<:AbstractVector} <: WorkParams
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

function WorkParamsMatrix(lmi::LocalMoranInteraction{N}) where N
    neighbor_idxs = Vector{Int64}(undef, N)
    payoffs = Vector{Float64}(undef, N)
    fitnesses = Vector{Float64}(undef, N)
    weights_float = Weights(ones(Float64, N))
    weights_int = Weights(ones(Int64, N))
    actions_matrix = Matrix{Int64}(undef,N,lmi.num_actions)
    actions_matrix_T = Matrix{Int64}(undef,lmi.num_actions,N)
    work_array_1 = Matrix{Float64}(undef,N,lmi.num_actions)
    work_array_2 = Matrix{Float64}(undef,N,N)
    work_array_3 = Matrix{Float64}(undef,N,N)
    ones_matrix = ones(Int64,N)
    return WorkParamsMatrix(neighbor_idxs, payoffs, fitnesses, weights_float, weights_int, actions_matrix, actions_matrix_T, work_array_1, work_array_2, work_array_3, ones_matrix)
end

function calc_payoffs!(aux::WorkParamsLoop,
		actions::AbstractVector,
		lmi::LocalMoranInteraction{N}) where N
    @inbounds for focal_idx = 1:length(actions)
        focal_strategy = actions[focal_idx]
        total = 0.0
        @inbounds for neighbor_idx = 1:length(actions)
            neighbor_strategy = actions[neighbor_idx]
            total += lmi.payoff_matrix_transpose[neighbor_strategy, focal_strategy] * lmi.interaction_adj_matrix[neighbor_idx, focal_idx]
        end
        aux.payoffs[focal_idx] = total
    end
end

function calc_payoffs!(aux::WorkParamsSingleUpdate,
		actions::AbstractVector,
		lmi::LocalMoranInteraction{N}) where N
    # Calculate payoffs
    sum!(aux.payoffs_transpose, aux.payoffs_player_pairwise_transpose)
    transpose!(aux.payoffs, aux.payoffs_transpose)
end

function update_aux!(aux::WorkParams,
		new_strategy::Integer,
		death_idx::Integer,
		lmi::LocalMoranInteraction{N}) where N
end

function update_aux!(aux::WorkParamsSingleUpdate,
		new_strategy::Integer,
		death_idx::Integer,
		lmi::LocalMoranInteraction{N}) where N
    # Update various matrices
    @. aux.payoffs_rowplayer_per_strategy[death_idx, :] = @view lmi.payoff_matrix_transpose[:, new_strategy]
    @. aux.payoffs_colplayer_per_strategy[death_idx, :] = @view lmi.payoff_matrix[:, new_strategy]
    @. aux.payoffs_player_pairwise_transpose[:, death_idx] = @view aux.payoffs_colplayer_per_strategy[:, new_strategy]
    for idx = 1:N
        aux.payoffs_player_pairwise_transpose[idx, death_idx] *= lmi.interaction_adj_matrix[idx, death_idx]
    end
    @. aux.payoffs_player_pairwise_transpose[death_idx, :] = @view aux.payoffs_rowplayer_per_strategy[:, new_strategy]
    for idx = 1:N
        aux.payoffs_player_pairwise_transpose[death_idx, idx] *= lmi.interaction_adj_matrix[death_idx,idx]
    end
end


 
function calc_payoffs!(aux::WorkParamsMatrix,
		actions::AbstractVector,
		lmi::LocalMoranInteraction{N}) where N
    @inbounds for col = 1:lmi.num_actions, row = 1:N
        aux.actions_matrix[row,col] = 0
    end
    @inbounds for player_idx = 1:N
       strategy_idx = actions[player_idx]
        aux.actions_matrix[player_idx,strategy_idx] = 1
    end
    transpose!(aux.actions_matrix_T, aux.actions_matrix)
    mul!(aux.work_array_1, aux.actions_matrix, lmi.players[1].payoff_array)
    mul!(aux.work_array_2, aux.work_array_1, aux.actions_matrix_T)
    aux.work_array_3 .= aux.work_array_2 .* lmi.interaction_adj_matrix_transpose
    mul!(aux.payoffs, aux.work_array_3, aux.ones_matrix)
end

function play!(actions::AbstractVector,
    rng::AbstractRNG,
    lmi::LocalMoranInteraction{N},
    aux::WorkParams) where N
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
    @. aux.neighbor_idxs = @view lmi.reproduction_adj_matrix_transpose[:,focal_idx]
    # Choose death node
    update_weights!(aux.weights_int, aux.neighbor_idxs)
    death_idx = sample(rng, 1:N, aux.weights_int)
    
    # Check for mutation
    if rand(rng) <= lmi.epsilon
	new_strategy = rand(rng, 1:lmi.num_actions)
    else
	new_strategy = actions[focal_idx]
    end
    
    # Apply spatial Moran process
    actions[death_idx] = new_strategy

    # Update aux variables
    update_aux!(aux, new_strategy, death_idx, lmi)

    return nothing
end

function count!(counts::Vector{N}, new_values::Vector{N}) where N
    for (i, _) in enumerate(counts)
        counts[i] += count(==(i), new_values)
    end
    return counts
end

function time_series(lmi::LocalMoranInteraction{N}, ts_length::Integer, payoff_counting_method::String="single-update", seed::Integer=12345) where N
    rng = Xoshiro(seed)
    actions = rand(rng, 1:lmi.num_actions, N)
    if payoff_counting_method == "loop"
        aux = WorkParamsLoop(lmi)
    elseif payoff_counting_method == "single-update"
        aux = WorkParamsSingleUpdate(lmi, actions)
    elseif payoff_counting_method == "matrix"
        aux = WorkParamsMatrix(lmi)
    else
        throw(ArgumentError("payoff_counting_method must be a string in set [\"loop\", \"single-update\", \"matrix\"]"))
    end
    out = Matrix{Int}(undef, N, ts_length+1)
    for i in 1:N
        out[i,1] = actions[i]
    end
    for t in 1:ts_length
	play!(actions, rng, lmi, aux)
        out[:,t+1] = actions
    end
    return out
end

function cumulative(lmi::LocalMoranInteraction{N}, ts_length::Integer, payoff_counting_method::String="single-update", seed::Integer=12345) where N
    rng = Xoshiro(seed)
    actions = rand(rng, 1:lmi.num_actions, N)
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
    B(phi) = B_0*(1 + cos(phi))/2
    beta(phi) = beta_0*(1+cos(phi))/2

    # Define fixation probabilities
    rho_CC(phi) = 1/(1+sum([exp(k*(k+1-n)*delta/(n-1)*(B(phi)-B_0)) for k in 1:(n-1)]))
    rho_NC(phi) = 1/(1+sum([exp(delta/(n-1)*
                    (j^2*(beta(phi)-1/2*B_0) + j*(1/2*B_0-(n-1)*beta(phi)+(n-1)*c)))
                for j in 1:(n-1)]))
    omega = 1/exp(delta*((n-1)*c-(n-2)/2*B_0))
    rho_CN(phi) = rho_NC(phi)/omega

    # Define matrices
    B_1 = [i != j ?
        rho_CC(floor(d/2) - abs(floor(d/2)-abs(i-j))) :
        1 - sum([rho_CC(floor(d/2) - abs(floor(d/2) - abs(delta_phi_prime))) for delta_phi_prime in 1:(d-1)])
        for i in 1:d, j in 1:d]
    B_2 = [rho_CN(floor(d/2) - abs(floor(d/2)-abs(i-j))) for i in 1:d, j in 1:d]
    B_3 = [rho_NC(floor(d/2) - abs(floor(d/2)-abs(i-j))) for i in 1:d, j in 1:d]
    B_4 = [1/n for i in 1:d, j in 1:d]

    # Define full Markov transition matrix
    M = mortar([[B_1,B_3] [B_2,B_4]])

    # Define stationary distribution
    # Tripp (2024) shows s_1, s_2 independent of q, so f
    q = 1
    s_1 = sum([M[r,q] for r in (d+1):(2d)])/(d*sum([M[q,r] + M[r,q] for r in (d+1):(2d)]))
    s_2 = sum([M[q,r] for r in (d+1):(2d)])/(d*sum([M[q,r] + M[r,q] for r in (d+1):(2d)]))
    s = hcat([s_1*ones(d), s_2*ones(d)]);

    ## Check that the explicitly calculated stationary distribution matches the eigenvector
    #eigen = eigvals(M)[1]
    #if s != eigen
    #    throw(ErrorException(\"Explicitly calculated stationary distribution differs from eigenvector\"))

    ## Check that fraction of communicative simplifies correctly
    #if not isequal(simplify(d*s_1, threaded=true), 1/(1+1/omega))
    #    throw(ErrorException(\"Cooperative fraction does not simplify to known correct answer 1/(1+1/omega)))
end

function extract_num_communicative(players_per_strategy::Vector{<:Integer})
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
    fraction_communicative = nb_communicative./(time_steps*nb_players)
    return fraction_communicative
end

@enum GameType all_communicative all_noncommunicative prisoners_dilemma snowdrift coordination mutualism deadlock unknown

function combine_communicative_noncommunicative(
        players_per_strategy::Vector{<:Integer},
        nb_phases::Integer)
    # Combine communicatative (value < nb_phases)
    # and non-communicatative (values >= nb_phases)
    combine_com_noncom_matrix = hcat(Matrix{Int}(I,nb_phases,nb_phases),Matrix{Int}(I,nb_phases,nb_phases))
    phase_idxs =  combine_com_noncom_matrix * players_per_strategy
    return phase_idxs
end

function extract_most_common_game_types(
        strategies_per_player::Vector{<:Integer},
        mutual_benefit_synchronous::Real,
        unilateral_benefit_synchronous::Real,
        cost::Real,
        nb_phases::Integer,
        nb_strategies::Integer,
        interaction_adj_matrix::AbstractMatrix{<:Integer})

    nb_players = length(strategies_per_player)
    players_per_strategy = extract_counts(strategies_per_player, nb_strategies)
    # Check if all players were communicative/noncommunicative
    nb_communicative = extract_num_communicative(players_per_strategy)
    if nb_communicative == 0
        return all_noncommunicative
    elseif nb_communicative == nb_players
        return all_communicative
    end

    players_per_phase = mod1.(strategies_per_player, nb_phases)

    # Define payoff response submatrices
    payoff = payoff_matrix(mutual_benefit_synchronous, unilateral_benefit_synchronous, cost, nb_phases)
    reward_submatrix = payoff[1:nb_phases, 1:nb_phases]
    sucker_submatrix = payoff[1:nb_phases, (nb_phases+1):(2*nb_phases)]
    temptation_submatrix = payoff[(nb_phases+1):(2*nb_phases), 1:nb_phases]
    punishment_submatrix = payoff[(nb_phases+1):(2*nb_phases), (nb_phases+1):(2*nb_phases)]

    # Define mask matrix for each game type
    prisoners_dilemma_matrix = (temptation_submatrix .>= reward_submatrix) .&
        (reward_submatrix .>= punishment_submatrix) .&
        (punishment_submatrix .>= sucker_submatrix)
    snowdrift_matrix = (temptation_submatrix .>= reward_submatrix) .&
        (reward_submatrix .>= sucker_submatrix) .&
        (sucker_submatrix .>= punishment_submatrix)
    coordination_matrix = (reward_submatrix .>= temptation_submatrix) .&
        (temptation_submatrix .>= sucker_submatrix) .&
        (sucker_submatrix .>= punishment_submatrix)
    mutualism_matrix = (reward_submatrix .>= temptation_submatrix) .&
        (temptation_submatrix .>= punishment_submatrix) .&
        (punishment_submatrix .>= sucker_submatrix)
    deadlock_matrix = (temptation_submatrix .>= punishment_submatrix) .&
        (punishment_submatrix .>= reward_submatrix) .&
        (reward_submatrix .>= sucker_submatrix)

    # Define game type per strategy pair
    game_types = fill(unknown, nb_phases, nb_phases)
    for idx in eachindex(game_types)
        if prisoners_dilemma_matrix[idx]
	    game_type = prisoners_dilemma
        elseif snowdrift_matrix[idx]
            game_type = snowdrift
        elseif coordination_matrix[idx]
            game_type = coordination
        elseif mutualism_matrix[idx]
            game_type = mutualism
        elseif deadlock_matrix[idx]
            game_type = deadlock
	else
	    throw(Exception("Unknown game types"))
        end
	game_types[idx] = game_type
    end

    # Count game types
    game_counts = Dict{GameType, Integer}(prisoners_dilemma => 0,
                   snowdrift => 0,
                   coordination => 0,
                   mutualism => 0,
                   deadlock => 0)
    for (cart_idx,value) in pairs(interaction_adj_matrix)
        if value == 0
            continue
	end
        (row, col) = Tuple(cart_idx)
	row_phase = players_per_phase[row]
	col_phase = players_per_phase[col]
	game_type = game_types[row_phase,col_phase]
        game_counts[game_type] += value
    end

    # Find most common game type
    most_common_game_type = findmax(game_counts)[2]

    return most_common_game_type
end;

function extract_phases(players_per_phase::Vector{<:Integer}, nb_phases::Integer)
    # Convert phase idxs to phases
    phases = 2*pi/nb_phases.*players_per_phase
    return phases
end

function extract_order_parameters(players_per_strategy::Vector{<:Integer}, nb_phases::Integer)
    # Combine communicative and noncommunicative
    phase_indxs = combine_communicative_noncommunicative(players_per_strategy, nb_phases)
    # Convert strategy idxs to phases
    phases = extract_phases(phase_indxs, nb_phases)
    order_parameters = abs(mean(exp.(im*phases)))
    return order_parameters
end;

function extract_counts(strategies_per_player::Vector{<:Integer}, nb_strategies::Integer)
    counts = zeros(Int, nb_strategies)
    count!(counts, strategies_per_player)
    return counts
end

function calc_cumulative(config::Dict)
	# Unpack values
	@unpack selection_strength, adj_matrix_source, payoff_update_method, time_steps = config

	# Define system
	cost = 0.1
	nb_phases = 20
	nb_strategies = nb_phases*2
	mutation_rate = 0.0001

	# Define interaction graph without loops
	if adj_matrix_source == "well-mixed"
	    interaction_adj_matrix = ones(Int64,20,20) - I
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
				*"\"c-elegans\", \"c-elegans-unweighted\", "
				*"\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\"]"))
        end
	# Define reproduction graph with loops
	# Specify number of players
	nb_players = size(interaction_adj_matrix)[1]

        # Define Bs on which to run
        B_crit = 2 * cost * (nb_players-1) / (nb_players-2)
        nb_Bs = 11
        step_size_Bs = 0.04
        Bs = B_crit .+ ((0:(nb_Bs-1)) .- (nb_Bs-1)/2) .* step_size_Bs

        # Run the model for weak selection strength
        @time cumulative_populations = [cumulative(
                LocalMoranInteraction(NormalFormGame(payoff_matrix(B, 0.95*B, cost, nb_phases)), interaction_adj_matrix, reproduction_adj_matrix, selection_strength, mutation_rate),
                time_steps, payoff_update_method) for B in Bs]

        # Plot fraction communcative
        nb_communicative = [extract_num_communicative(final_population) for final_population in cumulative_populations]
        fraction_communicative = nb_communicative./(time_steps*nb_players)

	# Package results
	return @strdict(Bs,nb_players,cost,fraction_communicative)
end

function plot_cumulative(selection_strength::Real, adj_matrix_source::String="well-mixed", payoff_update_method::String="single-update",time_steps::Integer=2_000_000)

	# Load results
	config = @strdict(adj_matrix_source,payoff_update_method,time_steps,selection_strength)
	data, _ = produce_or_load(calc_cumulative, config, datadir("cumulative"))

	# Disable printing plots to screen
	ENV["GKSwstype"]="nul"

	plt = scatter(data["Bs"], data["fraction_communicative"], title=L"Selection $\delta = $"*string(round(selection_strength,sigdigits=2)),
		label="Simulation", xlabel=L"Maximum benefit of mutual communication, $B(0)$", ylabel="Frequency of communicative strategies", ylims=(0,1))
	plot!(plt, B0 -> 1/(1+exp(selection_strength*(data["nb_players"]-1)*((data["nb_players"]-1)*data["cost"]-(data["nb_players"]-2)/2*B0))), label="Theory")

	# Save figure
	filename = plotsdir("cumulative", savename(config))
	png(plt, filename)

	return plt
end

function calc_timeseries(config::Dict)
	# Unpack variables
	@unpack B_factor, selection_strength, adj_matrix_source, payoff_update_method, time_steps = config

	# Define system
	cost = 0.1
	nb_phases = 20
	nb_strategies = nb_phases*2
	mutation_rate = 0.0001
	B = cost*B_factor

	# Define interaction graph without loops
	if adj_matrix_source == "well-mixed"
	    interaction_adj_matrix = ones(Int64,20,20) - I
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
				*"\"c-elegans\", \"c-elegans-unweighted\", "
				*"\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\"]"))
        end
	# Define reproduction graph with loops
	# Specify number of players
	nb_players = size(interaction_adj_matrix)[1]

        # Run the model for weak selection strength
        all_populations = time_series(
                LocalMoranInteraction(NormalFormGame(payoff_matrix(B, 0.95*B, cost, nb_phases)), interaction_adj_matrix, reproduction_adj_matrix, selection_strength, mutation_rate),
                time_steps, payoff_update_method)

        # Extract results
        most_common_game_types = dropdims(mapslices(x -> extract_most_common_game_types(x, B,
                0.9*B, cost, nb_phases, nb_strategies, interaction_adj_matrix), all_populations, dims=1), dims=1)
        counts = mapslices(x -> extract_counts(x, nb_strategies), all_populations, dims=1)
        nb_communicative = map(x -> extract_num_communicative(Vector(x)), eachslice(counts, dims=2))
        fraction_communicative = nb_communicative/nb_players
        order_parameters = dropdims(mapslices(x -> extract_order_parameters(x, nb_phases), counts, dims=1), dims=1)

	# Package results
	return @strdict(fraction_communicative,order_parameters,most_common_game_types)
end

function plot_timeseries(B_factor::Real, selection_strength::Real, adj_matrix_source::String="well-mixed", payoff_update_method::String="single-update",time_steps::Integer=80_000)
	# Load results
	config = @strdict(adj_matrix_source,payoff_update_method,time_steps,B_factor,selection_strength)
	data, _ = produce_or_load(calc_timeseries, config, datadir("timeseries"))

        # Create array of times
        # Note: the populations include the initial data, so we need one more than time-steps
        times = 1:(time_steps+1)

	# Disable printing plots to screen
	ENV["GKSwstype"]="nul"

	# Only plot subset of points to prevent large file sizes
	plot_times = 1:Int(floor(length(times)/1000)):length(times)
        # Plot fraction communcative
	colors = Dict(all_noncommunicative => :darkgrey,
		      all_communicative => :lightgrey,
		      mutualism => :green,
		      coordination => :blue,
		      snowdrift => :yellow,
		      prisoners_dilemma => :red)
	color_values = get.(Ref(colors), data["most_common_game_types"][plot_times], :purple)
        plt1 = plot(times[plot_times], data["fraction_communicative"][plot_times], color=color_values, label="frequency_communicative", title=raw"Strong selection $\delta = 0.2$", xlabel=raw"Time", ylabel="Frequency of communicative strategies", ylims=(0,1))
        plot!(plt1, times[plot_times], data["order_parameters"][plot_times], label="Order parameter")
        # Plot histogram of game types
        plt2 = bar(countmap(String.(Symbol.(data["most_common_game_types"]))), title=raw"Strong selection $\delta = 0.2$", legend=false)

	# Combine plots
        plt = plot(plt1, plt2, layout=(2,1))

	# Save figure
	filename = plotsdir("timeseries", savename(config))
	png(plt, filename)

	return plt
end

function save_heatmap()
	plt = heatmap(payoff_matrix(1, 0.5, 0.1, 20))

	filename = plotsdir("heatmap")
	png(plt, filename)

	return plt
end

function plot_payoff_regions()
    @variables B_on_c beta_on_c
    payoff_matrix = payoff_matrix_template(ones(Int64,1,1), B_on_c, beta_on_c, 1)

    R = build_function(payoff_matrix[1,1], B_on_c, beta_on_c, expression=Val{false}) # Reward
    S = build_function(payoff_matrix[1,2], B_on_c, beta_on_c, expression=Val{false}) # Sucker's payoff
    T = build_function(payoff_matrix[2,1], B_on_c, beta_on_c, expression=Val{false}) # Temptation
    P = build_function(payoff_matrix[2,2], B_on_c, beta_on_c, expression=Val{false}) # Punishment

    # Disable printing plots to screen
    ENV["GKSwstype"]="nul"

    # Plot regions
    plt = plot(Ge(T, R) & Ge(R, P) & Ge(P, S),
	       fillcolor = :red, label = "Prisoner's Dilemma",
	       xlims=(0,4), ylims=(0,4),
	       xlabel=L"Benefit of mutual communication $B(\delta \phi)/c$",
	       ylabel=L"Benefit of unilateral communication $\beta(\delta \phi)/c$")
    plot!(plt, Ge(T, R) & Ge(R, S) & Ge(S, P),
          fillcolor = :yellow, label = "Snowdrift")
    plot!(plt, Ge(R, T) & Ge(T, S) & Ge(S, P),
          fillcolor = :green, label = "Coordination")
    plot!(plt, Ge(R, T) & Ge(T, P) & Ge(P, S),
          fillcolor = :blue, label = "Mutualism")
    plot!(plt, Ge(T, P) & Ge(P, R) & Ge(R, S),
          fillcolor = :purple, label = "Deadlock")

    # Plot guide lines
    Plots.abline!(plt, 1, 0, color=:grey, label = nothing)
    Plots.abline!(plt, 1, -1, color=:grey, label = nothing)
    Plots.hline!(plt, [1], color=:grey, label = nothing)
    Plots.vline!(plt, [1], color=:grey, label = nothing)

    # Add legend
    plot!(legend = true)

    filename = plotsdir("payoff_regions")
    png(plt, filename)

    return plt
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
	row_labels = connectome_and_muscles_with_labels[:,1]
	col_labels = connectome_and_muscles_with_labels[1,:]
	connectome_and_muscles = connectome_and_muscles_with_labels[2:end,2:end]
	results = Dict("connectome" => connectome_and_muscles, "row_labels" =>
		       row_labels, "col_labels" => col_labels)
	return results
end

end

using .local_moran_interaction
