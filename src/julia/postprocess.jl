# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson
using Graphs
using Random
using SimpleWeightedGraphs
using DataFramesMeta
using Memoize
using DimensionalData
using DimensionalData.Lookups
using NetworkLayout
using CSV
using Statistics
using PlotUtils
using SplitApplyCombine
using JLD2
using YAXArrays
using NetCDF
using IGraphs
using StatsBase

include("game_taxonomy.jl")
include("utils.jl")
include("moran.jl")

# Declare zero type for sparse matrices
Base.zero(::Union{Missing,Type{GameType}}) = missing

function count_games(strategies_per_player::AbstractVector{<:Integer};
                                        game_types::AbstractMatrix{Union{GameType,Missing}},
                                        interaction_adj_matrix::AbstractMatrix{<:Integer},
                                        interaction_adj_matrix_transpose::AbstractMatrix{<:Integer}=transpose(interaction_adj_matrix),
                                        )
    # Initial run with no pre-calculated games

    # All entries are changed
    changes = fill(true, size(strategies_per_player))

    # Create an initial games matrix filled with neutral as a placeholder;
    # this will be replaced with the actual value
    games = sparse([],[],Union{Missing,GameType}[],size(interaction_adj_matrix)...)
    num_games = sum(interaction_adj_matrix)
    game_counts = Dict{Union{GameType,Missing},Integer}(instances(GameType) .=> 0)
    game_counts[zero(Union{Missing,GameType})] = num_games

    # Now actually run the extraction
    return count_games(strategies_per_player, games, game_counts, changes;
					  game_types, interaction_adj_matrix, interaction_adj_matrix_transpose)
end

function count_games(strategies_per_player::AbstractVector{<:Integer},
                                        games::AbstractMatrix{Union{GameType,Missing}},
																				game_counts::Dict{Union{GameType,Missing},Integer},
                                        changes::AbstractVector,
																				;
																				game_types::AbstractMatrix{Union{GameType,Missing}},
                                        interaction_adj_matrix::AbstractMatrix{<:Integer},
                                        interaction_adj_matrix_transpose::AbstractMatrix{<:Integer}=transpose(interaction_adj_matrix),
                                        )

    changed_idxs = findall(!iszero, changes .!== 0)

    # Count game types
    function update_counts(row, col, value, games, game_counts)
        row_phase = strategies_per_player[row]
        col_phase = strategies_per_player[col]
        game_type = game_types[row_phase, col_phase]
        old_game_type = games[row, col]
        game_counts[old_game_type] -= value
        games[row, col] = game_type
        game_counts[game_type] += value
				return games, game_counts
    end

    # Changed nodes are destination nodes
    for col in changed_idxs
      for row = findall(!iszero, interaction_adj_matrix[:, col])
				value = interaction_adj_matrix[row, col]
				games, game_counts = update_counts(row, col, value, games, game_counts)
      end
    end

    # Changed nodes are source nodes
    for row in changed_idxs
      for col = findall(!iszero, interaction_adj_matrix_transpose[:, row])
        value = interaction_adj_matrix_transpose[col, row]
        games, game_counts = update_counts(row, col, value, games, game_counts)
      end
    end

		return games, game_counts
end

function check_all_same_strategy(strategies_per_player::AbstractVector{<:Integer},
																 nb_phases::Integer)
		nb_players = length(strategies_per_player)
		players_per_strategy = extract_counts(strategies_per_player, nb_phases)
		# Check if all players were communicative/noncommunicative
		nb_communicative = extract_num_communicative(players_per_strategy)
		if nb_communicative == 0
				return allNoncommunicative
		elseif nb_communicative == nb_players
				return allCommunicative
		end

		return nothing
end

function check_disconnected_synchronized(
																game_counts::Dict{Union{GameType,Missing},Integer},
																)
		total_games = sum(values(game_counts))
		if get(game_counts,missing,0) == total_games
			# There were no mixed games;
			# however, the all-C/all-N check also shows the population is
			# not entirely communicative or entirely non-communicative
			# This means that the communicative and non-communicative
			# subpopulations are not connected (and therefore do not play any
			# games together)
			return disconnectedSynchronizedPopulations
		end

		return nothing
end

function extract_most_common_game_types(game_counts::Dict{Union{GameType,Missing},Integer})
    # Create game count dictionary without missing
    non_missing_game_counts = copy(game_counts)
    delete!(non_missing_game_counts, missing)
		# Find most common game type
		most_common_game_type = findmax(non_missing_game_counts)[2]

		return most_common_game_type
end

function extract_order_parameters(players_per_strategy::AbstractVector{<:Integer},
                                  nb_phases::Integer)
		function combine_communicative_noncommunicative(players_per_strategy::AbstractVector{<:Integer},
																										nb_phases::Integer)
				# Combine communicatative (value < nb_phases)
				# and non-communicatative (values >= nb_phases)
				combine_com_noncom_matrix = hcat(Matrix{Int}(I, nb_phases, nb_phases),
																				 Matrix{Int}(I, nb_phases, nb_phases))
				phase_idxs = combine_com_noncom_matrix * players_per_strategy
				return phase_idxs
		end

		function extract_phases(players_per_phase::AbstractVector{<:Integer}, nb_phases::Integer)
				# Convert phase idxs to phases
				phases = 2 * pi / nb_phases .* players_per_phase
				return phases
		end

    # Combine communicative and noncommunicative
    phase_counts = combine_communicative_noncommunicative(players_per_strategy, nb_phases)
    # Generate phase for each index
    phases = extract_phases(1:nb_phases, nb_phases)
    order_parameters = abs(sum(phase_counts .* exp.(im * phases)))/sum(phase_counts)
    return order_parameters
end

function extract_counts(strategies_per_player::AbstractVector{<:Integer}, nb_phases::Integer)
    nb_strategies = 2*nb_phases
    counts = zeros(Int, nb_strategies)
    count_actions!(counts, strategies_per_player)
    return counts
end

function generate_communities(graph::AbstractSimpleWeightedGraph, community_algorithm::String;
        covariance_cutoff_fraction::Union{Real,Nothing}=nothing,
        covariance_data::Union{DimArray,Nothing}=nothing,
        walktrap_steps::Union{Integer,Nothing}=nothing,
        community_resolution::Union{Real,Nothing}=nothing,
        community_beta::Union{Real,Nothing}=nothing,
        community_n_iter::Union{Integer,Nothing}=nothing,
        )
    if community_algorithm == "label-propagation"
        # Label propagation
        communities = label_propagation(graph; rng=Xoshiro(12345))[1]
    elseif community_algorithm == "strongly-connected"
        # Strongly connected components
        connected_components = strongly_connected_components(graph)
        communities = Vector{Int64}(undef, nv(graph))
        for (community, idxs) in pairs(connected_components)
            for idx in idxs
                communities[idx] = community
            end
        end
    elseif community_algorithm == "covariance"
        # Covariance
        if covariance_cutoff_fraction == nothing
            throw(ArgumentError("covariance_cutoff_fraction must be set if 'community_algorithm' == 'covariance'"))
        end
        if covariance_data == nothing
            throw(ArgumentError("covariance_data must be set if 'community_algorithm' == 'covariance'"))
        end

        covariances = cov(covariance_data, dims=:time_step)

        communities = 2*ones(Int64,nv(graph))
        communities[map(x -> x[2], findall(sum(covariances,dims=1) .>= covariance_cutoff_fraction*maximum(covariances)))] .= 1
    elseif community_algorithm == "walktrap"
        if walktrap_steps == nothing
            throw(ArgumentError("walktrap_steps must be set if 'community_algorithm' == 'walktrap'"))
        end

        communities_igraph = IGVectorInt(zeros(nv(graph)))
				LibIGraph.community_walktrap(
						IGraph(SimpleDiGraph(graph)),
						IGVectorFloat(weight.(edges(graph))),
						walktrap_steps,
						IGNull(),
						IGNull(),
						communities_igraph,
				)
				communities = Integer.(communities_igraph) .+ 1
				communities = Integer.(communities_igraph)
    elseif community_algorithm == "leiden"
        if community_resolution == nothing
            throw(ArgumentError("community_resolution must be set if 'community_algorithm' == 'leiden'"))
        end
        if community_beta == nothing
            throw(ArgumentError("community_beta must be set if 'community_algorithm' == 'leiden'"))
        end
        if community_n_iter == nothing
            throw(ArgumentError("community_n_iter must be set if 'community_algorithm' == 'leiden'"))
        end

        communities_igraph = IGVectorInt(zeros(nv(graph)))
        LibIGraph.community_leiden_simple(
				    IGraph(SimpleDiGraph(graph)),
						IGVectorFloat(weight.(edges(graph))),
						LibIGraph.IGRAPH_LEIDEN_OBJECTIVE_CPM,
						community_resolution, community_beta, false, community_n_iter,
						communities_igraph)
        communities = Integer.(communities_igraph) .+ 1
    else
        throw(ArgumentError("community_algorithm must be a string in set [\"label-propagation\", "
                            * "\"strongly-connected\", \"covariance\", \"walktrap\", \"leiden\"]"))
    end

    return communities
end

function order_parameters_by_community(
															initial_actions::AbstractVector{<:Integer},
															deltas::AbstractMatrix{<:Integer},
															communities::AbstractVector{<:Integer},
															nb_phases::Integer,
															)

   num_communities = length(unique(communities))
   num_times = size(deltas,2)
   order_parameters = Matrix{Float64}(undef,num_times+1,num_communities)

   current_actions = initial_actions

   player_indices_by_community = [communities .== community_idx
     for community_idx in 1:num_communities]

   # Populate initial phase parameters for each community
   for community_idx = 1:num_communities
	   order_parameters[0+1, community_idx] = extract_order_parameters(
	     extract_counts(current_actions[player_indices_by_community[community_idx]],
	       nb_phases),nb_phases)
   end

   # Only update phase parameter for affected community
   for time_idx = 1:num_times
	   order_parameters[time_idx+1,:] = order_parameters[time_idx,:]
	   current_actions += deltas[:,time_idx]

	   changed_idxs = findall(deltas[:,time_idx] .!= 0)
	   changed_community_idxs = unique(communities[changed_idxs])

	   for changed_community_idx in changed_community_idxs
		   order_parameters[time_idx+1, changed_community_idx] = extract_order_parameters(
		     extract_counts(current_actions[player_indices_by_community[changed_community_idx]],
		       nb_phases),nb_phases)
	   end
   end

	 return order_parameters
end

function get_chimera_indices(order_parameters::AbstractMatrix)
	 # Important: order_parameters should be m-by-n array
	 # with m time steps and n communities
   metastability = mean(var(order_parameters;dims=1))
   chimera_index = mean(var(order_parameters;dims=2))

   results = Dict("metastability_index"=>metastability, "chimera_index"=>chimera_index)
   return results
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

function generate_nodes_edges(graph::AbstractGraph;
                data::Union{Dict,Nothing}=nothing,
                time_step::Union{Integer,Nothing}=nothing)
    # Get vertex coordinates
    vertex_coordinates = stress(graph)

    # Split points into coordinates
    xs = (p -> p[1]).(vertex_coordinates)
    ys = (p -> p[2]).(vertex_coordinates)

    # Get edges
    graph_edges = collect(edges(graph))

		if isnothing(data)
            # Put into DataFrame
						nodes_df = DataFrame(index=1:nv(graph), x=xs, y=ys)
    else
	    # Decode data to get vertex labels
	    vertex_labels = decode_delta_encoded(data["initial_actions"],
	      data["deltas"],time_step)

            # Put into DataFrame
						nodes_df = DataFrame(index=1:nv(graph), x=xs, y=ys, strategyIndex=vertex_labels)
    end

    return nodes_df, graph_edges
end

function communicative_fraction_theory(graph::AbstractGraph;
                              type::String,
															selection_strength::Real,
                              symmetry_breaking::Real,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              nb_players::Integer=20,
			      )

		nb_effective = ( mean(indegree(graph))+1) # Add one since n=degree+1 for well-mixed case

		beta0(B0) = beta_to_B*B0
		if type == "theory"
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

				func = B0 -> analytic_frac_communicative(B0, beta0(B0);
						selection_strength, cost, nb_players=nb_effective, symmetry_breaking, nb_phases)
		elseif type == "approx"
				func = B0 -> 1 / (1 + exp(selection_strength * (nb_effective - 1) *
													((nb_effective - 1) * cost
							 + nb_effective * beta0(B0) * (1-2*symmetry_breaking)/2
													- (nb_effective - 2) / 2 * B0)))
		else
				throw(ArgumentError("type must be a string in set ["\"theory\", \"approx\"]"))
		end

		# Determine B range
		B_crit = 2 * cost * (nb_players - 1) / (nb_players - 2)
		nb_Bs = 11
		step_size_Bs = 0.04
		B_max = B_crit + (nb_Bs - 1)/2 * step_size_Bs
		B_min = B_crit - (nb_Bs - 1)/2 * step_size_Bs

		Bs, communicative_fraction = PlotUtils.adapted_grid(func, (B_min, B_max))
		return @strdict(Bs, communicative_fraction)
end

@memoize function game_types_per_strategy_pair(mutual_benefit_synchronous::Real,
																							 unilateral_benefit_synchronous::Real,
																							 cost::Real,
																							 symmetry_breaking::Real,
																							 nb_phases::Integer;
																							 only_mixed_games::Bool=false,
								 )
		# Choose game type by assuming either player can switch (strategy,phase) to the other player, giving a 2x2 payoff matrix

		function game_type(payoff_matrix::AbstractMatrix)
				# Schema: doi:10.3390/g6040495

				# Make ordinal using modified compete rank
				ordinal_payoffs = 5 .- competerank(-payoff_matrix)

				# Orient in left-up:
				# put highest row-player payoff in left column
				# and highest col-player payoff in top row
				#
				# (Since we are dealing with symmetric games,
				# we only need to ensure the row-player condition,
				# as the col-player condition is automatically satisfied)
				#
				# (side note: we use left-up convention,
				# so the col-player's matrix is the transpose of the row-player's matrix;
				# in right-up convention, the col-player's matrix is the
				# "anti-transpose" of the row-player's matrix,
				# ie a transpose across the anti-diagonal)
				function canonical_payoff!(matrix)
					function swap_strategies!(payoff_matrix::AbstractMatrix)
							return payoff_matrix[2, 2], payoff_matrix[2, 1], payoff_matrix[1, 2], payoff_matrix[1, 1] =
								payoff_matrix[1, 1], payoff_matrix[1, 2], payoff_matrix[2, 1], payoff_matrix[2, 2]
					end

					# Put in left-up orientation
					# by ensuring (at least one of)
					# the highest entry (4) # is in the left column
					if ! (4 in matrix[:,1])
						swap_strategies!(matrix)
					end
				end
				canonical_payoff!(ordinal_payoffs)

				# Determine ties
				binomial_nomenclature = game_taxonomy[ordinal_payoffs]

				# Drop tie information
				local game_type = binomial_nomenclature[2]

				return game_type
		end

		# Define payoff response submatrices
		payoff = payoff_matrix(nb_phases, mutual_benefit_synchronous,
													 unilateral_benefit_synchronous, cost; symmetry_breaking)

		# Define game type per strategy pair
		game_types = similar(payoff, Union{GameType,Missing})
		for idx in eachindex(IndexCartesian(), game_types)
				row_player_idx = idx[1]
				col_player_idx = idx[2]
				# Note: the indexing ensures we are in standard symmetric
				# bimatrix representation (i.e. the main diagonal has the same
				# payoffs for both players)
				game_types[idx] = game_type([
					 [payoff[row_player_idx,row_player_idx], payoff[col_player_idx,row_player_idx]] [payoff[row_player_idx,col_player_idx], payoff[col_player_idx,col_player_idx]]
					])
		end

		if only_mixed_games
			game_types[1:nb_phases,1:nb_phases] .= missing
			game_types[nb_phases+1:end,nb_phases+1:end] .= missing
		end
		return game_types
end

function calc_timeseries_statistics(initial_actions::AbstractVector{<:Integer}, deltas::AbstractMatrix{<:Integer}, nb_phases::Integer,
  symmetry_breaking::Real, B_to_c::Real, beta_to_B::Real, cost::Real, interaction_adj_matrix::AbstractMatrix{<:Integer};
  only_mixed_games::Bool=false)
    # Extract results
		nb_players = size(interaction_adj_matrix,2)
    time_steps = size(deltas,2)
    fraction_communicative = Vector{Float64}(undef,time_steps+1)
    order_parameters = Vector{Float64}(undef,time_steps+1)
    most_common_game_types = Vector{GameType}(undef,time_steps+1)

    current_actions = initial_actions
    interaction_adj_matrix_transpose = copy(transpose(interaction_adj_matrix))

    # Initial statistics
    time_idx = 0
    game_types = game_types_per_strategy_pair(B_to_c*cost,
                                              beta_to_B*B_to_c*cost, cost,
                                              symmetry_breaking, nb_phases; only_mixed_games)
    games, game_counts = count_games(current_actions;
      game_types, interaction_adj_matrix, interaction_adj_matrix_transpose)
    # Check if all cooperative, all non-cooperative, or disconnected C/N
    # Do this after the loop so games and game_counts are updated
    all_same_strategy = check_all_same_strategy(current_actions, nb_phases)
    disconnected_synchronized = check_disconnected_synchronized(game_counts)
    if !isnothing(all_same_strategy)
      most_common_game_types[time_idx+1] = all_same_strategy
		elseif !isnothing(disconnected_synchronized)
      most_common_game_types[time_idx+1] = disconnected_synchronized
		else
    	most_common_game_types[time_idx+1] = extract_most_common_game_types(game_counts)
    end
    action_counts_per_timestep = extract_counts(current_actions, nb_phases)
    nb_communicative = extract_num_communicative(action_counts_per_timestep)
    fraction_communicative[time_idx+1] = nb_communicative / nb_players
    order_parameters[time_idx+1] = extract_order_parameters(action_counts_per_timestep, nb_phases)

    for time_idx in 1:time_steps
      current_actions += deltas[:,time_idx]
	    games, game_counts = count_games(current_actions, games, game_counts, deltas[:,time_idx];
				game_types, interaction_adj_matrix, interaction_adj_matrix_transpose)
      # Check if all cooperative, all non-cooperative, or disconnected C/N
			# Do this after the loop so games and game_counts are updated
			all_same_strategy = check_all_same_strategy(current_actions, nb_phases)
      disconnected_synchronized = check_disconnected_synchronized(game_counts)
			if !isnothing(all_same_strategy)
				most_common_game_types[time_idx+1] = all_same_strategy
			elseif !isnothing(disconnected_synchronized)
				most_common_game_types[time_idx+1] = disconnected_synchronized
			else
	    	most_common_game_types[time_idx+1] = extract_most_common_game_types(game_counts)
			end
	    action_counts_per_timestep = extract_counts(current_actions, nb_phases)
	    nb_communicative = extract_num_communicative(action_counts_per_timestep)
	    fraction_communicative[time_idx+1] = nb_communicative / nb_players
	    order_parameters[time_idx+1] = extract_order_parameters(action_counts_per_timestep, nb_phases)
    end

    # Package results
    return @strdict(fraction_communicative, order_parameters, most_common_game_types)
end

function downsample_data(data::AbstractDataFrame, num_samples::Integer=1000)
    # Only save a subset of points to prevent large file sizes
    downsample_ratio = Int(floor((size(data,1) + 1) / num_samples))
    if downsample_ratio == 0
        throw(ErrorException("Stride between samples is less than 1; either decrease `num_samples`, increase `time_steps`, or increase `early_cutoff_fraction`"))
    end

    # Downsample
    # Note: the populations include the initial statistics, so we need one more than time-steps
    data = data[1:downsample_ratio:end,:]

		return data
end

function calc_number_unidirection_bidirectional_edges(graph::AbstractGraph)
    # Calculate self-loops
    loops = num_self_loops(graph)

    # Convert to unweighted, undirected graph
    undirected_graph = SimpleGraph(graph)

		# Each bidirectional edge *pair* adds a single edge to the directed
		# graph
		bidirectional_edge_pairs = ne(graph) - ne(undirected_graph)

    # Degree counts both in- and out-neighbors,
    # so number of unidirectional edges is half the degree
		unidirectional_edges = ne(graph) - loops - 2*bidirectional_edge_pairs

    # Calculate bidirectional edges
    results = Dict("self-loops" => loops,
		  "unidirectional_edges" => unidirectional_edges,
		  "bidirectional_edge_pairs" => bidirectional_edge_pairs)

		return results
end
