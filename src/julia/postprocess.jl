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
using CondaPkg
using PythonCall
using JLD2
using YAXArrays
using NetCDF

include("moran.jl")
include("game_taxonomy.jl")
include("utils.jl")

function load_all_cumulative(time_steps::Integer=200_000_000)
    # Load dataframe
    df_raw = collect_results(datadir("raw","cumulative"); rinclude = [Regex("time_steps=$time_steps[._]")])

    # Add path
    df = transform(df_raw, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)

    return df
end

function load_all_timeseries(time_steps::Integer=80_000)
    # Load dataframe
    df_raw = collect_results(datadir("raw","timeseries"); rinclude = [Regex("time_steps=$time_steps[._]")])

    # Add path
    df = transform(df_raw, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)

    return df
end

# Source: doi:10.3390/g6040495
# In left-up convention (modified from right-up convention by swapping columns)
# = means exactly the same payoff matrix, \approx means equivalent up to swap_strategies!
const game_taxonomy = Dict(
			 [4 2;3 1] => (missing, concord),
			 [4 3;2 1] => (missing, harmony),
			 [4 3;1 2] => (missing, peace),
			 [4 2;1 3] => (missing, coordination),
			 [4 1;2 3] => (missing, assurance),
			 [4 1;3 2] => (missing, staghunt),
			 [3 1;4 2] => (missing, dilemma),
			 [2 1;4 3] => (missing, deadlock),
			 [1 2;4 3] => (missing, compromise),
			 [1 3;4 2] => (missing, hero),
			 [2 3;4 1] => (missing, battle),
			 [3 2;4 1] => (missing, chicken),
			 [2 3;4 2] => (low, battle),
			 [2 2;4 3] => (low, deadlock),
			 [3 2;4 2] => (low, dilemma),
			 [4 2;2 3] => (low, coordination),
			 [4 3;2 2] => (low, harmony),
			 [4 2;3 2] => (low, concord),
			 [3 3;4 1] => (mid, battle),
			 [1 3;4 3] => (mid, compromise),
			 [3 1;4 3] => (mid, deadlock),
			 [4 1;3 3] => (mid, staghunt),
			 [4 3;1 3] => (mid, peace),
			 [4 3;3 1] => (mid, harmony),
			 [1 4;4 2] => (high, hero),
			 [2 4;4 1] => (high, hero), # high battle \approx high hero
			 [4 2;4 1] => (high, concord), # = high chicken
			 [4 1;4 2] => (high, staghunt), # = high dilemma
			 [4 2;1 4] => (high, coordination),
			 [4 1;2 4] => (high, coordination), # high assurance \approx high coord
			 [4 4;1 2] => (high, peace),
			 [2 1;4 4] => (high, peace), # high deadlock \approx high peace
			 [4 4;2 1] => (high, harmony),
			 [1 2;4 4] => (high, harmony), # high compromise \approx high harmony
			 [4 2;4 2] => (double, staghunt), # = double dilemma; note: bruns2015 claims double dilemma = double dilemma, which is tautologically true, but (correctly) states double dilemma = double staghunt elsewhere
			 [4 2;2 4] => (double, coordination),
			 [2 4;4 2] => (double, coordination), # double hero \approx double coordination; note: this is the only equivalent pair that is not related by swap_strategies!;
			 # Instead, they are related by swapping only the columns or rows, but not both. Usually, doing this produces a non-symmetric game, but it this case, it is still symmetric
			 [4 4;2 2] => (double, harmony),
			 [2 2;4 4] => (double, harmony), # double compromise \approx double harmony
			 [4 4;1 4] => (triple, deadlock),
			 [4 1;4 4] => (triple, deadlock), # note: not explicitly given in bruns2015, but equivalent under swap_strategies!
			 [4 4;4 1] => (triple, harmony),
			 [1 4;4 4] => (triple, harmony), # note: not explicitly given in bruns2015, but equivalent under swap_strategies!
			 [3 3;4 3] => (basic, dilemma),
			 [4 3;3 3] => (basic, harmony),
			 [4 4;4 4] => (zero, neutral),
			 )

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
    game_type = binomial_nomenclature[2]
end

@memoize function game_types_per_strategy_pair(mutual_benefit_synchronous::Real,
                                               unilateral_benefit_synchronous::Real,
                                               cost::Real,
                                               symmetry_breaking::Real,
                                               nb_phases::Integer)
    # Choose game type by assuming either player can switch (strategy,phase) to the other player, giving a 2x2 payoff matrix

    # Define payoff response submatrices
    payoff = payoff_matrix(nb_phases, mutual_benefit_synchronous,
                           unilateral_benefit_synchronous, cost; symmetry_breaking)

    # Define game type per strategy pair
    game_types = similar(payoff, GameType)
    for idx in eachindex(IndexCartesian(), game_types)
	row_player_idx = idx[1]
	col_player_idx = idx[2]
        game_types[idx] = game_type([
			     [payoff[row_player_idx,col_player_idx], payoff[col_player_idx,col_player_idx]] [payoff[row_player_idx,row_player_idx], payoff[col_player_idx,row_player_idx]]
			    ])
    end
    return game_types
end

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
        row_phase = strategies_per_player[row]
        col_phase = strategies_per_player[col]
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
    count_actions!(counts, strategies_per_player)
    return counts
end

function generate_communities(graph::AbstractSimpleWeightedGraph, community_algorithm::String;
        covariance_cutoff::Union{Real,Nothing}=nothing,
        covariance_data::Union{DimArray,Nothing}=nothing)
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
    elseif community_algorithm == "infomap"
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
    elseif community_algorithm == "covariance"
        # Covariance
        if covariance_cutoff == nothing
            throw(ArgumentError("covariance_cutoff must be set if 'community_algorithm' == 'covariance'"))
        end
        if covariance_data == nothing
            throw(ArgumentError("covariance_data must be set if 'community_algorithm' == 'covariance'"))
        end

        covariances = cov(covariance_data, dims=:time_step)

        communities = 2*ones(Int64,nv(graph))
        communities[map(x -> x[2], findall(sum(covariances,dims=1) .>= covariance_cutoff))] .= 1
    else
        throw(ArgumentError("community_algorithm must be a string in set [\"label-propagation\", "
                            * "\"strongly-connected\", \"infomap\", "
                            * "\"covariance\"]"))
    end

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

function calc_coalescence_times(adj_matrix_source::String="well-mixed")
    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
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

function export_graph_nodes_edges(; time_step::Union{Real,Nothing}=nothing,
                              B_to_c::Union{Real,Nothing}=nothing, selection_strength::Union{Real,Nothing}=nothing,
                              symmetry_breaking::Union{Real,Nothing}=nothing,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              nb_players::Integer=20,
			      )
    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)

    if time_step == nothing
            # Generate configuration
            config = @strdict(adj_matrix_source)

            # Calculate nodes and edges
            nodes_df, graph_edges = generate_nodes_edges(graph)
    else
            # Generate configuration
            config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                              selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
            if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
                config["nb_players"] = nb_players
            end

            # Get data
            data = wload(datadir("raw", "timeseries", savename(config,"jld2")))

            # Calculate nodes and edges
            nodes_df, graph_edges = generate_nodes_edges(graph; data, time_step)

            # Add time_step to config dictionary
            config["time_step"] = time_step
    end

    # Add index number
    insertcols!(nodes_df, 1, :index => 1:nrow(nodes_df))

    # Write out vertices and edges
    CSV.write(datadir("processed", "graph_structure", savename("vertices",config,"csv")), nodes_df)
    CSV.write(datadir("processed", "graph_structure", "edges_adj_matrix_source=$(adj_matrix_source).csv"), graph_edges)

    return graph_edges
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

    if data == nothing
            # Put into DataFrame
            df = DataFrame(x=xs, y=ys)
    else
            # Get vertex labels
            vertex_labels = data["all_populations"][:,time_step]

            # Put into DataFrame
            df = DataFrame(x=xs, y=ys, strategyIndex=vertex_labels)
    end

    return df, graph_edges
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

function extract_cumulative(; type::String, selection_strength::Real,
                              symmetry_breaking::Real,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              nb_players::Integer=20,
			      )

    # Generate configuration
    config = @strdict(adj_matrix_source, time_steps, beta_to_B,
                      selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end

    # Get data
    data = wload(datadir("raw","cumulative",savename(config,"jld2")))

    if type == "simulation"
        df = DataFrame(B0=data["Bs"],
                       communicative_fraction=data["fraction_communicative"])
    elseif type == "theory" || type == "approx"
        # Generate graph
        interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
        graph = SimpleDiGraph(interaction_adj_matrix)
        nb_effective = ( mean(indegree(graph))+1) # Add one since n=degree+1 for well-mixed case

        beta0(B0) = beta_to_B*B0
        if type == "theory"
            func = B0 -> analytic_frac_communicative(B0, beta0(B0);
                selection_strength=config["selection_strength"], cost=config["cost"],
                nb_players=nb_effective, symmetry_breaking=config["symmetry_breaking"],
                nb_phases=config["nb_phases"])
        else
            func = B0 -> 1 / (1 + exp(config["selection_strength"] * (nb_effective - 1) *
                              ((nb_effective - 1) * config["cost"]
                   + nb_effective * beta0(B0) * (1-2*config["symmetry_breaking"])/2
                              - (nb_effective - 2) / 2 * B0)))
        end

        Bs, frac = PlotUtils.adapted_grid(func, (minimum(data["Bs"]), maximum(data["Bs"])))
        df = DataFrame(B0 = Bs, communicative_fraction=frac)
    else
        throw(ArgumentError("type must be a string in set [\"simulation\", "
                            * "\"theory\", \"approx\"]"))
    end

    # Add type to config dictionary
    config["type"] = type

    # Write out data
    CSV.write(datadir("processed","cumulative",savename(config,"csv")), df)
end

function calc_timeseries_statistics(all_populations::AbstractMatrix{<:Integer}, nb_phases::Integer, nb_players::Integer,
  symmetry_breaking::Real, B_to_c::Real, beta_to_B::Real, cost::Real, interaction_adj_matrix::AbstractMatrix{<:Integer})
    # Extract results
    most_common_game_types = dropdims(mapslices(x -> extract_most_common_game_types(x,
                                                                                    B_to_c *
                                                                                    cost,
                                                                                    beta_to_B *
                                                                                    B_to_c *
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
    return @strdict(fraction_communicative, order_parameters, most_common_game_types, strategy_parity)
end

function extract_timeseries_statistics(; B_to_c::Real, selection_strength::Real,
                              symmetry_breaking::Real,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              num_samples::Integer=1000,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              nb_players::Integer=20,
			      )

    # Create config dict for saving filename
    config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                      selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
      config["nb_players"] = nb_players
    end

    # Create a dictionary with only the parameters that affect the RNG
    rng_config = Dict(key => config[key] for key in ["adj_matrix_source", "time_steps"])

    # Get data and load into dataframe
    data_dict = wload(datadir("raw", "timeseries", savename(config, "jld2")))
    # Add configuration
    data_dict = merge(data_dict, config)
    # Wrap all elements in a list to allow for matrices in individual
    # DataFrame elements
    data_dict = Dict(k => [v] for (k,v) in data_dict)
    # Convert to dataframe
    df_all_asymm = DataFrame(data_dict)

    # Generate statistics
    statistics = transform(df_all_asymm, [:all_populations, :nb_phases, :nb_players, :symmetry_breaking,
        :B_to_c, :beta_to_B, :cost, :interaction_adj_matrix] => ByRow(calc_timeseries_statistics) => AsTable)

    # Replace game type if all communicative/noncommunicative
    transform!(statistics, [:strategy_parity, :most_common_game_types] => ByRow((strategy_parity, game_type) -> strategy_parity != mixed ? strategy_parity :
      game_type) => :game_type_or_parity)

    # Write out mutation times
    CSV.write(datadir("processed","mutation_timesteps", savename(rng_config,"csv")),
      rename(select(statistics, :steps_following_mutation), Dict(:steps_following_mutation => :mutation_timesteps)))

    # Select subset of columns and rename
    rename!(select!(statistics,
      [:fraction_communicative, :order_parameters , :game_type_or_parity]),
      Dict(:fraction_communicative => :communicative_fraction, :order_parameters => :order_parameter, :game_type_or_parity => :game_type))

    # Add time
    insertcols!(statistics, 1, :time => 1:nrow(statistics))

    # Only plot subset of points to prevent large file sizes
    downsample_ratio = Int(floor((time_steps + 1) / num_samples))

    # Downsample
    # Note: the populations include the initial statistics, so we need one more than time-steps
    statistics = statistics[1:downsample_ratio:end,:]

    # Write out statistics
    CSV.write(datadir("processed","timeseries_statistics",savename(config,"csv")), statistics)
end

function extract_chimera_indices(; community_algorithm::String,
                              B_to_c::Real, selection_strength::Real,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              covariance_cutoff::Real,
                              nb_players::Integer=20,
			      )

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
    graph = SimpleWeightedDiGraph(interaction_adj_matrix)

    # Generate configuration
    config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                      selection_strength, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end

    df_all_asymm = DataFrame()
    for symmetry_breaking in [0,0.25,0.5,0.75,1.0]
      config_with_asymmetry = merge(config,Dict("symmetry_breaking" => symmetry_breaking))
      # Get data and load into dataframe
      data_dict = wload(datadir("raw", "timeseries",
          savename(config_with_asymmetry, "jld2")))
      # Add configuration
      data_dict = merge(data_dict, config_with_asymmetry)
      # Wrap all elements in a list to allow for matrices in individual
      # DataFrame elements
      data_dict = Dict(k => [v] for (k,v) in data_dict)
      # Append to dataframe
      append!(df_all_asymm, data_dict)
    end

    # Get communities
    if community_algorithm == "covariance"
        # Choose asymmetry=0.75 for reference covariance
        symmetry_breaking_ref = 0.75
        df_ref = @rsubset(df_all_asymm, :symmetry_breaking == symmetry_breaking_ref)

        if nrow(df_ref) < 1
            throw(ErrorException("Did not find any timeseries data with `symmetry_breaking`=$(symmetry_breaking_ref)"))
        elseif nrow(df_ref) > 1
            throw(ErrorException("Found multiple timeseries data with `symmetry_breaking`=$(symmetry_breaking_ref)"))
        end

        # Convert to DimArray
        data_ref = DimArray(df_ref.all_populations[1], (:player_index, :time_step))

        # Only use subset 1% of data to calculate communities
        data_ref = data_ref[time_step = At(1:ceil(time_steps*0.01))]

        communities = generate_communities(graph, community_algorithm;
                                           covariance_cutoff=covariance_cutoff, covariance_data = data_ref)

        # Add covariance_cutoff to config dictionary
        config["covariance_cutoff"] = covariance_cutoff
    else
        communities = generate_communities(graph, community_algorithm)
    end

    # Convert to DimArray
    transform!(df_all_asymm,
               :all_populations => ByRow(array  -> DimArray(array, (:player_index, :time_step))) => :all_populations)

    # Get chimera indices
    transform!(df_all_asymm, :all_populations => ByRow(pop -> get_chimera_indices(pop, communities, nb_phases)) => AsTable)

    # Only keep columns we're interested in
    df = select(df_all_asymm, :symmetry_breaking => :asymmetry,
                :chimera_index, :metastability_index)

    # Add community_algorithm to config dictionary
    config["community_algorithm"] = community_algorithm

    # Write out data
    CSV.write(datadir("processed","chimeraindex",savename(config,"csv")), df)
end

function extract_game_types(; B_to_c::Real, selection_strength::Real,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              nb_players::Integer=20,
			      )

    # Generate configuration
    config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                      selection_strength, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end

    df_all_asymm = DataFrame()
    for symmetry_breaking in [0,0.25,0.5,0.75,1.0]
      config_with_asymmetry = merge(config,Dict("symmetry_breaking" => symmetry_breaking))
      # Get data and load into dataframe
      data_dict = wload(datadir("raw", "timeseries",
          savename(config_with_asymmetry, "jld2")))
      # Add configuration
      data_dict = merge(data_dict, config_with_asymmetry)
      # Wrap all elements in a list to allow for matrices in individual
      # DataFrame elements
      data_dict = Dict(k => [v] for (k,v) in data_dict)
      # Append to dataframe
      append!(df_all_asymm, data_dict)
    end

    # Generate statistics
    transform!(df_all_asymm, [:all_populations, :nb_phases, :nb_players, :symmetry_breaking,
        :B_to_c, :beta_to_B, :cost, :interaction_adj_matrix] => ByRow(calc_timeseries_statistics) => AsTable)

    # Combine parity and game type
    transform!(df_all_asymm,
	       [:strategy_parity, :most_common_game_types] => ByRow((parity_col,game_col) ->
								    [parity != mixed ? parity : game for
								     (parity, game) in zip(parity_col, game_col)])
	       => :parity_or_game_type)
    game_types = proportionmap.(df_all_asymm.parity_or_game_type)

    # Convert Dict keys from GameType to Symbol
    game_type_symbols = (dict -> Dict(Symbol(k) => v for (k,v) in pairs(dict))).(game_types)

    # Combine asymmetries into a single data frame
    df_missing = vcat(DataFrame.(game_type_symbols)...; cols=:union)

    # Sort columns alphabetically
    select!(df_missing, sort(names(df_missing)))

    # Add asymmetry
    insertcols!(df_missing, 1, :asymmetry => df_all_asymm.symmetry_breaking)

    # Replace missing data (i.e. game types that do not appear for a particular asymmetry) with zero
    df = coalesce.(df_missing, 0.0)

    # Write out data
    CSV.write(datadir("processed","gametype",savename(config,"csv")), df)
end

function calc_number_unidirection_bidirectional_edges(; adj_matrix_source::String)
    # Note: removes edge weights and self-loops
    interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
    # Form graph
    graph = SimpleWeightedDiGraph(interaction_adj_matrix)
    # Count total number of connections
    total_edges = sum(indegree(graph))
    # Calculate self-loops
    loops = num_self_loops(graph)
    # Remove self-loops
    non_loop_edges = total_edges - loops

    # Remove edge weighting
    unweighted_adj_matrix = collect(round.(interaction_adj_matrix) .!= 0)
    # Xor with transpose to only keep unidirectional edges
    unidirectional_edge_adj_matrix = xor.(unweighted_adj_matrix,
					  transpose(unweighted_adj_matrix))
    # Convert to graph
    unidirectional_edge_graph = SimpleGraph(unidirectional_edge_adj_matrix)
    # Degree counts both in- and out-neighbors,
    # so number of unidirectional edges is half the degree
    unidirectional_edges = sum(degree(unidirectional_edge_graph))/2

    # Bidirectional edges is total number of edges minus unidirectional edges
    # Bidirectional edge *pairs* is half this number
    bidirectional_edge_pairs = (non_loop_edges - unidirectional_edges)/2

    # Calculate bidirectional edges
    results = Dict("self-loops" => loops,
		  "unidirectional_edges" => unidirectional_edges,
		  "bidirectional_edge_pairs" => bidirectional_edge_pairs)

    config = @strdict(adj_matrix_source)
    CSV.write(datadir("processed", "graph_loop_edge_number", savename(config,"csv")), results)
end

function create_netcdf(;data_type::String, adj_matrix_source::String, time_steps::Integer, decimation_factor::Union{Nothing,Integer}=nothing)
	function get_properties(df,adj_matrix_source; include_time_steps::Bool=true)
		property_list = ["nb_phases", "adj_matrix_source", "cost", "mutation_rate"]
		if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
			append!(property_list, ["nb_players"])
		end
		if include_time_steps
			append!(property_list, ["time_steps"])
		end
		properties_df = unique(select(df, property_list))
		if nrow(properties_df) != 1
			throw(ErrorException("Datasets with adj_matrix_source==$adj_matrix_source do not have identical properties=$property_list"))
		end
		properties_dict = Dict(names(properties_df[1,:]) .=> values(properties_df[1,:]))
		return properties_dict
	end

  if data_type == "cumulative"
    # Load data
    df_cumulative = @rsubset(load_all_cumulative(time_steps), :matrix_source == adj_matrix_source)

    # Ensure column name consistency
    rename!(df_cumulative, :matrix_source => :adj_matrix_source)

    properties_dict_cumulative = get_properties(df_cumulative,adj_matrix_source; include_time_steps=true)

    # Combine into a YAXArray
    axes = (
      Dim{:maximum_joint_benefit}(df_cumulative.Bs[1]),
      Dim{:symmetry_breaking}(unique(df_cumulative.symmetry_breaking)),
      Dim{:selection_strength}(unique(df_cumulative.selection_strength)),
      )
    cumulative = YAXArray(axes,
              stack(only(@rsubset(df_cumulative, :symmetry_breaking == alpha,  :selection_strength == delta)).fraction_communicative
              for alpha in axes[2], delta in axes[3]),
              properties_dict_cumulative,
              )
    savecube(cumulative, datadir("processed","netcdf","cumulative_matrixSource=$(adj_matrix_source)_timesteps=$time_steps.nc"), driver=:netcdf, overwrite=true)
  elseif data_type == "timeseries-statistics"
    # Load data
    df_timeseries = @rsubset(load_all_timeseries(time_steps), :adj_matrix_source == adj_matrix_source)
    transform!(df_timeseries, [:all_populations, :nb_phases, :nb_players, :symmetry_breaking,
          :to_c, :beta_to_B, :cost, :interaction_adj_matrix] => ByRow(calc_timeseries_statistics) => AsTable)

    properties_dict_timeseries = get_properties(df_timeseries,adj_matrix_source; include_time_steps=false)
    transform!(df_timeseries, :to_c => ByRow(x -> x*properties_dict_timeseries["cost"]) => :maximum_joint_benefit)

    nb_players = size(df_timeseries.all_populations[1])[1]

    symmetry_breaking_timeseries_vals = unique(df_timeseries.symmetry_breaking)
    selection_strength_timeseries_vals = unique(df_timeseries.selection_strength)
    maximum_joint_benefit_timeseries_vals = unique(df_timeseries.maximum_joint_benefit)
    axes_timeseries = (
      Dim{:time_step}(0:time_steps),
      Dim{:symmetry_breaking}(symmetry_breaking_timeseries_vals, span=Regular(0.25)),
      Dim{:maximum_joint_benefit}(round.(maximum_joint_benefit_timeseries_vals; digits=5)),
      Dim{:selection_strength}(selection_strength_timeseries_vals, span=Regular(4.8)),
      )
    axes_timeseries_with_players = (Dim{:player_index}(1:nb_players), axes_timeseries...)

    # Combine into a YAXArray
    transform!(df_timeseries, :most_common_game_types => ByRow(x -> Integer.(x)) => :most_common_game_types)
    timeseries = YAXArray(axes_timeseries_with_players,
              stack(only(
                   begin
                     x = @rsubset(df_timeseries, :symmetry_breaking == alpha,  :selection_strength == delta, :maximum_joint_benefit == B_0)
                     !isempty(x) ? x : DataFrame(all_populations = [fill(missing, size(axes_timeseries[[1,2]]))])
                   end
             ).all_populations
              for alpha in symmetry_breaking_timeseries_vals, B_0 in maximum_joint_benefit_timeseries_vals,
              delta in selection_strength_timeseries_vals), # Use *_vals instead elements of axes_timeseries because axes rounds maximum_joint_benefit
              properties_dict_timeseries,
              )
    most_common_game_types = YAXArray(
                axes_timeseries,
                stack(only(
                     begin
                       x = @rsubset(df_timeseries, :symmetry_breaking == alpha,  :selection_strength == delta, :maximum_joint_benefit == B_0)
                       !isempty(x) ? x : DataFrame(most_common_game_types = [fill(missing, size(axes_timeseries[1]))])
                     end
                     ).most_common_game_types
                for alpha in symmetry_breaking_timeseries_vals, B_0 in maximum_joint_benefit_timeseries_vals,
                delta in selection_strength_timeseries_vals), # Use *_vals instead elements of axes_timeseries because axes rounds maximum_joint_benefit
                merge(properties_dict_timeseries, Dict("enum_lookup:" .* string.(Integer.(instances(GameType)))
                     .=> String.(Symbol.(instances(GameType))))),

               )
    order_parameters = YAXArray(
                axes_timeseries,
                stack(only(
                     begin
                       x = @rsubset(df_timeseries, :symmetry_breaking == alpha,  :selection_strength == delta, :maximum_joint_benefit == B_0)
                       !isempty(x) ? x : DataFrame(order_parameters = [fill(missing, size(axes_timeseries[1]))])
                     end
                     ).order_parameters
                for alpha in symmetry_breaking_timeseries_vals, B_0 in maximum_joint_benefit_timeseries_vals,
                delta in selection_strength_timeseries_vals), # Use *_vals instead elements of axes_timeseries because axes rounds maximum_joint_benefit
                properties_dict_timeseries,
               )
    transform!(df_timeseries, :strategy_parity => ByRow(x -> Integer.(x)) => :strategy_parity)
    strategy_parity = YAXArray(
                axes_timeseries,
                stack(only(
                     begin
                       x = @rsubset(df_timeseries, :symmetry_breaking == alpha,  :selection_strength == delta, :maximum_joint_benefit == B_0)
                       !isempty(x) ? x : DataFrame(strategy_parity = [fill(missing, size(axes_timeseries[1]))])
                     end
                     ).strategy_parity
                for alpha in symmetry_breaking_timeseries_vals, B_0 in maximum_joint_benefit_timeseries_vals,
                delta in selection_strength_timeseries_vals), # Use *_vals instead elements of axes_timeseries because axes rounds maximum_joint_benefit
                merge(properties_dict_timeseries, Dict("enum_lookup:" .* string.(Integer.(instances(StrategyParity)))
                     .=> String.(Symbol.(instances(StrategyParity))))),
               )
    fraction_communicative = YAXArray(
                axes_timeseries,
                stack(only(
                     begin
                       x = @rsubset(df_timeseries, :symmetry_breaking == alpha,  :selection_strength == delta, :maximum_joint_benefit == B_0)
                       !isempty(x) ? x : DataFrame(fraction_communicative = [fill(missing, size(axes_timeseries[1]))])
                     end
                     ).fraction_communicative
                for alpha in symmetry_breaking_timeseries_vals, B_0 in maximum_joint_benefit_timeseries_vals,
                delta in selection_strength_timeseries_vals), # Use *_vals instead elements of axes_timeseries because axes rounds maximum_joint_benefit
                properties_dict_timeseries,
               )
    timeseries_statistics = Dataset(; Dict(:most_common_game_types => most_common_game_types,
                   :order_parameters => order_parameters,
                   :strategy_parity => strategy_parity,
                   :fraction_communicative => fraction_communicative,
                   :timeseries => timeseries,
                   )...)
    if isnothing(decimation_factor)
      savedataset(timeseries_statistics,
            path=datadir("processed","netcdf", "timeseries-statistics_matrixSource=$(adj_matrix_source)_timesteps=$time_steps.nc"),
            driver=:netcdf, overwrite=true, compress=9)
    else
      timeseries_statistics_decimated = timeseries_statistics[time_step = 1:decimation_factor:time_steps]
      savedataset(timeseries_statistics_decimated,
            path=datadir("processed","netcdf", "timeseries-statistics_decimationFactor=$(decimation_factor)_matrixSource=$(adj_matrix_source)_timesteps=$time_steps.nc"),
            driver=:netcdf, overwrite=true)
    end
  else
        throw(ArgumentError("data_type must be a string in set [\"cumulative\", \"timeseries-statistics\"]"))
  end

	return nothing
end
