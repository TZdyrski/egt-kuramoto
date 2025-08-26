using DataToolkit
using LinearAlgebra

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
    elseif adj_matrix_source == "drosophilia"
        interaction_adj_matrix = round.(get_drosophilia_connectome())
        reproduction_adj_matrix = interaction_adj_matrix + I
    else
        throw(ArgumentError("adj_matrix_source must be a string in set [\"well-mixed\", "
                            * "\"c-elegans\", \"c-elegans-unweighted\", "
                            * "\"c-elegans-undirected\", \"c-elegans-undirected-unweighted\", "
                            * "\"drosophilia\", "
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

function get_drosophilia_connectome()
    connectome = get_drosophilia_connectome_labelled()
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

function get_drosophilia_connectome_labelled()
    # Get list of neurons
    classification = read(dataset("drosophilia-classification"), DataFrame)
    # Add index as a column
    classification.index = 1:nrow(classification)
    # Create dictionary from root ID to index
    root_id_dictionary = Dict(eachrow(classification[!,["root_id", "index"]]))
    # Get adjacency list
    adjacency_list_full = read(dataset("drosophilia-connectome"), DataFrame)
    # Replace root IDs with indices
    adjacency_list = transform(adjacency_list_full[!,["pre_root_id","post_root_id","syn_count"]],
			       :pre_root_id => ByRow(x -> root_id_dictionary[x]) => :pre_root_id,
			       :post_root_id => ByRow(x -> root_id_dictionary[x]) => :post_root_id)
    # Convert adjacency list to weighted digraph
    connectome_graph = SimpleWeightedDiGraph(adjacency_list.pre_root_id, adjacency_list.post_root_id,
				       adjacency_list.syn_count)
    # Convert weighted digraph to adjacency matrix
    connectome = adjacency_matrix(connectome_graph)
    return connectome
end
