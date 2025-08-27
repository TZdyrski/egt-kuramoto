using DrWatson
using BlockArrays

include("moran.jl")
include("postprocess.jl")

unilateral_to_mutual_benefit = 0.95 # \beta_0/B_0

function payoff_matrix(nb_phases::Integer,
                       mutual_benefit_synchronous::Real,
                       unilateral_benefit_synchronous::Real, cost::Real;
                       symmetry_breaking::Real=1 / 2)
    benefit_scaling = [(1 + cos(2 * pi * (phi_j - phi_i))) / 2
                       for phi_i in (0:(nb_phases - 1)) / nb_phases,
                           phi_j in (0:(nb_phases - 1)) / nb_phases]
    payoff_matrix = Matrix(mortar([[mutual_benefit_synchronous*ones(size(benefit_scaling)) .- cost,
                                    2 * (1 - symmetry_breaking) *
                                    unilateral_benefit_synchronous * benefit_scaling];;
                                   [2 * symmetry_breaking * unilateral_benefit_synchronous *
                                    benefit_scaling .- cost,
                                    zeros(eltype(benefit_scaling), size(benefit_scaling))]]) .+
                           cost)
    return payoff_matrix
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

function calc_cumulative(config::Dict)
    # Unpack values
    @unpack selection_strength, symmetry_breaking, adj_matrix_source, payoff_update_method, time_steps, nb_phases, cost, mutation_rate = config

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
    @time cumulative_populations = [cumulative(Moran(NormalFormGame(payoff_matrix(nb_phases,
                                                                                                  B,
                                                                                                  unilateral_to_mutual_benefit *
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
    return @strdict(Bs, nb_players, fraction_communicative)
end

function calc_timeseries(config::Dict)
    # Unpack variables
    @unpack B_to_c, selection_strength, symmetry_breaking, adj_matrix_source, payoff_update_method, time_steps, nb_phases, cost, mutation_rate = config

    # Define system
    B = cost * B_to_c

    # Define interaction graph and reproduction graphs
    interaction_adj_matrix, reproduction_adj_matrix = get_adj_matrices(adj_matrix_source)
    # Specify number of players
    nb_players = size(interaction_adj_matrix)[1]

    # Run the model for weak selection strength
    all_populations, steps_following_mutation = time_series(Moran(NormalFormGame(payoff_matrix(nb_phases,
                                                                                     B,
                                                                                     unilateral_to_mutual_benefit *
                                                                                     B,
                                                                                     cost;
                                                                                     symmetry_breaking)),
                                                        interaction_adj_matrix,
                                                        reproduction_adj_matrix,
                                                        selection_strength, mutation_rate),
                                  time_steps, payoff_update_method)

    # Package results
    return @strdict(all_populations, steps_following_mutation, nb_phases, nb_players, interaction_adj_matrix)
end

