using Symbolics
using LaTeXStrings
using DrWatson
using DataFramesMeta
using CairoMakie
using GraphMakie
using Colors
using ColorBrewer
using Polyhedra
using NetworkLayout

include("game_taxonomy.jl")
include("utils.jl")

const paired_colors = ColorBrewer.palette("Paired", 12)
const game_type_colors = Dict(harmony => paired_colors[7],
                              chicken => paired_colors[6], # exponential fixation time
                              battle => paired_colors[5],  # exponential fixation time
                              hero => paired_colors[8], # exponential fixation time
                              compromise => paired_colors[3],
                              concord => paired_colors[4],
                              staghunt => paired_colors[2],
                              dilemma => paired_colors[11],
                              deadlock => paired_colors[1],
                              assurance => paired_colors[10],
                              coordination => paired_colors[9],
                              peace => paired_colors[12],
                              neutral => :grey,
                             )

function plot_cumulative(selection_strength::Real, symmetry_breaking::Real,
                         adj_matrix_source::String="well-mixed",
                         payoff_update_method::String="single-update",
                         time_steps::Integer=2_000_000,
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 mutation_rate::Real=0.0001,
			 )

    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps,
                      selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = 20
    end
    data, _ = produce_or_load(calc_cumulative, config, datadir("raw","cumulative"))

    # Produce figure
    fig = generate_cumulative_plot(data, config)

    # Save figure
    filename = plotsdir("cumulative", savename(config, "png"))
    save(filename, fig)

    return fig
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
    beta0(B0) = unilateral_to_mutual_benefit*B0

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(config["adj_matrix_source"])
    graph = SimpleDiGraph(interaction_adj_matrix)
    nb_effective = ( mean(indegree(graph))+1) # Add one since n=degree+1 for well-mixed case

    lines!(ax, data["Bs"][begin] .. data["Bs"][end],
	   B0 -> analytic_frac_communicative(B0, beta0(B0);
			selection_strength=config["selection_strength"], cost=config["cost"], nb_players=nb_effective,
			symmetry_breaking=config["symmetry_breaking"], nb_phases=config["nb_phases"]); label="Theory",
           color=:orange)
    lines!(ax, data["Bs"][begin] .. data["Bs"][end],
           B0 -> 1 / (1 + exp(config["selection_strength"] * (nb_effective - 1) *
                              ((nb_effective - 1) * config["cost"]
			       + nb_effective * beta0(B0) * (1-2*config["symmetry_breaking"])/2
                              - (nb_effective - 2) / 2 * B0))); label="Approx. Theory",
           color=:purple, linestyle=:dash)


    # Add legend
    axislegend(ax; position=:lt)

    return fig
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
			 nb_phases::Integer=20,
			 cost::Real=0.1,
			 mutation_rate::Real=0.0001,
			 )
    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps, B_factor,
                      symmetry_breaking, selection_strength, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = 20
    end
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
                              nb_phases::Integer=20,
			      cost::Real=0.1,
			      mutation_rate::Real=0.0001,
			      )
    # Load results
    config = @strdict(adj_matrix_source, payoff_update_method, time_steps, B_factor,
                      selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = 20
    end
    data, _ = produce_or_load(calc_timeseries, config, datadir("raw","timeseries"))

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)

    # Create animation
    animation = generate_graph_evolution(data, graph)

    # Save animation
    filename = plotsdir("animations", savename(config, "gif"))
    save(filename, animation)

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
                                 concord => [R <= T, T <= S, S < P],
				 )
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

	non_neutral_games = [game for game in instances(GameType) if game != neutral]
        for game_type in non_neutral_games
            poly = polyhedron(intersect(inequality_to_hrep.(map(x -> x.val,
                                                                [game_type_inequalities(R,
                                                                                        S,
                                                                                        T,
                                                                                        P)[game_type]...,
                                                                 image_borders...]))...))
            if Polyhedra.volume(poly) == 0
                continue
            end
            mesh!(ax, Polyhedra.Mesh(poly); color=game_type_colors[game_type])
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
