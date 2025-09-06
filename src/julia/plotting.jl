using DrWatson
using CairoMakie
using Graphs
using GraphMakie
using Colors
using NetworkLayout

include("utils.jl")

function plot_graph_evolution(;B_to_c::Real, selection_strength::Real,
                              symmetry_breaking::Real,
                              adj_matrix_source::String="well-mixed",
                              time_steps::Integer=80_000,
                              nb_phases::Integer=20,
                              cost::Real=0.1,
                              beta_to_B::Real=0.95,
                              mutation_rate::Real=0.0001,
                              nb_players::Integer=20,
			      )
    # Load results
    config = @strdict(adj_matrix_source, time_steps, B_to_c, beta_to_B,
                      selection_strength, symmetry_breaking, nb_phases, cost, mutation_rate)
    if adj_matrix_source == "well-mixed" || adj_matrix_source == "random-regular-graph" || adj_matrix_source == "random-regular-digraph"
	    config["nb_players"] = nb_players
    end
    data = wload(datadir("raw", "timeseries", savename(config,"jld2")))

    # Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(; adj_matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)

    # Create animation
    animation = generate_graph_evolution(data, graph)

    # Save animation
    filename = plotsdir("animations", savename(config, "mp4"))
    save(filename, animation)

end

function generate_graph_evolution(data, graph::Graphs.SimpleGraphs.AbstractSimpleGraph)
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

    # Create plot
    layout = Stress()
    color_observable = @lift(colors[:, $time])
    fig, _ = graphplot(graph; node_color=color_observable, layout=layout,
                                    arrow_show=false, edge_color=(:black, 0.05),
                                    edge_plottype=:linesegments)

    recording = CairoMakie.Makie.Record(fig, 1:time_stride:num_times; framerate=framerate_Hz) do t
        return time[] = t
    end
    return recording
end
