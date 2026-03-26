# SPDX-License-Identifier: GPL-3.0-or-later

using CairoMakie
using Graphs
using GraphMakie
using Colors
using NetworkLayout

include("utils.jl")

function plot_graph_evolution(data::Dict, interaction_graph::AbstractGraph)
		# Decode data
		time_steps = size(data["deltas"])[2]
    data["all_populations"] = decode_delta_encoded_all(data["initial_actions"], data["deltas"], time_steps)

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
		titleStr = @lift("Time step = "*string($time))
    fig = Figure()
    ax = Axis(fig[1,1]; title=titleStr)
    graphplot!(ax, interaction_graph; node_color=color_observable, layout=layout,
                                    arrow_show=false, edge_color=(:black, 0.05),
                                    edge_plottype=:linesegments)

    recording = CairoMakie.Makie.Record(fig, 1:time_stride:num_times; framerate=framerate_Hz) do t
        return time[] = t
    end
    return recording

end
