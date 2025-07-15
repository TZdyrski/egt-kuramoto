### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 69e8f96d-55ad-4bf3-86d6-be9a0f448c41
import Pkg

# ╔═╡ d1ed094e-dd46-4199-916e-8fabcd24660c
# ╠═╡ show_logs = false
# This cell disables Pluto's package manager and activates the global environment. Click on ? inside the bubble next to Pkg.activate to learn more.
# (added automatically because a sysimage is used)
Pkg.activate("/home/jovyan/binder")

# ╔═╡ 50ccf7ba-55f7-11f0-1752-a5d12fed803c
using CairoMakie, WGLMakie, Downloads, YAXArrays, NetCDF, PlutoUI, LaTeXStrings, Graphs, StatsBase, XLSX, MetaGraphsNext, SimpleWeightedGraphs, NetworkLayout, GraphMakie, Colors, ColorBrewer

# ╔═╡ 3bb479c9-6dc1-4abf-88ed-f95b7f555972
md"# Graph Structure"

# ╔═╡ 4ba8db21-5ef0-44a1-9d58-efa99aa6db91
md"### Matrix Source"

# ╔═╡ 5c0fbece-68ba-4af5-81e9-459f0c70a383
@bind matrix_source Select([
	"well-mixed",
	"c-elegans",
	"c-elegans-undirected",
	"c-elegans-unweighted",
])

# ╔═╡ a9fdbf60-477e-4e55-a970-f3e12e04382f
md"## Strongly Connected Components"

# ╔═╡ 23c42826-463f-4b3a-b03f-5e9f5456898d
md"# Simulation Results"

# ╔═╡ c0a8e4f7-6e24-4f3c-b6e4-982f3f7794e7
md"### Selection Strength"

# ╔═╡ 1013bc9f-45ab-43b0-b121-aab2d12e284d
@bind selection_strength PlutoUI.Slider([0.005,0.2,5]; default=0.2)

# ╔═╡ dabbdc97-cade-4643-9b6e-a2b30c1c7639
selection_strength

# ╔═╡ af6fc3f6-cb76-49ca-ad14-1ac18121d477
md"### Symmetry Breaking"

# ╔═╡ 26ea5ad3-76e3-4f58-a441-707fc9995be4
@bind symmetry_breaking PlutoUI.Slider(0:0.25:1; default=0.5)

# ╔═╡ f9dc48b9-3654-4570-b0b4-014cde8ffe37
symmetry_breaking

# ╔═╡ 64b0ba25-8a2e-427e-baf0-984389802336
md"## Communicative Fraction"

# ╔═╡ 974538d6-19da-4251-8441-a0ff70a85e39
md"# Time Dependent Features"

# ╔═╡ 9166474d-67ae-47d2-9bf0-04eda2fcce3d
md"### Time step"

# ╔═╡ 60a09648-d129-4005-9d06-df6844d28955
md"### Maximum Joint Benefit"

# ╔═╡ cb324cb3-48ea-4a3a-b4c7-cf3d84a9030a
@bind maximum_joint_benefit PlutoUI.Slider([0.0005,0.15,0.25]; default=0.15)

# ╔═╡ 507c9c3f-7fd9-4f1a-ad7d-743697a5c580
maximum_joint_benefit

# ╔═╡ fa0030dd-4cb3-4edb-9067-119552d3d0ab
md"## Time snapshot"

# ╔═╡ 4d0adc28-0eaa-481b-934a-a4344a638c22
md"## Time Series"

# ╔═╡ 74c12270-0058-4c09-817a-ac46f799518f
md"## Game Types"

# ╔═╡ 14ded9ce-8568-4df9-8def-1ce144efedbf
md"## Time Evolution Animation"

# ╔═╡ b6d6c690-46d6-44fc-b627-b66a1c8290e2
# Use latex theme in figures
set_theme!(theme_latexfonts())

# ╔═╡ ee4e891d-b61b-424a-8542-2d6dbada13d1
# Download data
begin
	function open_or_download(filename::String,url::String)
		if isfile(filename)
			# If file is already downloaded to current directory, use it
			file = filename
		else
			# Download file to current directory
			file = Downloads.download(url, filename)
		end
		return file
	end

	local cumulativeDataUrls = Dict(
		"well-mixed" => "https://www.dropbox.com/scl/fi/fz313w3d25u5e4na7ynv8/cumulative_matrixSource-well-mixed_timesteps-200000000.nc?rlkey=ybxwciz5espqibpeyasjn0ob7&st=qy3yag25&dl=0",
		"c-elegans" => "https://www.dropbox.com/scl/fi/jhfucubszew1qcrmxl7ng/cumulative_matrixSource-c-elegans_timesteps-200000000.nc?rlkey=0bbngxqr7mgwjrrj3l8cz3b65&st=rsrtax2g&dl=0",
		"c-elegans-undirected" => "https://www.dropbox.com/scl/fi/rku7qod93ap3exxtwhnh2/cumulative_matrixSource-c-elegans-undirected_timesteps-200000000.nc?rlkey=e80jmmbb24658yew7l8w2va9e&st=rn5fzbap&dl=0",
		"c-elegans-unweighted" => "https://www.dropbox.com/scl/fi/40hz3gd8f9lzpfihj3ggp/cumulative_matrixSource-c-elegans-unweighted_timesteps-200000000.nc?rlkey=cnd7peqz3sk00dzirjryp1thm&st=owps8tza&dl=0",
	)
	local cumulativeDataFilenames = Dict(
		"well-mixed" => "cumulative_matrixSource=well-mixed_timesteps=200000000.nc",
		"c-elegans" => "cumulative_matrixSource=c-elegans_timesteps=200000000.nc",
		"c-elegans-undirected" => "cumulative_matrixSource=c-elegans-undirected_timesteps=200000000.nc",
		"c-elegans-unweighted" => "cumulative_matrixSource=c-elegans-unweighted_timesteps=200000000.nc",
	)
	cumulativeDatasets = Dict(k => open_or_download(cumulativeDataFilenames[k],v) for (k,v) in cumulativeDataUrls)

	local timeseriesDataUrls = Dict(
		"well-mixed" => "https://www.dropbox.com/scl/fi/0jn3ysnbfenfjc5wccrfi/timeseries-statistics_decimationFactor-1000_matrixSource-well-mixed_timesteps-800000.nc?rlkey=blznt21inkcxvcuco0rp47t5j&st=dhea3wdp&dl=0",
		"c-elegans" => "https://www.dropbox.com/scl/fi/8zs03ogez74nziz8koe6c/timeseries-statistics_decimationFactor-1000_matrixSource-c-elegans_timesteps-800000.nc?rlkey=cu1427hlxen620786msbw4bas&st=jgkn63qh&dl=0",
		"c-elegans-undirected" => "https://www.dropbox.com/scl/fi/x5gwv3eawng2nv2bbvxee/timeseries-statistics_decimationFactor-1000_matrixSource-c-elegans-undirected_timesteps-800000.nc?rlkey=iv75xqtiegnsre6pyjk0txuzv&st=eme7yl67&dl=0",
		"c-elegans-unweighted" => "https://www.dropbox.com/scl/fi/8p14eq33izckbpqdr5f5v/timeseries-statistics_decimationFactor-1000_matrixSource-c-elegans-unweighted_timesteps-800000.nc?rlkey=5jrc1yv2h01596ebjlpr3bjls&st=3mrmq94v&dl=0",
	)
	local timeseriesDataFilenames = Dict(
		"well-mixed" => "timeseries-statistics_decimationFactor=1000_matrixSource=well-mixed_timesteps=800000.nc",
		"c-elegans" => "timeseries-statistics_decimationFactor=1000_matrixSource=c-elegans_timesteps=800000.nc",
		"c-elegans-undirected" => "timeseries-statistics_decimationFactor=1000_matrixSource=c-elegans-undirected_timesteps=800000.nc",
		"c-elegans-unweighted" => "timeseries-statistics_decimationFactor=1000_matrixSource=c-elegans-unweighted_timesteps=800000.nc",
	)
	timeseriesDatasets = Dict(k => open_or_download(timeseriesDataFilenames[k],v) for (k,v) in timeseriesDataUrls)
end;

# ╔═╡ 86f70955-e3a7-402c-9a2a-945271d667be
begin
	cumulative_da = open_dataset(cumulativeDatasets[matrix_source], driver=:netcdf).layer
	timeseries_ds = open_dataset(timeseriesDatasets[matrix_source], driver=:netcdf)

	nb_phases = cumulative_da.properties["nb_phases"]
	cost = cumulative_da.properties["cost"]
end;

# ╔═╡ 074771dd-e442-484d-94b1-8d023d7ae74d
@bind time_step PlutoUI.Slider(collect(timeseries_ds.time_step); default=554000)

# ╔═╡ 57277407-6b80-40a5-87de-ef026c8b2c04
time_step

# ╔═╡ a0682454-e0f1-4e8c-b4fc-0066facbd79a
# Get graphs
begin
	local connectomeUrl = "https://github.com/openworm/ConnectomeToolbox/raw/05c50d5d197bdb7b1607c480a4f9e4aef9287232/cect/data/SI%205%20Connectome%20adjacency%20matrices.xlsx"
	local connectomeFilename = "SI 5 Connectome adjacency matrices.xlsx"
	local connectomeXlsx = open_or_download(connectomeFilename,connectomeUrl)
	local connectomeAndMuscleTable = XLSX.readxlsx(connectomeXlsx)["hermaphrodite chemical"]["C3:QN303"]
	local connectomeTable = connectomeAndMuscleTable[:, [1:21..., 52:323..., 439:446...]]

	local row_labels = connectomeTable[2:end,1]
	local col_labels = connectomeTable[1,2:end]
	# Verify that the row and column labels match
	if row_labels != col_labels
		throw(ErrorException("Rows and columns in c-elegans adjacency table do not match up"))
	end

	local connectomeAdjMatrix = coalesce.(connectomeTable[2:end,2:end],0)
	local weightedDiGraph = SimpleWeightedDiGraph(connectomeAdjMatrix)
	local weightedGraph = SimpleWeightedGraph(connectomeAdjMatrix
			+ transpose(connectomeAdjMatrix))
	local connectomeGraph = MetaGraph(
		weightedDiGraph,
		[row_labels[v] => nothing for v in collect(vertices(weightedDiGraph))],
		[(row_labels[e.src], row_labels[e.dst]) => e.weight for e in  edges(weightedDiGraph)],
		nothing,
		identity,
	)
	local connectomeGraphUndirected = MetaGraph(
		weightedGraph,
		[row_labels[v] => nothing for v in collect(vertices(weightedGraph))],
		[(row_labels[e.src], row_labels[e.dst]) => e.weight for e in  edges(weightedGraph)],
		nothing,
		identity,
	)
	local connectomeGraphUnweighted = MetaGraph(
		weightedDiGraph,
		[row_labels[v] => nothing for v in collect(vertices(weightedDiGraph))],
		[(row_labels[e.src], row_labels[e.dst]) => nothing for e in  edges(weightedDiGraph)],
	)

	function get_graph(adj_matrix_source::String, properties::Dict)
		if adj_matrix_source == "well-mixed"
			return complete_digraph(properties["nb_players"])
		elseif adj_matrix_source == "c-elegans"
			return connectomeGraph
		elseif adj_matrix_source == "c-elegans-undirected"
			return connectomeGraphUndirected
		elseif adj_matrix_source == "c-elegans-unweighted"
			return connectomeGraphUnweighted
		end
	end

	# Precalculate strongly connected components
	local connectomeConnectedComponents = strongly_connected_components(connectomeGraph)

	function get_strongly_connected_component(adj_matrix_source::String, properties::Dict)
		if adj_matrix_source == "well-mixed"
			return strongly_connected_components(complete_digraph(properties["nb_players"]))
		elseif adj_matrix_source == "c-elegans"
			return connectomeConnectedComponents
		elseif adj_matrix_source == "c-elegans-undirected"
			return connectomeConnectedComponents # We don't use the weights or directedness in this notebook
		elseif adj_matrix_source == "c-elegans-unweighted"
			return connectomeConnectedComponents # We don't use the weights or directedness in this notebook
		end
	end
end;

# ╔═╡ 1051a191-22b3-4161-9ff9-96f1ed27c68f
graph = get_graph(matrix_source, cumulative_da.properties);

# ╔═╡ a50ec494-b320-4912-9fbe-61fdabd7230c
begin
    local conn_comp = get_strongly_connected_component(matrix_source, cumulative_da.properties);

    # Calculate connected components
    local conn_comp_index = zeros(Int32, nv(graph))
    for (conn_comp_indx, elem) in pairs(conn_comp)
        for indx in elem
            conn_comp_index[indx] = conn_comp_indx
        end
    end
    local num_conn_comp = length(conn_comp)

    # Plot graph
    local fig = graphplot(graph;
					layout=Stress(),
					edge_color=(:black, 0.05),
					node_color=distinguishable_colors(num_conn_comp)[conn_comp_index],
					edge_plottype=:linesegments)
	fig
end

# ╔═╡ 2f94f5b7-7c43-4083-93d5-d21e6c9a31c9
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
end;

# ╔═╡ d2829848-20fe-4f52-a71b-976ca9d0131a
begin
	# Plot
    local fig = Figure()
    local ax = Axis(fig[1, 1];
		title = L"Selection $\delta = %$selection_strength$",
	    xlabel = L"Maximum benefit of mutual communication, $B(0)$",
	    ylabel = "Frequency of communicative strategies",
	    limits = (nothing, nothing, 0-0.05f0, 1+0.05f0),
	)
    scatter!(ax, cumulative_da[
		selection_strength = At(selection_strength),
		symmetry_breaking = At(symmetry_breaking)];
			 label="Simulation")

	if matrix_source == "well-mixed"
		# Add one since n=degree+1 for well-mixed case
		local nb_effective = mean(indegree(graph))+1
	    local beta0 = B0 -> 0.95*B0
		Bs = lookup(cumulative_da, :maximum_joint_benefit)
		lines!(ax, minimum(Bs)..maximum(Bs),
			   B0 -> analytic_frac_communicative(B0, beta0(B0);
					selection_strength, cost, nb_players=nb_effective,
					symmetry_breaking, nb_phases);
			   label="Theory",
	           color=:orange)
	end

    # Add legend
    axislegend(ax; position=:lt)

	# Return figure
    fig
end

# ╔═╡ eaa6bf62-daea-450c-85f4-b4dd7ad4048b
begin
	@enum GameType harmony chicken battle hero compromise concord staghunt dilemma deadlock assurance coordination peace
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
	                              peace => paired_colors[12]
	                             )
	const enum_to_strategy = Dict(parse(Int,match(r"(?<=enum_lookup:)(\d+)",k).match)=> eval(Symbol(v)) for (k,v) in timeseries_ds.most_common_game_types.properties
		if startswith(k,"enum_lookup"))

	@enum StrategyParity all_communicative all_noncommunicative mixed
    const strategy_parity_colors = Dict(all_communicative => :lightgrey,
										all_noncommunicative => :darkgrey)
	const enum_to_parity = Dict(parse(Int,match(r"(?<=enum_lookup:)(\d+)",k).match)=> eval(Symbol(v)) for (k,v) in timeseries_ds.strategy_parity.properties
		if startswith(k,"enum_lookup"))

	const game_type_or_parity_colors = merge(game_type_colors, strategy_parity_colors)
end;

# ╔═╡ 0b9e9e44-6c46-4481-bf89-07cd11c76adb
begin
    # Plot fraction communicative
    local fig = Figure()
    local ax1 = Axis(fig[1, 1];
               title=L"Strong selection $\delta = 0.2$",
               xlabel="Time",
               ylabel="Frequency of communicative strategies",
			   yticklabelcolor=:blue,
			   ylabelcolor=:blue,
               limits=(nothing, nothing, -0.05, 1.05))

	local ds = timeseries_ds[
		selection_strength = At(selection_strength),
		maximum_joint_benefit = At(maximum_joint_benefit),
		symmetry_breaking = At(symmetry_breaking),
	]

	if ! any(ismissing.(ds.most_common_game_types))
		game_parity_or_type = [enum_to_parity[strategy_parity] != mixed ? enum_to_parity[strategy_parity] : enum_to_strategy[game_type]
	                           for (strategy_parity, game_type) in
	                               zip(ds.strategy_parity,
									   ds.most_common_game_types)]
		colors = [game_type_or_parity_colors[key] for key in game_parity_or_type]
	    local li1 = lines!(ax1,
						   lookup(ds.fraction_communicative, :time_step),
						   collect(ds.fraction_communicative);
						   color=collect(colors),
						   linewidth=5,
						  )

	    # Plot order parameter
	    local ax2 = Axis(fig[1, 1];
	               ylabel="Order parameter",
	               limits=(nothing, nothing, -0.05, 1.05),
	               yaxisposition=:right,
	               yticklabelcolor=:magenta,
				   ylabelcolor=:magenta)
	    hidespines!(ax2)
	    hidexdecorations!(ax2)
	    local li2 = lines!(ax2,
						   lookup(ds.order_parameters, :time_step),
						   collect(ds.order_parameters);
						   color=:magenta)

	    # Add legend
		local game_list = unique(game_parity_or_type)
		local legend_items = [PolyElement(color=game_type_or_parity_colors[game]) for game in game_list]
		local name_list = [String(Symbol(x)) for x in game_list]
		replace!(name_list, "all_communicative" => "all-C",
								"all_noncommunicative" => "all-N")
		append!(legend_items, [PolyElement(color=:magenta)])
		append!(name_list, ["Order parameter"])
	    axislegend(ax1, legend_items, name_list;
	               position=:rb)

		# Return figure
		fig
	else
		md"No game-type data for this set of parameters"
	end
end

# ╔═╡ b868b861-2962-4a8b-af81-bf85838b1ab9
begin
	local ds = timeseries_ds[
		selection_strength = At(selection_strength),
		maximum_joint_benefit = At(maximum_joint_benefit),
		symmetry_breaking = At(symmetry_breaking),
	]

	if !any(ismissing.(ds.most_common_game_types))
		# Plot histogram of game types
	    local game_parity_or_type = [enum_to_parity[strategy_parity] != mixed ? enum_to_parity[strategy_parity] : enum_to_strategy[game_type]
	                           for (strategy_parity, game_type) in
	                               zip(ds.strategy_parity,
									   ds.most_common_game_types)]
	    local hist_data = proportionmap(game_parity_or_type)

	    local fig = Figure()
		local name_list = String.(Symbol.(collect(keys(hist_data))))
		replace!(name_list, "all_communicative" => "all-C",
								"all_noncommunicative" => "all-N")
	    local ax = Axis(fig[1,1];
	               title=L"Strong selection $\delta = 0.2$",
	               limits=(nothing, nothing, 0, nothing),
	               xticks=(1:length(keys(hist_data)), name_list))
	    barplot!(ax, collect(values(hist_data));
				color = [game_type_or_parity_colors[game] for game in keys(hist_data)])

		# Return figure
		fig
	else
		md"No game-type data for this set of parameters"
	end
end

# ╔═╡ 0e599cc9-42a9-4c32-b441-4ef15346f24e
begin
	local ds = timeseries_ds[
		selection_strength = At(selection_strength),
		maximum_joint_benefit = At(maximum_joint_benefit),
	]

	if ! any(ismissing.(ds.most_common_game_types))
		local game_parity_or_type = [enum_to_parity[strategy_parity] != mixed ? enum_to_parity[strategy_parity] : enum_to_strategy[game_type]
	                           for (strategy_parity, game_type) in
	                               zip(ds.strategy_parity,
									   ds.most_common_game_types)]
		local game_type_distribution = [proportionmap(slice) for slice in eachslice(game_parity_or_type, dims=timeseries_ds.symmetry_breaking)]

	    # Generate plot
	    local fig = Figure()
	    local ax = Axis(fig[1, 1];
						xlabel = L"Asymmetry $\alpha$",
						title = "Proportion of Game Types by Asymmetry")
	    for (indx, row) in enumerate(game_type_distribution)
	        local asymm = ds.symmetry_breaking[indx]
	        barplot!(repeat([asymm],length(row)),
					 collect(values(row)),
					 stack = repeat([indx], length(row)),
					 color = [game_type_or_parity_colors[game] for game in keys(row)],
					 width=0.2)
	    end

	    # Add legend
		local game_list = unique(collect(Iterators.flatten(map(keys,game_type_distribution))))
		local legend_items = [PolyElement(color=game_type_or_parity_colors[game]) for game in game_list]
		local name_list = [String(Symbol(x)) for x in game_list]
		replace!(name_list, "all_communicative" => "all-C",
								"all_noncommunicative" => "all-N")
		Legend(fig[2,1], legend_items,name_list, orientation=:horizontal)

		# Return figure
		fig
	else
		md"No game-type data for this set of parameters"
	end
end

# ╔═╡ baad79d8-539e-4e71-8001-06809527202c
begin
	# Create colormap
	cooperative_colors = range(colorant"navyblue";
							   stop=colorant"paleturquoise1",
							   length=nb_phases)
	noncooperative_colors = range(colorant"darkred";
								  stop=colorant"lightsalmon1",
								  length=nb_phases)
	cooperative_colormap = [cooperative_colors; noncooperative_colors]
end;

# ╔═╡ 16b202df-7050-4505-8a4a-2f32f6ade20b
begin
	# Get timeseries data
	local timeseries_data = timeseries_ds.timeseries[
		time_step = At(time_step),
		selection_strength = At(selection_strength),
		maximum_joint_benefit = At(maximum_joint_benefit),
		symmetry_breaking = At(symmetry_breaking),
	]

	if ! any(ismissing.(timeseries_data))
		# Apply colormap
		local colors = cooperative_colormap[timeseries_data]

		# Create plot
		GC.gc()
		local fig = graphplot(graph;
				  node_color=colors,
				  layout=Stress(),
				  arrow_show=false, edge_color=(:black, 0.05),
				  edge_plottype=:linesegments)
		GC.gc()
		fig
	end
end

# ╔═╡ 54ba97a0-7da0-4227-8427-f54701df0c5c
# ╠═╡ disabled = true
#=╠═╡
begin
	# Get timeseries data
	local timeseries_data = timeseries_ds.timeseries[
		selection_strength = At(selection_strength),
		maximum_joint_benefit = At(maximum_joint_benefit),
		symmetry_breaking = At(symmetry_breaking),
	]

	if ! any(ismissing.(timeseries_data))
		# It seems WGLMakie animations aren't working yet
		CairoMakie.activate!()

		# Apply colormap
		local colors = cooperative_colormap[timeseries_data]

	    # Set animation parameters
	    local num_times = size(timeseries_data, :time_step)
	    local total_time_s = 10
	    local framerate_Hz = 30
	    local time_stride = round(num_times / (framerate_Hz * total_time_s))

	    # Generate animation
	    local time = Observable(1)
	    local layout = Stress()
		local titleStr = @lift("Time step = "*string(timeseries_data.time_step[$time]))

	    # Create plot
	    local color_observable = @lift(colors[:, $time])
		local fig = Figure()
	    local ax = Axis(fig[1,1];
						title=titleStr)
	    local graph_plot = graphplot!(ax, graph;
										node_color=color_observable,
										layout=layout,
	                                    arrow_show=false,
										edge_color=(:black, 0.05),
	                                    edge_plottype=:linesegments)

	    recording = Record(fig, 1:time_stride:num_times; framerate=framerate_Hz) do t
	        return time[] = t
	    end
		WGLMakie.activate!()

		# Return recording
		recording
	else
		md"No timeseries data for this set of parameters"
	end
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─3bb479c9-6dc1-4abf-88ed-f95b7f555972
# ╟─4ba8db21-5ef0-44a1-9d58-efa99aa6db91
# ╟─5c0fbece-68ba-4af5-81e9-459f0c70a383
# ╟─a9fdbf60-477e-4e55-a970-f3e12e04382f
# ╟─a50ec494-b320-4912-9fbe-61fdabd7230c
# ╟─23c42826-463f-4b3a-b03f-5e9f5456898d
# ╟─c0a8e4f7-6e24-4f3c-b6e4-982f3f7794e7
# ╟─1013bc9f-45ab-43b0-b121-aab2d12e284d
# ╟─dabbdc97-cade-4643-9b6e-a2b30c1c7639
# ╟─af6fc3f6-cb76-49ca-ad14-1ac18121d477
# ╟─26ea5ad3-76e3-4f58-a441-707fc9995be4
# ╟─f9dc48b9-3654-4570-b0b4-014cde8ffe37
# ╟─64b0ba25-8a2e-427e-baf0-984389802336
# ╟─d2829848-20fe-4f52-a71b-976ca9d0131a
# ╟─974538d6-19da-4251-8441-a0ff70a85e39
# ╟─9166474d-67ae-47d2-9bf0-04eda2fcce3d
# ╟─074771dd-e442-484d-94b1-8d023d7ae74d
# ╟─57277407-6b80-40a5-87de-ef026c8b2c04
# ╟─60a09648-d129-4005-9d06-df6844d28955
# ╟─cb324cb3-48ea-4a3a-b4c7-cf3d84a9030a
# ╟─507c9c3f-7fd9-4f1a-ad7d-743697a5c580
# ╟─fa0030dd-4cb3-4edb-9067-119552d3d0ab
# ╟─16b202df-7050-4505-8a4a-2f32f6ade20b
# ╟─4d0adc28-0eaa-481b-934a-a4344a638c22
# ╟─0b9e9e44-6c46-4481-bf89-07cd11c76adb
# ╟─74c12270-0058-4c09-817a-ac46f799518f
# ╟─b868b861-2962-4a8b-af81-bf85838b1ab9
# ╟─0e599cc9-42a9-4c32-b441-4ef15346f24e
# ╟─14ded9ce-8568-4df9-8def-1ce144efedbf
# ╟─54ba97a0-7da0-4227-8427-f54701df0c5c
# ╟─69e8f96d-55ad-4bf3-86d6-be9a0f448c41
# ╟─d1ed094e-dd46-4199-916e-8fabcd24660c
# ╟─50ccf7ba-55f7-11f0-1752-a5d12fed803c
# ╟─b6d6c690-46d6-44fc-b627-b66a1c8290e2
# ╟─ee4e891d-b61b-424a-8542-2d6dbada13d1
# ╟─86f70955-e3a7-402c-9a2a-945271d667be
# ╟─a0682454-e0f1-4e8c-b4fc-0066facbd79a
# ╟─1051a191-22b3-4161-9ff9-96f1ed27c68f
# ╟─2f94f5b7-7c43-4083-93d5-d21e6c9a31c9
# ╟─eaa6bf62-daea-450c-85f4-b4dd7ad4048b
# ╟─baad79d8-539e-4e71-8001-06809527202c
