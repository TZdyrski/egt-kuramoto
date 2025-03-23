### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 3e050942-505f-4328-be97-d6e72b260be2
using DrWatson

# ╔═╡ f10a29cf-c5ff-4e1d-8864-6efab1634b48
# ╠═╡ show_logs = false
@quickactivate "2024_EGT_Kuramoto"

# ╔═╡ 0ecb3c44-a2d0-11ef-1d79-d3e37535b2dc
begin
	using CairoMakie
	using PlutoUI
	using DataFrames
	using DataFramesMeta
	using SimplePlutoInclude
end

# ╔═╡ fdf8f135-7dda-45ef-b8f1-08276e1e0375
@plutoinclude srcdir("julia", "local_moran_interaction.jl")

# ╔═╡ f5d72e8a-c624-411d-8b8f-a721c44d401f
md"""
# Communicative Fraction vs. $B(0)$
"""

# ╔═╡ 93340069-644f-4da1-b0af-824bf61a77ca
# ╠═╡ show_logs = false
begin
	time_steps_cumulative = Int(2E8)
	
	df_raw_cumulative = collect_results(datadir("raw","cumulative"); rinclude = [Regex("time_steps=$time_steps_cumulative[._]")])
	
	df_cumulative = transform(df_raw_cumulative, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)
end;

# ╔═╡ b9040b04-c6a4-428f-9ac2-518760c46ce2
begin
	selection_strength_bond = @bind selection_strength PlutoUI.Slider([0.005,0.2,5]);
	symmetry_breaking_bond = @bind symmetry_breaking PlutoUI.Slider(0:0.25:1);
	matrix_source_bond = @bind matrix_source Select(["well-mixed", "c-elegans-undirected-unweighted", "c-elegans-undirected", "c-elegans-unweighted", "c-elegans", "random-regular-graph", "random-regular-digraph"]);
	B_factor_bond = @bind B_factor PlutoUI.Slider([1.5,2.5]);
end;

# ╔═╡ 49e07eb6-c2fc-4223-85a9-aa9c986e34da
begin
	using Graphs
	using StatsBase
	# Generate graph
    interaction_adj_matrix, _ = get_adj_matrices(matrix_source)
    graph = SimpleDiGraph(interaction_adj_matrix)
    nb_effective =  mean(indegree(graph))+1 # Add one since n=degree+1 for well-mixed case
end;

# ╔═╡ 1740d8de-7c44-4a69-8711-5c25149c560a
md"""
### Selection Strength
"""

# ╔═╡ 26864283-0e1d-44e9-819f-44fc3e4c3aec
selection_strength_bond

# ╔═╡ 7c019aea-78ee-46b5-b9ee-3489dfa2570f
selection_strength

# ╔═╡ 5b6bdcd4-f608-4ea9-baec-d0a62d841bc0
md"""
### Symmetry Breaking
"""

# ╔═╡ adbf6ba7-420c-4d03-a92d-29f1c85929ca
symmetry_breaking_bond

# ╔═╡ ebf8ec31-a391-4a6c-864f-efd3f732d0c2
symmetry_breaking

# ╔═╡ 7e0e13aa-524f-4137-a6f6-3d4ee46a420e
md"""
### Matrix Source
"""

# ╔═╡ 3cc72edb-1a4b-4ac6-a18f-2ad9296b2da2
matrix_source_bond

# ╔═╡ 2f8d5984-acd3-4958-afe6-acc136d1c70a
begin
	df_cumulative_selected = @rsubset(df_cumulative, :selection_strength == selection_strength, :matrix_source == matrix_source, :symmetry_breaking == symmetry_breaking)

	dict_cumulative_selected = Dict()
	if nrow(df_cumulative_selected) == 1
		dict_cumulative_selected = Dict(names(df_cumulative_selected[1,:]) .=> values(df_cumulative_selected[1,:]))

		nothing
	elseif nrow(df_cumulative_selected) < 1
		println("No datasets found")

		println("Available datasets:")
		df_cumulative
	else
		println("More than one dataset matching conditions found: please fully specify")
		df_cumulative_selected
	end
end

# ╔═╡ 9e64dcc2-574e-4690-a10d-9b45e0ba5219
!isempty(dict_cumulative_selected) ? local_moran_interaction.generate_cumulative_plot(dict_cumulative_selected, @strdict(selection_strength,symmetry_breaking,nb_phases=20,adj_matrix_source=matrix_source)) : nothing

# ╔═╡ 0ca5e4f8-faab-4057-a823-956c3711b5d8
md"""
### Strongly Connected Components
"""

# ╔═╡ 2dabcbb3-7dc6-4bc4-af4f-8e48407130f8
local_moran_interaction.plot_connected_components(matrix_source)

# ╔═╡ 32510bf3-a932-4394-a4fb-8101580a6140
md"""
# Time Series
"""

# ╔═╡ 656dd0d7-4e47-422c-a104-9b2473a51306
# ╠═╡ show_logs = false
begin
	time_steps_timeseries = Int(8E5)
	df_timeseries_statistics = local_moran_interaction.load_all_timeseries_statistics(time_steps_timeseries)
end;

# ╔═╡ 4fefdb80-5168-44e3-ae7e-e92ef29a30e5
# ╠═╡ show_logs = false
begin
	df_raw_timeseries = collect_results(datadir("raw","timeseries"); rinclude = [Regex("time_steps=$time_steps_timeseries[._]")])

	df_timeseries = transform(df_raw_timeseries, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)
end;

# ╔═╡ abd83b14-d699-42cc-b0f0-6c1703c94ed3
begin
	df_timeseries_selected = @rsubset(df_timeseries, :selection_strength == selection_strength, :adj_matrix_source == matrix_source, :symmetry_breaking == symmetry_breaking)

	dict_timeseries_selected = Dict()
	if nrow(df_timeseries_selected) == 1
		dict_timeseries_selected = Dict(names(df_timeseries_selected[1,:]) .=> values(df_timeseries_selected[1,:]))

		nothing
	elseif nrow(df_timeseries_selected) < 1
		println("No datasets found")

		println("Available datasets:")
		df_timeseries
	else
		println("More than one dataset matching conditions found: please fully specify")
		df_timeseries_selected
	end
end

# ╔═╡ e8890847-13d7-447c-9d65-f660e189a154
md"""
### Selection Strength
"""

# ╔═╡ 0bc23086-dd82-49ca-bdc7-2c844b0f7a7c
selection_strength_bond

# ╔═╡ 57fa35f2-4705-499f-8163-806884448651
selection_strength

# ╔═╡ 654f8f76-a48b-4efa-939a-f0a201dcae18
md"""
### Symmetry Breaking
"""

# ╔═╡ 8b7f4a77-66dd-4c97-897f-747cfe4e7e91
symmetry_breaking_bond

# ╔═╡ 80cd1f63-32d9-4007-aefa-8c475b5eb7a9
symmetry_breaking

# ╔═╡ a862bbea-c3a1-46dc-8178-2ae2aa76273c
md"""
### Matrix Source
"""

# ╔═╡ 28605e7e-368a-41e6-99f9-2b592adb086e
matrix_source_bond

# ╔═╡ 68d724a6-ec93-450f-a9cb-523e7400fb21
md"""
### B-Factor
"""

# ╔═╡ a454013f-60fc-43fc-b862-86618ee220c0
B_factor_bond

# ╔═╡ d237ab6a-09d8-4010-aa7b-6e627354f116
B_factor

# ╔═╡ a0c268e0-0066-4456-a0bd-687d570ee805
begin
	df_timeseries_statistics_selected = @rsubset(df_timeseries_statistics, :selection_strength == selection_strength, :adj_matrix_source == matrix_source, :symmetry_breaking == symmetry_breaking, :factor == B_factor)

	dict_timeseries_statistics_selected = Dict()
	if nrow(df_timeseries_statistics_selected) == 1
		dict_timeseries_statistics_selected = Dict(names(df_timeseries_statistics_selected[1,:]) .=> values(df_timeseries_statistics_selected[1,:]))

		nothing
	elseif nrow(df_timeseries_statistics_selected) < 1
		println("No datasets found")

		println("Available datasets:")
		df_timeseries
	else
		println("More than one dataset matching conditions found: please fully specify")
		
		df_timeseries_statistics_selected
	end
end

# ╔═╡ 6f015d85-ee58-4a7f-b0a5-ddd846f9a6cd
!isempty(dict_timeseries_statistics_selected) ? local_moran_interaction.get_frame(dict_timeseries_selected, 560000; adj_matrix_source = matrix_source)[1] : nothing

# ╔═╡ 656e9cda-4600-4396-9183-0663a56a80a1
!isempty(dict_timeseries_statistics_selected) ? local_moran_interaction.generate_timeseries_plot(dict_timeseries_statistics_selected; time_steps = time_steps_timeseries) : nothing

# ╔═╡ 4e47651d-6e33-46a9-82b3-e1e614f19d62
md"""
# Game Types
"""

# ╔═╡ 9fb6798e-2bb6-4931-91aa-6ffacc695f56
# ╠═╡ show_logs = false
!isempty(dict_timeseries_statistics_selected) ? local_moran_interaction.generate_game_type_distribution_vs_asymmetry_plot(df_timeseries_statistics; B_factor = B_factor, selection_strength = selection_strength, adj_matrix_source = matrix_source ) : nothing

# ╔═╡ 9670bb93-9d11-488a-98c3-6e7fca2dd621
md"""
# Time evolution animation
"""

# ╔═╡ 921fcfc1-bf04-4901-bce9-0cacfaff0a1c
md"""
Note: the following cell is disabled because generating an animation is the slowest part when updating parameters (up to 10 seconds each time).

If you would like to enable the animation click the three dots (...) to the bottom-right of the image and select `Enable cell`"""

# ╔═╡ 0b2f5bdc-3799-4152-a59f-7bafb29bb1fd
# ╠═╡ disabled = true
#=╠═╡
!isempty(dict_timeseries_selected) ? generate_graph_evolution(dict_timeseries_selected, graph) : nothing
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═3e050942-505f-4328-be97-d6e72b260be2
# ╠═0ecb3c44-a2d0-11ef-1d79-d3e37535b2dc
# ╠═f10a29cf-c5ff-4e1d-8864-6efab1634b48
# ╠═fdf8f135-7dda-45ef-b8f1-08276e1e0375
# ╟─f5d72e8a-c624-411d-8b8f-a721c44d401f
# ╟─93340069-644f-4da1-b0af-824bf61a77ca
# ╟─b9040b04-c6a4-428f-9ac2-518760c46ce2
# ╟─1740d8de-7c44-4a69-8711-5c25149c560a
# ╟─26864283-0e1d-44e9-819f-44fc3e4c3aec
# ╟─7c019aea-78ee-46b5-b9ee-3489dfa2570f
# ╟─5b6bdcd4-f608-4ea9-baec-d0a62d841bc0
# ╟─adbf6ba7-420c-4d03-a92d-29f1c85929ca
# ╟─ebf8ec31-a391-4a6c-864f-efd3f732d0c2
# ╟─7e0e13aa-524f-4137-a6f6-3d4ee46a420e
# ╟─3cc72edb-1a4b-4ac6-a18f-2ad9296b2da2
# ╟─2f8d5984-acd3-4958-afe6-acc136d1c70a
# ╟─9e64dcc2-574e-4690-a10d-9b45e0ba5219
# ╟─49e07eb6-c2fc-4223-85a9-aa9c986e34da
# ╠═6f015d85-ee58-4a7f-b0a5-ddd846f9a6cd
# ╟─0ca5e4f8-faab-4057-a823-956c3711b5d8
# ╟─2dabcbb3-7dc6-4bc4-af4f-8e48407130f8
# ╟─4fefdb80-5168-44e3-ae7e-e92ef29a30e5
# ╟─abd83b14-d699-42cc-b0f0-6c1703c94ed3
# ╟─32510bf3-a932-4394-a4fb-8101580a6140
# ╟─656dd0d7-4e47-422c-a104-9b2473a51306
# ╟─e8890847-13d7-447c-9d65-f660e189a154
# ╟─0bc23086-dd82-49ca-bdc7-2c844b0f7a7c
# ╟─57fa35f2-4705-499f-8163-806884448651
# ╟─654f8f76-a48b-4efa-939a-f0a201dcae18
# ╟─8b7f4a77-66dd-4c97-897f-747cfe4e7e91
# ╟─80cd1f63-32d9-4007-aefa-8c475b5eb7a9
# ╟─a862bbea-c3a1-46dc-8178-2ae2aa76273c
# ╟─28605e7e-368a-41e6-99f9-2b592adb086e
# ╟─68d724a6-ec93-450f-a9cb-523e7400fb21
# ╟─a454013f-60fc-43fc-b862-86618ee220c0
# ╟─d237ab6a-09d8-4010-aa7b-6e627354f116
# ╟─a0c268e0-0066-4456-a0bd-687d570ee805
# ╟─656e9cda-4600-4396-9183-0663a56a80a1
# ╟─4e47651d-6e33-46a9-82b3-e1e614f19d62
# ╟─9fb6798e-2bb6-4931-91aa-6ffacc695f56
# ╟─9670bb93-9d11-488a-98c3-6e7fca2dd621
# ╟─921fcfc1-bf04-4901-bce9-0cacfaff0a1c
# ╠═0b2f5bdc-3799-4152-a59f-7bafb29bb1fd
