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
	using PGFPlotsX
end

# ╔═╡ fdf8f135-7dda-45ef-b8f1-08276e1e0375
@plutoinclude srcdir("julia", "local_moran_interaction.jl")

# ╔═╡ 1d7658e1-0484-448c-a03b-117ca5ac7e84
begin
	function show_tikz(pic::PGFPlotsX.TikzPicture)
		game_colors_sorted = [x.second for x in sort(collect(local_moran_interaction.game_type_colors), by=x->(Int(x[1])))]

		game_colors_and_parity_sorted = [game_colors_sorted...,
			colorant"lightgrey", colorant"darkgrey"]

		push!(pic.options, "scale=2")
		push!(pic.options, "background rectangle/.style={fill=white}")
		push!(pic.options, "show background rectangle")
		return @pgf PGFPlotsX.TikzDocument(
		preamble = [raw"\usetikzlibrary{backgrounds}",
			("game_colors", game_colors_sorted),
			("game_colors_and_parity", game_colors_and_parity_sorted),
			],
			pic
		)
	end

	function show_tikz(ax::PGFPlotsX.AxisLike)
		return @pgf show_tikz(PGFPlotsX.TikzPicture(ax))
	end
end;

# ╔═╡ f5d72e8a-c624-411d-8b8f-a721c44d401f
md"""
# Communicative Fraction vs. $B(0)$
"""

# ╔═╡ 93340069-644f-4da1-b0af-824bf61a77ca
# ╠═╡ show_logs = false
begin
	time_steps_cumulative = Int(2E8)
	
	df_raw_cumulative = collect_results(datadir("cumulative"); rinclude = [Regex("time_steps=$time_steps_cumulative[._]")])
	
	df_cumulative = transform(df_raw_cumulative, :path => (x-> DataFrame(map(y -> parse_savename(y)[2], x))) => AsTable)
end;

# ╔═╡ b9040b04-c6a4-428f-9ac2-518760c46ce2
begin
	selection_strength_bond = @bind selection_strength PlutoUI.Slider([0.005,0.2,5]);
	symmetry_breaking_bond = @bind symmetry_breaking PlutoUI.Slider(0:0.25:1);
	matrix_source_bond = @bind matrix_source Select(["well-mixed", "c-elegans-undirected-unweighted", "c-elegans-undirected", "c-elegans-unweighted", "c-elegans", "random-regular-graph", "random-regular-digraph"]);
	B_factor_bond = @bind B_factor PlutoUI.Slider([1.5,2.5]);
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
!isempty(dict_cumulative_selected) ? show_tikz(local_moran_interaction.generate_cumulative_plot(dict_cumulative_selected, @strdict(selection_strength,symmetry_breaking,nb_phases=20,adj_matrix_source=matrix_source))) : nothing

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
	df_timeseries = local_moran_interaction.load_all_timeseries(time_steps_timeseries)
end;

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
	df_timeseries_selected = @rsubset(df_timeseries, :selection_strength == selection_strength, :adj_matrix_source == matrix_source, :symmetry_breaking == symmetry_breaking, :factor == B_factor)

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

# ╔═╡ 656e9cda-4600-4396-9183-0663a56a80a1
!isempty(dict_timeseries_selected) ? show_tikz(local_moran_interaction.generate_timeseries_plot(dict_timeseries_selected; time_steps = time_steps_timeseries)) : nothing

# ╔═╡ 4e47651d-6e33-46a9-82b3-e1e614f19d62
md"""
# Game Types
"""

# ╔═╡ 9fb6798e-2bb6-4931-91aa-6ffacc695f56
# ╠═╡ show_logs = false
plot_game_type_distribution_vs_asymmetry(B_factor,selection_strength,matrix_source,time_steps_timeseries)

# ╔═╡ Cell order:
# ╠═3e050942-505f-4328-be97-d6e72b260be2
# ╠═0ecb3c44-a2d0-11ef-1d79-d3e37535b2dc
# ╠═f10a29cf-c5ff-4e1d-8864-6efab1634b48
# ╠═fdf8f135-7dda-45ef-b8f1-08276e1e0375
# ╠═1d7658e1-0484-448c-a03b-117ca5ac7e84
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
# ╟─0ca5e4f8-faab-4057-a823-956c3711b5d8
# ╟─2dabcbb3-7dc6-4bc4-af4f-8e48407130f8
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
