# SPDX-License-Identifier: GPL-3.0-or-later

adj_matrix_source_vals = ["well-mixed", "c-elegans", "c-elegans-undirected", "c-elegans-unweighted"]
symmetry_breaking_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
selection_strength_vals = [0.005, 0.2, 5.0]
B_to_c_vals = [1.5, 2.5]

beta_to_B_netcdf = 0.95
cost_netcdf = 0.1
mutation_rate_netcdf = 0.0001
nb_phases_netcdf = 20
nb_players_netcdf = 20
time_steps_cumulative_netcdf = 200000000
time_steps_timeseries_netcdf = 8000000
num_seed_vals_netcdf=2
seed_vals_netcdf = range(1,num_seed_vals_netcdf+1)

rule all:
  input:
    "papers/primary_manuscript/Report.pdf",
    "papers/primary_manuscript/SupplementaryInformation.pdf",
    # Animations for supplemental information
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.mp4",
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.mp4",
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=1.0_time_steps=8000000.mp4",
    expand(["data/processed/netcdf/cumulative_adj_matrix_source={adj_matrix_source}_time_steps={time_steps_cumulative}.nc",
      "data/processed/netcdf/timeseries_statistics_adj_matrix_source={adj_matrix_source}_time_steps={time_steps_timeseries}.nc"], adj_matrix_source=adj_matrix_source_vals,
      time_steps_cumulative=time_steps_cumulative_netcdf,time_steps_timeseries=time_steps_timeseries_netcdf),


rule manuscript:
  priority: 50 # Increase priority relative to netcdf to it generates first
  input:
    "papers/primary_manuscript/preamble.sty",
    "papers/primary_manuscript/custom-definitions.tex",
    "papers/primary_manuscript/references.bib",
    "papers/primary_manuscript/plos2025.bst",
    "papers/primary_manuscript/tikz/preamble.tex",
    "papers/primary_manuscript/tikz/c-elegans.pdf",
    "papers/primary_manuscript/tikz/chimera-states.pdf",
    "papers/primary_manuscript/tikz/game-types.pdf",
    "papers/primary_manuscript/tikz/model-setup.pdf",
    "papers/primary_manuscript/tikz/phase-diagram.pdf",
    "papers/primary_manuscript/tikz/well-mixed.pdf",
    # Info on number of edges and self loops
    "data/processed/graph_loop_edge_number/adj_matrix_source=c-elegans.csv",
    # Info on mutation times
    "data/processed/mutation_timesteps/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
    "data/processed/mutation_timesteps/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
    tex="papers/primary_manuscript/Report.tex",
  output:
    "papers/primary_manuscript/Report.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_c_elegans:
  input:
    "data/processed/graph_structure_snapshot/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.75_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=c-elegans.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans-unweighted_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans-undirected_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.csv",
    tex="papers/primary_manuscript/tikz/c-elegans.tex",
  output:
    "papers/primary_manuscript/tikz/c-elegans.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_chimera_states:
  input:
    "data/processed/chimeraindex/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_community_algorithm=leiden_community_beta=0.01_community_n_iter=2_community_resolution=0.1_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/graph_structure_snapshot/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=0.0_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure_snapshot/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_seed=12345_selection_strength=0.2_symmetry_breaking=1.0_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=c-elegans.csv",
    tex="papers/primary_manuscript/tikz/chimera-states.tex",
  output:
    "papers/primary_manuscript/tikz/chimera-states.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_game_payoffs:
  input:
    tex="papers/primary_manuscript/tikz/game-payoffs.tex",
  output:
    "papers/primary_manuscript/tikz/game-payoffs.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_game_types:
  input:
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_only_mixed_games=false_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_only_mixed_games=false_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_time_steps=8000000.csv",
    tex="papers/primary_manuscript/tikz/game-types.tex",
  output:
    "papers/primary_manuscript/tikz/game-types.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_model_setup:
  input:
    "data/processed/graph_structure/vertices_adj_matrix_source=well-mixed.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=well-mixed.csv",
    tex="papers/primary_manuscript/tikz/model-setup.tex",
  output:
    "papers/primary_manuscript/tikz/model-setup.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_phase_diagram:
  input:
    tex="papers/primary_manuscript/tikz/phase-diagram.tex",
  output:
    "papers/primary_manuscript/tikz/phase-diagram.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule tikz_well_mixed:
  input:
    "data/processed/cumulative_theory/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_nb_players=20_selection_strength=0.2_symmetry_breaking=0.0_type=theory.csv",
    "data/processed/cumulative_theory/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_nb_players=20_selection_strength=0.2_symmetry_breaking=0.25_type=theory.csv",
    "data/processed/cumulative_theory/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_nb_players=20_selection_strength=0.2_symmetry_breaking=0.5_type=theory.csv",
    "data/processed/cumulative_theory/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_nb_players=20_selection_strength=0.2_symmetry_breaking=0.75_type=theory.csv",
    "data/processed/cumulative_theory/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_nb_players=20_selection_strength=0.2_symmetry_breaking=1.0_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.0_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.25_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.5_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_selection_strength=0.2_symmetry_breaking=1.0_time_steps=200000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_num_seeds=10_only_mixed_games=true_selection_strength=0.2_symmetry_breaking=1.0_time_steps=8000000.csv",
    tex="papers/primary_manuscript/tikz/well-mixed.tex",
  output:
    "papers/primary_manuscript/tikz/well-mixed.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

rule supplementary_information:
  priority: 50 # Increase priority relative to netcdf to it generates first
  input:
    "papers/primary_manuscript/preamble.sty",
    "papers/primary_manuscript/tikz/preamble.tex",
    "papers/primary_manuscript/tikz/game-payoffs.tex",
    tex="papers/primary_manuscript/SupplementaryInformation.tex",
  output:
    "papers/primary_manuscript/SupplementaryInformation.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -interaction=nonstopmode ../../{input.tex}"

wildcard_constraints:
    adj_matrix_source=r"[a-z\-]+",
    beta_to_b=r"[\d.]+",
    community_algorithm=r"[a-z\-]+",
    cost=r"[\d.]+",
    data_type=r"[a-z\-]+",
    early_cutoff_fraction_flag=r"_early_cutoff_fraction=|",
    early_cutoff_fraction=r"[\d.]*", # use * instead of + because this is an optional parameter
    only_mixed_games=r"[a-z]+",
    fixed_aspect_flag=r"_fixed_aspect=|",
    fixed_aspect=r"[a-z]*", # use * instead of + because this is an optional parameter
    covariance_cutoff_fraction_flag=r"_covariance_cutoff_fraction=|",
    covariance_cutoff_fraction=r"[\d.]*", # use * instead of + because some algorithms don't use this
    walktrap_steps_flag=r"_walktrap_steps=|",
    walktrap_steps=r"\d*", # use * instead of + because some algorithms don't use this
    community_resolution_flag=r"_community_resolution=|",
    community_resolution=r"[\d.]*", # use * instead of + because some algorithms don't use this
    community_beta_flag=r"_community_beta=|",
    community_beta=r"[\d.]*", # use * instead of + because some algorithms don't use this
    community_n_iter_flag=r"_community_n_iter=|",
    community_n_iter=r"\d*", # use * instead of + because some algorithms don't use this
    mutation_rate=r"[\d.]+",
    nb_phases=r"\d+",
    nb_players_flag=r"_nb_players=|",
    nb_players=r"\d*", # use * instead of + because this will be empty for non "well-mixed" datasets
    num_seeds_flag=r"_num_seeds=|",
    num_seeds=r"\d*", # use * instead of + because this will be empty for theory cumulative processing
    seed=r"\d+",
    selection_strength=r"[\d.]+",
    symmetry_breaking=r"[\d.]+",
    time_steps=r"\d*",

def external_data_filename(wildcards):
  if wildcards.adj_matrix_source in ["c-elegans",
    "c-elegans-undirected", "c-elegans-unweighted",
    "c-elegans-undirected-unweighted"]:
    return "data/external/store/k12:33f530402e3f4ac495d8fafe8515269f.xlsx"
  return []

safety_factor_memory = 2.5
safety_factor_runtime = 2.0

rule external_celegans:
  resources:
    mem_mb = 3000,
    runtime = 10,
  params:
    adj_matrix_source="c-elegans",
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    "scripts/snakemake/cache_external_data.jl",
  output:
    protected("data/external/store/k12:33f530402e3f4ac495d8fafe8515269f.xlsx"),
  script:
    "scripts/snakemake/cache_external_data.jl"

def num_players(wildcards):
  if wildcards.adj_matrix_source in ["c-elegans",
    "c-elegans-undirected", "c-elegans-unweighted",
    "c-elegans-undirected-unweighted"]:
    return 300
  elif wildcards.adj_matrix_source  == "well-mixed":
    return int(wildcards.nb_players)
  else:
    raise("Unknown `adj_matrix_source`")

def memory_mb_timeseries(wildcards):
  nb_players = num_players(wildcards)
  time_steps = int(wildcards.time_steps)
  float_size_mb = 8/1E6
  # mem_mb = 1300, # well-mixed
  # mem_mb = 19500, # c-elegans
  memory_mb = (time_steps + nb_players)*float_size_mb*safety_factor_memory
  return memory_mb

def runtime_min_cumulative(wildcards):
  nb_players = num_players(wildcards)
  time_steps = int(wildcards.time_steps)
  # runtime = 23 # well-mixed time_steps=2E8
  # runtime = 533 # c-elegans time_steps=2E8
  runtime_min_per_player_per_timestep = 1.8/2E8
  runtime_min = runtime_min_per_player_per_timestep*nb_players*time_steps*safety_factor_runtime
  return runtime_min

rule raw_cumulative_data:
  resources:
    mem_mb = 2500, # well-mixed and c-elegans
    runtime = runtime_min_cumulative,
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/moran.jl",
    "src/julia/run_simulations.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/run_cumulative_simulations.jl",
  output:
    "data/raw/cumulative/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  script:
    "scripts/snakemake/run_cumulative_simulations.jl"

def raw_cumulative_filename(wildcards):
	return expand(["data/raw/cumulative/adj_matrix_source={{adj_matrix_source}}_beta_to_B={{beta_to_B}}_cost={{cost}}{{fixed_aspect_flag}}{{fixed_aspect}}_mutation_rate={{mutation_rate}}_nb_phases={{nb_phases}}{{nb_players_flag}}{{nb_players}}_seed={seed}_selection_strength={{selection_strength}}_symmetry_breaking={{symmetry_breaking}}_time_steps={{time_steps}}.jld2"], seed = range(1,int(wildcards.num_seeds)+1))

rule processed_cumulative_data:
  resources:
    mem_mb = 2500, # well-mixed and c-elegans
    runtime = 10, # well-mixed and c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_cumulative.jl",
    raw_cumulative_filename,
  output:
    "data/processed/cumulative/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}{num_seeds_flag}{num_seeds}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_cumulative.jl"

rule cumulative_theory:
  resources:
    mem_mb = 4000, # well-mixed and c-elegans
    runtime = 10, # well-mixed and c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/postprocess.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_cumulative_theory.jl",
  output:
    "data/processed/cumulative_theory/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_type={type}.csv",
  script:
    "scripts/snakemake/postprocess_cumulative_theory.jl"

rule raw_timeseries_data:
  resources:
    mem_mb = memory_mb_timeseries,
    # runtime = 2, # well-mixed
    runtime = 10, # c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/moran.jl",
    "src/julia/run_simulations.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/run_timeseries_simulations.jl",
  output:
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  script:
    "scripts/snakemake/run_timeseries_simulations.jl"

def runtime_min_timeseries_statistics(wildcards):
  nb_players = num_players(wildcards)
  time_steps = int(wildcards.time_steps)
  # runtime = 8 # well-mixed time_steps=8E6
  # runtime = 110 # c-elegans time_steps=8E6
  runtime_min_per_player_per_timestep = 0.4/8E6
  runtime_min = runtime_min_per_player_per_timestep*nb_players*time_steps*safety_factor_runtime
  return runtime_min

def raw_timeseries_filename(wildcards):
  return expand(["data/raw/timeseries/B_to_c={{B_to_c}}_adj_matrix_source={{adj_matrix_source}}_beta_to_B={{beta_to_B}}_cost={{cost}}{{fixed_aspect_flag}}{{fixed_aspect}}_mutation_rate={{mutation_rate}}_nb_phases={{nb_phases}}{{nb_players_flag}}{{nb_players}}_seed={seed}_selection_strength={{selection_strength}}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"], seed = range(1,int(wildcards.num_seeds)+1),
	symmetry_breaking = wildcards.symmetry_breaking if hasattr(wildcards, "symmetry_breaking") else symmetry_breaking_vals)

rule processed_timeseries_statistics:
  resources:
    mem_mb = memory_mb_timeseries,
    runtime = runtime_min_timeseries_statistics,
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_timeseries_statistics.jl",
    raw_timeseries_filename,
  output:
    "data/processed/timeseries_statistics/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{early_cutoff_fraction_flag}{early_cutoff_fraction}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_num_seeds={num_seeds}_only_mixed_games={only_mixed_games}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_timeseries_statistics.jl"

rule processed_mutation_timesteps:
  resources:
    # mem_mb = 155, # well-mixed
    mem_mb = 10000, # c-elegans
    runtime = 10, # well-mixed and c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_mutation_timesteps.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/mutation_timesteps/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_mutation_timesteps.jl"

rule processed_animations:
  resources:
    mem_mb = memory_mb_timeseries,
    runtime = 30, # Guess
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/plot_animations.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "plots/animations/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.mp4",
  script:
    "scripts/snakemake/plot_animations.jl"

rule processed_graph_structure:
  resources:
    mem_mb = 4000, # well-mixed and c-elegans
    runtime = 10, # well-mixed and c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/postprocess.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_graph_structure.jl",
  output:
    "data/processed/graph_structure/edges_adj_matrix_source={adj_matrix_source}.csv",
    "data/processed/graph_structure/vertices_adj_matrix_source={adj_matrix_source}.csv",
  script:
    "scripts/snakemake/postprocess_graph_structure.jl"

rule processed_graph_structure_snapshot:
  resources:
    mem_mb = 6000, # c-elegans
    runtime = 10, # c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_graph_structure_snapshot.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/graph_structure_snapshot/vertices_B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_step={time_step}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_graph_structure_snapshot.jl"

rule processed_graph_loop_edge_number:
  resources:
    mem_mb = 4000, # c-elgans
    runtime = 10, # c-elgans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_graph_loop_edge_number.jl",
  output:
    "data/processed/graph_loop_edge_number/adj_matrix_source={adj_matrix_source}.csv",
  script:
    "scripts/snakemake/postprocess_graph_loop_edge_number.jl"

def runtime_min_game_types(wildcards):
  nb_players = num_players(wildcards)
  time_steps = int(wildcards.time_steps)
  # runtime = 30 # well-mixed time_steps=8E6
  # runtime = 540 # c-elegans time_steps=8E6
  runtime_min_per_player_per_timestep = 1.8/8E6
  runtime_min = runtime_min_per_player_per_timestep*nb_players*time_steps*safety_factor_runtime
  return runtime_min

def memory_mb_game_types(wildcards):
  nb_players = num_players(wildcards)
  time_steps = int(wildcards.time_steps)
  float_size_mb = 8/1E6
  # mem_mb = 2300, # well-mixed
  # mem_mb = 17000, # c-elegans
  fudge_factor = 2
  memory_mb = fudge_factor*float_size_mb*nb_players*time_steps*safety_factor_memory
  return memory_mb

rule processed_game_types:
  resources:
    mem_mb = memory_mb_game_types,
    runtime = runtime_min_game_types,
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_game_types.jl",
    raw_timeseries_filename,
  output:
    "data/processed/gametype/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_num_seeds={num_seeds}_only_mixed_games={only_mixed_games}_selection_strength={selection_strength}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_game_types.jl"

rule processed_chimera_indices:
  resources:
    mem_mb = 25000, # c-elegans
    runtime = 20, # c-elegans
  input:
    "Project.toml",
    "Manifest.toml",
    "Data.toml",
    external_data_filename,
    "scripts/snakemake/cache_external_data.jl",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_chimera_indices.jl",
    raw_timeseries_filename,
  output:
    "data/processed/chimeraindex/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_community_algorithm={community_algorithm}{community_beta_flag}{community_beta}{community_n_iter_flag}{community_n_iter}{community_resolution_flag}{community_resolution}_cost={cost}{covariance_cutoff_fraction_flag}{covariance_cutoff_fraction}{fixed_aspect_flag}{fixed_aspect}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_num_seeds={num_seeds}_selection_strength={selection_strength}_time_steps={time_steps}{walktrap_steps_flag}{walktrap_steps}.csv",
  script:
    "scripts/snakemake/postprocess_chimera_indices.jl"

def netcdf_timeseries_inputs(wildcards):
  return expand(["data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={{adj_matrix_source}}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"],
    selection_strength=selection_strength_vals, symmetry_breaking=symmetry_breaking_vals,
    seed=seed_vals_netcdf,
    B_to_c=B_to_c_vals,
    beta_to_B=beta_to_B_netcdf, cost=cost_netcdf,
    mutation_rate=mutation_rate_netcdf, nb_phases=nb_phases_netcdf,
    time_steps=time_steps_cumulative_netcdf,
    nb_players_flag=("_nb_players=" if wildcards.adj_matrix_source=="well-mixed" else ""),
    nb_players=(nb_players_netcdf if wildcards.adj_matrix_source=="well-mixed" else ""),
    )

def netcdf_cumulative_inputs(wildcards):
  return expand(["data/raw/cumulative/adj_matrix_source={{adj_matrix_source}}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_seed={seed}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"],
    selection_strength=selection_strength_vals, symmetry_breaking=symmetry_breaking_vals,
    seed=seed_vals_netcdf,
    beta_to_B=beta_to_B_netcdf, cost=cost_netcdf,
    mutation_rate=mutation_rate_netcdf, nb_phases=nb_phases_netcdf,
    time_steps=time_steps_cumulative_netcdf,
    nb_players_flag=("_nb_players=" if wildcards.adj_matrix_source=="well-mixed" else ""),
    nb_players=(nb_players_netcdf if wildcards.adj_matrix_source=="well-mixed" else ""),
    )

rule netcdf_timeseries:
  resources:
    mem_mb = 2500, # cumulative well-mixed and c-elegans
    runtime = 3, # cumulative well-mixed and c-elegans
  input:
    netcdf_timeseries_inputs,
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_netcdf_timeseries_statistics.jl",
  output:
    "data/processed/netcdf/timeseries_statistics_adj_matrix_source={adj_matrix_source}_time_steps={time_steps}.nc",
  script:
    "scripts/snakemake/postprocess_netcdf_timeseries_statistics.jl"

rule netcdf_cumulative:
  resources:
    mem_mb = 2500, # cumulative well-mixed and c-elegans
    runtime = 3, # cumulative well-mixed and c-elegans
  input:
    netcdf_cumulative_inputs,
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_netcdf_cumulative.jl",
  output:
    "data/processed/netcdf/cumulative_adj_matrix_source={adj_matrix_source}_time_steps={time_steps}.nc",
  script:
    "scripts/snakemake/postprocess_netcdf_cumulative.jl"
