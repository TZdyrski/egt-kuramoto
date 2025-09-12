adj_matrix_source_vals = ["well-mixed", "c-elegans", "c-elegans-undirected", "c-elegans-unweighted"]
symmetry_breaking_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
selection_strength_vals = [0.005, 0.2, 5.0]

B_to_c_netcdf = 1.5,
beta_to_B_netcdf = 0.95,
cost_netcdf = 0.1,
mutation_rate_netcdf = 0.0001,
nb_phases_netcdf = 20,
nb_players_netcdf = 20,
time_steps_cumulative_netcdf = 200000000,
time_steps_timeseries_netcdf = 8000000,

rule all:
  input:
    "papers/primary_manuscript/Report.pdf",
    # Animations for supplemental information
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.mp4",
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.mp4",
    "plots/animations/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=1.0_time_steps=8000000.mp4",
    expand(["data/processed/netcdf/cumulative_matrixSource={adj_matrix_source}_timesteps={time_steps_cumulative}.nc",
      "data/processed/netcdf/timeseries-statistics_decimationFactor=1000_matrixSource={adj_matrix_source}_timesteps={time_steps_timeseries}.nc"], adj_matrix_source=adj_matrix_source_vals,
      time_steps_cumulative=time_steps_cumulative_netcdf,time_steps_timeseries=time_steps_timeseries_netcdf,),

rule manuscript:
  priority: 50 # Increase priority relative to netcdf to it generates first
  input:
    # Manuscript
    "papers/primary_manuscript/Report.tex",
    "papers/primary_manuscript/preamble.sty",
    "papers/primary_manuscript/custom-definitions.tex",
    "papers/primary_manuscript/references.bib",
    "papers/primary_manuscript/sn-jnl.cls",
    "papers/primary_manuscript/sn-nature.bst",
    "papers/primary_manuscript/sections/appendix.tex",
    "papers/primary_manuscript/tikz/preamble.tex",
    # tikz/c-elegans.tex
    "papers/primary_manuscript/tikz/c-elegans.tex",
    "data/processed/graph_structure/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=c-elegans.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans-unweighted_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=c-elegans-undirected_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000_type=simulation.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.csv",
    # tikz/chimera-states.tex
    "papers/primary_manuscript/tikz/chimera-states.tex",
    "data/processed/chimeraindex/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_community_algorithm=covariance_cost=0.1_covariance_cutoff=1500_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/graph_structure/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.0_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure/vertices_B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=1.0_time_step=560000_time_steps=8000000.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=c-elegans.csv",
    # tikz/game-payoffs.tex
    "papers/primary_manuscript/tikz/game-payoffs.tex",
    # tikz/game-types.tex
    "papers/primary_manuscript/tikz/game-types.tex",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_exclude_CC_NN=true_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_time_steps=8000000.csv",
    "data/processed/gametype/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_exclude_CC_NN=true_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_time_steps=8000000.csv",
    # tikz/model-setup.tex
    "papers/primary_manuscript/tikz/model-setup.tex",
    "data/processed/graph_structure/vertices_adj_matrix_source=well-mixed.csv",
    "data/processed/graph_structure/edges_adj_matrix_source=well-mixed.csv",
    # tikz/phase-diagram.tex
    "papers/primary_manuscript/tikz/phase-diagram.tex",
    # tikz/well-mixed.tex
    "papers/primary_manuscript/tikz/well-mixed.tex",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=200000000_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.25_time_steps=200000000_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.5_time_steps=200000000_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=1.0_time_steps=200000000_type=theory.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.25_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.5_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=200000000_type=simulation.csv",
    "data/processed/cumulative/adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=1.0_time_steps=200000000_type=simulation.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.75_time_steps=8000000.csv",
    "data/processed/timeseries_statistics/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_early_cutoff_fraction=0.1_exclude_CC_NN=false_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=1.0_time_steps=8000000.csv",
    # Info on number of edges and self loops
    "data/processed/graph_loop_edge_number/adj_matrix_source=c-elegans.csv",
    # Info on mutation times
    "data/processed/mutation_timesteps/B_to_c=1.5_adj_matrix_source=well-mixed_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_nb_players=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
    "data/processed/mutation_timesteps/B_to_c=1.5_adj_matrix_source=c-elegans_beta_to_B=0.95_cost=0.1_mutation_rate=0.0001_nb_phases=20_selection_strength=0.2_symmetry_breaking=0.0_time_steps=8000000.csv",
  output:
    "papers/primary_manuscript/Report.pdf",
  shell:
    "cd papers/primary_manuscript && latexmk -lualatex -auxdir=latex.aux -interaction=nonstopmode Report.tex"

rule netcdf_datasets:
  input:
    expand(["data/processed/netcdf/cumulative_matrixSource={adj_matrix_source}_timesteps={time_steps_cumulative}.nc",
      "data/processed/netcdf/timeseries-statistics_decimationFactor=1000_matrixSource={adj_matrix_source}_timesteps={time_steps_timeseries}.nc"], adj_matrix_source=adj_matrix_source_vals,
      time_steps_cumulative=time_steps_cumulative_netcdf,time_steps_timeseries=time_steps_timeseries_netcdf,),

wildcard_constraints:
    adj_matrix_source=r"[a-z\-]+",
    beta_to_b=r"[\d.]+",
    community_algorithm=r"[a-z\-]+",
    cost=r"[\d.]+",
    data_type=r"[a-z\-]+",
    early_cutoff_fraction_flag=r"_early_cutoff_fraction=|",
    early_cutoff_fraction=r"[\d.]*", # use * instead of + because this is an optional parameter
    exclude_CC_NN=r"[a-z]+",
    covariance_cutoff_flag=r"_covariance_cutoff=|",
    covariance_cutoff=r"\d*", # use * instead of + because some algorithms don't use covariance
    mutation_rate=r"[\d.]+",
    nb_phases=r"\d+",
    nb_players_flag=r"_nb_players=|",
    nb_players=r"\d*", # use * instead of + because this will be empty for non "well-mixed" datasets
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
  memory_mb = float_size_mb*nb_players*time_steps*safety_factor_memory
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
    "data/raw/cumulative/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  script:
    "scripts/snakemake/run_cumulative_simulations.jl"

rule processed_cumulative_data:
  resources:
    mem_mb = 2500, # well-mixed and c-elegans
    runtime = 10, # well-mixed and c-elegans
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
    "scripts/snakemake/postprocess_cumulative.jl",
    "data/raw/cumulative/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/cumulative/adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}_type={type}.csv",
  script:
    "scripts/snakemake/postprocess_cumulative.jl"

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
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
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

rule processed_timeseries_statistics:
  resources:
    mem_mb = memory_mb_timeseries,
    runtime = runtime_min_timeseries_statistics,
  params:
    num_samples=1000,
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_timeseries_statistics.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/timeseries_statistics/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}{early_cutoff_fraction_flag}{early_cutoff_fraction}_exclude_CC_NN={exclude_CC_NN}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_timeseries_statistics.jl"

rule processed_mutation_timesteps:
  resources:
    # mem_mb = 1500, # well-mixed
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
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/mutation_timesteps/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_mutation_timesteps.jl"

rule processed_animations:
  input:
    "Project.toml",
    "Manifest.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/plot_animations.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "plots/animations/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.mp4",
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
    "src/julia/game_taxonomy.jl",
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
    "scripts/snakemake/postprocess_graph_structure.jl",
    "data/raw/timeseries/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={time_steps}.jld2",
  output:
    "data/processed/graph_structure/vertices_B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_step={time_step}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_graph_structure.jl"

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
    expand(["data/raw/timeseries/B_to_c={{B_to_c}}_adj_matrix_source={{adj_matrix_source}}_beta_to_B={{beta_to_B}}_cost={{cost}}_mutation_rate={{mutation_rate}}_nb_phases={{nb_phases}}{{nb_players_flag}}{{nb_players}}_selection_strength={{selection_strength}}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"], symmetry_breaking=symmetry_breaking_vals),
  output:
    "data/processed/gametype/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_cost={cost}_exclude_CC_NN={exclude_CC_NN}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_time_steps={time_steps}.csv",
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
    "CondaPkg.toml",
    "src/julia/postprocess.jl",
    "src/julia/game_taxonomy.jl",
    "src/julia/utils.jl",
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_chimera_indices.jl",
    expand(["data/raw/timeseries/B_to_c={{B_to_c}}_adj_matrix_source={{adj_matrix_source}}_beta_to_B={{beta_to_B}}_cost={{cost}}_mutation_rate={{mutation_rate}}_nb_phases={{nb_phases}}{{nb_players_flag}}{{nb_players}}_selection_strength={{selection_strength}}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"], symmetry_breaking=symmetry_breaking_vals),
  output:
    "data/processed/chimeraindex/B_to_c={B_to_c}_adj_matrix_source={adj_matrix_source}_beta_to_B={beta_to_B}_community_algorithm={community_algorithm}_cost={cost}{covariance_cutoff_flag}{covariance_cutoff}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_time_steps={time_steps}.csv",
  script:
    "scripts/snakemake/postprocess_chimera_indices.jl"

def inputs_dependent_on_adj_matrix_source(wildcards):
  return expand(["data/raw/{data_type_input}/{B_to_c_flag}{B_to_c}{B_to_c_postunderscore}adj_matrix_source={{adj_matrix_source}}_beta_to_B={beta_to_B}_cost={cost}_mutation_rate={mutation_rate}_nb_phases={nb_phases}{nb_players_flag}{nb_players}_selection_strength={selection_strength}_symmetry_breaking={symmetry_breaking}_time_steps={{time_steps}}.jld2"],
    selection_strength=selection_strength_vals, symmetry_breaking=symmetry_breaking_vals,
    B_to_c_flag=("B_to_c=" if wildcards.data_type!="cumulative" else ""),
    B_to_c=(B_to_c_netcdf if wildcards.data_type!="cumulative" else ""),
    B_to_c_postunderscore=("_" if wildcards.data_type!="cumulative" else ""),
    beta_to_B=beta_to_B_netcdf, cost=cost_netcdf,
    mutation_rate=mutation_rate_netcdf, nb_phases=nb_phases_netcdf,
    time_steps=time_steps_cumulative_netcdf,
    nb_players_flag=("_nb_players=" if wildcards.adj_matrix_source=="well-mixed" else ""),
    nb_players=(nb_players_netcdf if wildcards.adj_matrix_source=="well-mixed" else ""),
    # We want to convert timeseries-statistics to timeseries, so split before "-"
    data_type_input=wildcards.data_type.split("-")[0],
    )

rule netcdf_dataset:
  resources:
    mem_mb = 20000, # Guess
    runtime = 300, # Guess
  input:
    inputs_dependent_on_adj_matrix_source,
    "scripts/snakemake/snakemake_preamble.jl",
    "scripts/snakemake/postprocess_netcdf.jl",
  output:
    "data/processed/netcdf/{data_type}_matrixSource={adj_matrix_source}_timesteps={time_steps}.nc",
  script:
    "scripts/snakemake/postprocess_netcdf.jl"

use rule netcdf_dataset as netcdf_dataset_decimated with:
  output:
    "data/processed/netcdf/{data_type}_decimationFactor={decimation_factor}_matrixSource={adj_matrix_source}_timesteps={time_steps}.nc",
