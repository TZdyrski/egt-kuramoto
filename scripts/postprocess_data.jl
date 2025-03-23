using DrWatson
@quickactivate "2024_EGT_Kuramoto"

# Here you may include files from the source directory
include(srcdir("julia", "local_moran_interaction.jl"))

# Save specific graph timestep
time_step = 560000
local_moran_interaction.export_graph_nodes_edges(
    time_step;B_factor=1.5, selection_strength=0.2,
    symmetry_breaking=0.75, adj_matrix_source="c-elegans",
    time_steps=8000000)
