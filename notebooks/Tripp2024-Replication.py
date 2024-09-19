# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: myenv
#     language: python
#     name: myenv
# ---

# %% [markdown]
# # Replicating [Tripp (2024)](https://doi.org/10.1098/rspb.2022.0999)

# %%
import numpy as np
import egttools as egt

# %% [markdown]
# ## Setup

# %% [markdown]
# ### Define Game

# %%
from typing import Callable

# Generate payoff matrix given parameters
def payoff_matrix_template(mutual_benefit: Callable,
                  unilateral_benefit: Callable,
                  cost: float, nb_phases: int) -> np.array:
    # Generate phases
    phi = 2*np.pi*np.arange(nb_phases)/nb_phases

    # Create array of phi_i and phi_j for phases of two players
    phi_i, phi_j = np.meshgrid(phi, phi)

    # Calculate phase difference
    delta_phi = phi_j - phi_i

    # Generate payoff matrix of the form
    #                 Cooperate         Defect
    #              ---------------------------------
    #   Cooperate  | mutual-cost | unilateral-cost |
    #              +-------------+-----------------+
    #   Defect     | unilateral  |        0        |
    #              ---------------------------------
    payoff_matrix = np.block([[mutual_benefit(delta_phi)-cost,
        unilateral_benefit(delta_phi)-cost],
        [unilateral_benefit(delta_phi), np.zeros([nb_phases, nb_phases])]])

    # Ensure payoff matrix is non-negative
    payoff_matrix = payoff_matrix - np.min(payoff_matrix)
    
    return payoff_matrix


# %%
import sympy as sp

# Plot game types
B_on_c, beta_on_c = sp.symbols(r'B(0)/c \beta(\phi)/c')
payoff_matrix = np.array([[B_on_c-1, beta_on_c-1],[beta_on_c, 0]])

R = payoff_matrix[0,0] # Reward
S = payoff_matrix[0,1] # Sucker's payoff
T = payoff_matrix[1,0] # Temptation
P = payoff_matrix[1,1] # Punishment

pd = sp.plot_implicit(sp.And(T >= R, R >= P, P >= S),
              x_var = (B_on_c, 0, 4), y_var = (beta_on_c, 0, 4),
              line_color = 'red', label = "Prisoner's Dilemma",
              axis_center = (0,0), show = False)

sd = sp.plot_implicit(sp.And(T >= R, R >= S, S >= P),
              line_color = 'yellow', label = "Snowdrift",
              show = False)
pd.append(sd[0])
co = sp.plot_implicit(sp.And(R >= T, T >= S, S >= P),
              line_color = 'green', label = "Coordination",
              show = False)
pd.append(co[0])
mu = sp.plot_implicit(sp.And(R >= T, T >= P, P >= S),
              line_color = 'blue', label = "Mutualism",
              show = False)
pd.append(mu[0])
unit_line = sp.plot_implicit(sp.Eq(B_on_c, beta_on_c),
              line_color='grey', show=False)
pd.append(unit_line[0])
hline = sp.plot_implicit(sp.Eq(beta_on_c, 1),
              x_var = (B_on_c), y_var = (beta_on_c),
              line_color='grey', show = False)
pd.append(hline[0])
vline = sp.plot_implicit(sp.Eq(B_on_c, 1),
              x_var = (B_on_c), y_var = (beta_on_c),
              line_color='grey', show = False)
pd.append(vline[0])
pd.show()

# %%
# Define parameters
mutual_benefit_synchronous = 1
unilateral_benefit_synchronous = 0.95*mutual_benefit_synchronous
cost = 0.1
nb_phases = 20

selection_strength = 0.005
mutation_rate = 0.0001
time_steps = 2


# %%
# Calculate payoff matrix
def create_payoff_matrix(mutual_benefit_synchronous: float, unilateral_benefit_synchronous: float, **kwargs) -> np.array:
    mutual_benefit = lambda delta_phi: mutual_benefit_synchronous*(1 + np.cos(delta_phi))/2
    unilateral_benefit = lambda delta_phi: unilateral_benefit_synchronous*(1 + np.cos(delta_phi))/2
    payoff_matrix = payoff_matrix_template(mutual_benefit, unilateral_benefit, cost, nb_phases)
    return payoff_matrix


import matplotlib.pylab as plt
# Plot payoff matrix
plt.imshow(create_payoff_matrix(mutual_benefit_synchronous=mutual_benefit_synchronous, unilateral_benefit_synchronous=unilateral_benefit_synchronous))
plt.colorbar()
plt.title("Payoff matrix")
plt.xlabel(r"$C[\phi_i], N[\phi_i]$")
plt.ylabel(r"$C[\phi_j], N[\phi_j]$");


# %%
# Define fitness function
def fitness(selection_strength: float, scores: np.array) -> np.array:
    fitness = np.exp(selection_strength*scores)
    return fitness


# %% [markdown]
# ### Define Population

# %%
# Define parameters
nb_players = 20

# %%
from collections import Counter

# Define initial population
nb_strategies = 2*nb_phases
np.random.seed(0)
initial_strategies = np.random.randint(nb_strategies, size=nb_players)
# Convert array of strategies indexed by player to array of (number of) players indexed by strategy
initial_population = np.zeros(nb_strategies)
for key, value in Counter(initial_strategies).items():
    initial_population[key] = value

# %% [markdown]
# ### Define Dynamics

# %%
# Define parameters
selection_strength_weak = 0.005
selection_strength_strong = 0.2
mutation_rate = 0.0001

# %%
from functools import partial

# Create Moran process
def create_moran_game(nb_players, nb_strategies, selection_strength, cache_size = 10, **kwargs) -> iter:
    payoff_matrix = create_payoff_matrix(**kwargs)
    game = egt.games.Matrix2PlayerGameHolder(nb_strategies = nb_strategies,
                                             payoff_matrix = payoff_matrix)
    moran = egt.numerical.PairwiseComparisonNumerical(
        pop_size = nb_players,
        game = game,
        cache_size = cache_size,
    )
    return moran


# %% [markdown]
# ## Analysis Tools

# %%
def extract_num_communicative(players_per_strategy: np.array, nb_phases: int, **kwargs) -> float:
    # Define matrix that extracts number of communicative players
    communicative_matrix = np.hstack([np.ones(nb_phases), np.zeros(nb_phases)])
    # communicatative have value < nb_phases
    num_communicative = np.matmul(players_per_strategy, communicative_matrix)
    return num_communicative


# %%
def combine_communicative_noncommunicative(players_per_strategy: np.array, nb_phases: int, **kwargs) -> np.array:
    # Combine communicatative (value < nb_phases)
    # and non-communicatative (values >= nb_phases)
    combine_com_noncom_matrix = np.block([[np.eye(nb_phases)],[np.eye(nb_phases)]])
    phase_idxs = np.matmul(players_per_strategy, combine_com_noncom_matrix)
    return phase_idxs


# %%
from enum import Enum
GameType = Enum('GameType',
                ["all communicative", "all noncommunicative",
                 "prisoner's dilemma", "snowdrift", "coordination", "mutualism",
                 "unknown"])

def extract_most_common_game_types(players_per_strategy: np.array,
                   mutual_benefit_synchronous: float,
                   unilateral_benefit_synchronous: float,
                   cost: float,
                   nb_phases: int,
                   **kwargs) -> np.array:
    # Check if all players were communicative/noncommunicative
    nb_communicative = extract_num_communicative(players_per_strategy, nb_phases=nb_phases, **kwargs)
    if nb_communicative == 0:
        return GameType["all noncommunicative"]
    elif nb_communicative == nb_players:
        return GameType["all communicative"]

    # Combine communicative and noncommunicative strategies
    players_per_phase = combine_communicative_noncommunicative(players_per_strategy, nb_phases)

    # Define payoff response submatrices
    phi = np.arange(nb_phases)
    phi_i, phi_j = np.meshgrid(phi,phi)
    delta_phi = phi_j - phi_i
    reward_submatrix = mutual_benefit_synchronous*(1-np.cos(delta_phi))/2 - cost
    sucker_submatrix = unilateral_benefit_synchronous*(1-np.cos(delta_phi))/2 - cost
    temptation_submatrix = unilateral_benefit_synchronous*(1-np.cos(delta_phi))/2
    punishment_submatrix = np.zeros([nb_phases, nb_phases])
    # Define mask matrix for each game type
    prisoners_dilemma_matrix = (temptation_submatrix >= reward_submatrix) \
        & (reward_submatrix >= punishment_submatrix) \
        & (punishment_submatrix >= sucker_submatrix)
    snowdrift_matrix = (temptation_submatrix >= reward_submatrix) \
        & (reward_submatrix >= sucker_submatrix) \
        & (sucker_submatrix >= punishment_submatrix)
    coordination_matrix = (reward_submatrix >= temptation_submatrix) \
        & (temptation_submatrix >= sucker_submatrix) \
        & (sucker_submatrix >= punishment_submatrix)
    mutualism_matrix = (reward_submatrix >= temptation_submatrix) \
        & (temptation_submatrix >= punishment_submatrix) \
        & (punishment_submatrix >= sucker_submatrix)

    # Calculate number of games for each game type
    nb_prisoners_dilemma = np.linalg.multi_dot(
        [players_per_phase, prisoners_dilemma_matrix, players_per_phase])
    nb_snowdrift = np.linalg.multi_dot(
        [players_per_phase, snowdrift_matrix, players_per_phase])
    nb_coordination = np.linalg.multi_dot(
        [players_per_phase, coordination_matrix, players_per_phase])
    nb_mutualism = np.linalg.multi_dot(
        [players_per_phase, mutualism_matrix, players_per_phase])
    
    # Find most common game type
    game_counts = {GameType["prisoner's dilemma"]: nb_prisoners_dilemma,
                   GameType["snowdrift"]: nb_snowdrift,
                   GameType["coordination"]: nb_coordination,
                   GameType["mutualism"]: nb_mutualism}
    most_common_game_type = min(game_counts, key=game_counts.get)
    
    return most_common_game_type


# %%
def extract_phases(players_per_phase: np.array, nb_phases: int, **kwargs) -> np.array:
    # Convert phase idxs to phases
    phases = 2*np.pi/nb_phases*players_per_phase
    return phases

def extract_order_parameters(players_per_strategy: np.array, **kwargs) -> np.array:
    # Combine communicative and noncommunicative
    phase_indxs = combine_communicative_noncommunicative(players_per_strategy, **kwargs)
    # Convert strategy idxs to phases
    phases = extract_phases(phase_indxs, **kwargs)
    order_parameters = np.abs(np.mean(np.exp(1j*phases)))
    return order_parameters


# %% [markdown]
# ## Play the Game

# %%
def run_and_get_all_populations(time_steps, selection_strength, initial_population, **kwargs) -> list:
    mp = create_moran_game(turns=time_steps, selection_strength=selection_strength, **kwargs)
    egt.Random.init(1)
    results = mp.run(nb_generations = time_steps,
           beta = selection_strength,
           init_state = initial_population)
    return results


# %%
def run_and_get_cumulative_populations(**kwargs) -> np.array:
    results = run_and_get_all_populations(**kwargs)
    cumulative_populations = np.sum(results, 0)
    return cumulative_populations


# %% [markdown]
# ## Long time limit

# %%
# Define parameters
time_steps = int(2E9)
time_steps = int(2E4)

# %%
# Define Bs on which to run
B_crit = 2 * cost * (nb_players-1) / (nb_players-2)
nb_Bs = 11
step_size_Bs = 0.04
Bs = B_crit + (np.arange(nb_Bs) - (nb_Bs-1)/2)*step_size_Bs

# %%
# Run the model for weak selection strength
cumulative_populations_weak = [run_and_get_cumulative_populations(
        initial_population = initial_population,
        mutual_benefit_synchronous = B,
        unilateral_benefit_synchronous = 0.9*B,
        selection_strength = selection_strength_weak,
        time_steps = time_steps,
        nb_players = nb_players,
        nb_strategies = nb_strategies,
        cost = cost)
    for B in Bs]

# %%
# Plot fraction communcative
nb_communicative_weak = [extract_num_communicative(final_population, nb_phases) for final_population in cumulative_populations_weak]
fraction_communicative_weak = np.array(nb_communicative_weak)/(time_steps*nb_players)
plt.plot(Bs, fraction_communicative_weak, label="Simulation")
plt.xlabel(r"Maximum benefit of mutual communication, $B(0)$")
plt.ylabel("Frequency of communicative strategies")
plt.title(r"Weak selection $\delta = 0.005$")
plt.show()

# %%
# Run the model for weak selection strength
cumulative_populations_strong = [run_and_get_cumulative_populations(
        initial_population = initial_population,
        mutual_benefit_synchronous = B,
        unilateral_benefit_synchronous = 0.9*B,
        selection_strength = selection_strength_strong,
        time_steps = time_steps,
        nb_players = nb_players,
        nb_strategies = nb_strategies,
        cost = cost)
    for B in Bs]

# %%
# Plot fraction communcative
nb_communicative_strong = [extract_num_communicative(final_population, nb_phases) for final_population in cumulative_populations_strong]
fraction_communicative_strong = np.array(nb_communicative_strong)/(time_steps*nb_players)
plt.plot(Bs, fraction_communicative_strong, label="Simulation")
plt.xlabel(r"Maximum benefit of mutual communication, $B(0)$")
plt.ylabel("Frequency of communicative strategies")
plt.title(r"Strong selection $\delta = 0.2$")
plt.show()

# %% [markdown]
# ## Time-Evolution

# %%
# Source: https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html

from matplotlib.collections import LineCollection
def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_color(c)  # set the colors of each segment

    return ax.add_collection(lc)


# %%
from matplotlib.colors import to_rgb
# Define color map
color_map = {GameType["all noncommunicative"]: to_rgb("darkgrey"),
             GameType["all communicative"]: to_rgb("lightgrey"),
             GameType["mutualism"]: to_rgb("green"),
             GameType["coordination"]: to_rgb("blue"),
             GameType["snowdrift"]: to_rgb("yellow"),
             GameType["prisoner's dilemma"]: to_rgb("red")}

# %%
# Run the model for strong selection strength
B_timeseries = 1.5*cost
all_population_strong = run_and_get_all_populations(
        mutual_benefit_synchronous = B_timeseries,
        unilateral_benefit_synchronous = 0.9*B_timeseries,
        selection_strength = selection_strength_strong,
        time_steps = time_steps,
        nb_players = nb_players,
        nb_strategies = nb_strategies,
        initial_population = initial_population,
        cost = cost)

# %%
# Extract results
nb_communcative = np.apply_along_axis(extract_num_communicative, -1, all_population_strong, nb_phases=nb_phases)
fraction_communicative = np.array(nb_communcative)/(time_steps*nb_players)
most_common_game_types = np.apply_along_axis(extract_most_common_game_types, -1, all_population_strong, mutual_benefit_synchronous = B_timeseries, \
        unilateral_benefit_synchronous = 0.9*B_timeseries, cost = cost, nb_phases = nb_phases)
order_parameters = np.apply_along_axis(extract_order_parameters, -1, all_population_strong, nb_phases = nb_phases)

# Create array of times
# Note: the populations include the initial data, so we need one more than time-steps
times = np.arange(time_steps+1)

# Plot fraction communcative
# Create a figure and plot the line on it
fig1, ax1 = plt.subplots()
lines = colored_line(times, fraction_communicative, [color_map[game_type] for game_type in most_common_game_types], ax1, label="frequency communicative")
ax1.plot(times, order_parameters, label="Order parameter")
fig1.colorbar(lines)  # add a color legend
ax1.set_xlabel(r"Maximum benefit of mutual communication, $B(0)$")
ax1.set_ylabel("Frequency of communicative strategies")
ax1.set_title(r"Strong selection $\delta = 0.2$")
fig1.legend()
plt.show()

import pandas as pd
# Convert GameType to string
most_common_game_names = [game.name for game in most_common_game_types]
# Plot game type histogram
df = pd.DataFrame.from_dict(Counter(most_common_game_names), orient="index")
df.plot(kind="bar")
plt.title(r"Strong selection $\delta = 0.2$")
plt.show()

# %%
# Run the model for very strong selection strength
B_timeseries = 1.5*cost
selection_strength_very_strong = 5
all_population_very_strong = run_and_get_all_populations(
        mutual_benefit_synchronous = B_timeseries,
        unilateral_benefit_synchronous = 0.9*B_timeseries,
        selection_strength = selection_strength_very_strong,
        time_steps = time_steps,
        nb_players = nb_players,
        nb_strategies = nb_strategies,
        initial_population = initial_population,
        cost = cost)

# %%
# Extract results
nb_communcative = np.apply_along_axis(extract_num_communicative, -1, all_population_very_strong, nb_phases=nb_phases)
fraction_communicative = np.array(nb_communcative)/(time_steps*nb_players)
most_common_game_types = np.apply_along_axis(extract_most_common_game_types, -1, all_population_very_strong, mutual_benefit_synchronous = B_timeseries, \
        unilateral_benefit_synchronous = 0.9*B_timeseries, cost = cost, nb_phases = nb_phases)
order_parameters = np.apply_along_axis(extract_order_parameters, -1, all_population_very_strong, nb_phases = nb_phases)

# Create array of times
# Note: the populations include the initial data, so we need one more than time-steps
times = np.arange(time_steps+1)

# Plot fraction communcative
# Create a figure and plot the line on it
fig1, ax1 = plt.subplots()
lines = colored_line(times, fraction_communicative, [color_map[game_type] for game_type in most_common_game_types], ax1, label="frequency communicative")
ax1.plot(times, order_parameters, label="Order parameter")
fig1.colorbar(lines)  # add a color legend
ax1.set_xlabel(r"Maximum benefit of mutual communication, $B(0)$")
ax1.set_ylabel("Frequency of communicative strategies")
ax1.set_title(r"Very strong selection $\delta = 5$")
fig1.legend()
plt.show()

import pandas as pd
# Convert GameType to string
most_common_game_names = [game.name for game in most_common_game_types]
# Plot game type histogram
df = pd.DataFrame.from_dict(Counter(most_common_game_names), orient="index")
df.plot(kind="bar")
plt.title(r"Very strong selection $\delta = 5$")
plt.show()

# %% [markdown]
# # Compare to Theory

# %%
from sympy.abc import i, j, k, q, r

# Define symbols
B_0, beta_0, phi, delta, c = sp.symbols(r'B_0 \beta_0 \phi \delta c')

# Define constants
d = nb_phases
n = nb_players

# Define functions
B = B_0*(1+sp.cos(phi))/2
beta = beta_0*(1+sp.cos(phi))/2
rho_CC = 1/(1 + sp.Sum(sp.exp(delta/(n-1)*(
        k*(k+1)/2*(2*B-2*B_0) + k*(B_0*n-B*n)))
    , (k,1,n-1)))

# Note we've fixed typos in the coefficient of j*beta_0 was fixed: (n-1) -> n and 1/2*j*B_0 -> j*B_0
rho_NC = 1/(1 + sp.Sum(sp.exp(delta/(n-1)
                             *(k**2*(beta-1/2*B_0)
                               +k*(B_0-n*beta+(n-1)*c))), (k, 1, n-1)))
rho_CN = rho_NC*(sp.exp(delta*(n-1)*c-(n-2)/2*B_0))

# Define matrices
B_1 = sp.FunctionMatrix(d, d,
                     sp.Lambda((i,j),
                            rho_CC.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2)-(i-j)))
                                   if i != j
                                   else
                                   1 - sp.Sum(rho_CC.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2) - delta_phi_prime)), (delta_phi_prime, 1, nb_players) )
                           ))
B_2 = sp.FunctionMatrix(d, d,
                     sp.Lambda((i,j),
                            rho_CN.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2)-(i-j)))
                                   if i != j
                                   else
                                   1 - Sum(rho_CN.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2) - delta_phi_prime)), (delta_phi_prime, 1, nb_players) )
                           ))
B_3 = sp.FunctionMatrix(d, d,
                     sp.Lambda((i,j),
                            rho_NC.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2)-(i-j)))
                                   if i != j
                                   else
                                   1 - Sum(rho_NC.subs(phi, sp.floor(d/2) - sp.functions.Abs(sp.floor(d/2) - delta_phi_prime)), (delta_phi_prime, 1, nb_players) )
                           ))
B_4 = sp.FunctionMatrix(d, d,
                     sp.Lambda((i,j), 1/n))

# Define full Markov transition matrix
M = sp.BlockMatrix([[B_1,B_2],[B_3,B_4]]).as_explicit()

# Define stationary distribution
s_1 = sp.Sum(M[r,q], (r, d+1 - 1, 2*d - 1))/(d*sp.Sum(M[q,r] + M[r,q], (r, d+1 - 1, 2*d - 1)))
s_2 = sp.Sum(M[q,r], (r, d+1 - 1, 2*d - 1))/(d*sp.Sum(M[q,r] + M[r,q], (r, d+1 - 1, 2*d - 1)))
s = sp.Matrix([s_1*sp.ones(d,1), s_2*sp.ones(d,1)])

## Check that the explicitly calculated stationary distribution matches the eigenvector
#eigen = M.eigenvects()[0][2][0]
#if s != eigen:
#    raise("Explicitly calculated stationary distribution differs from eigenvector")

from  sympy.plotting.plot import MatplotlibBackend, Plot
# Source: https://stackoverflow.com/a/70486636
def get_sympy_subplots(plot: Plot):
    backend = MatplotlibBackend(plot)
    backend.process_series()
    backend.fig.tight_layout()
    return backend.fig, backend.ax[0]
# Plot communicative frequency vs maximum benefit
# Note: it can be shown that s is independent of q, so substitute any value
s_1_eval = s_1.subs(beta_0, 0.95*B_0).subs(c, cost).subs(q, 0)

# Weak selection
p = sp.plot(s_1_eval.subs(delta, selection_strength_weak).doit(),
     xlim = (0, 0.5), ylim = (0, 1), axis_center = (0,0), title='Weak selection',
     label = "Theory", show = False)

fig, ax = get_sympy_subplots(p)
ax.plot(Bs, fraction_communicative_weak, label="Simulation")
fig.legend()
fig.show()
# Strong selection
p = sp.plot(s_1_eval.subs(delta, selection_strength_strong).doit(),
     xlim = (0, 0.5), ylim = (0, 1), axis_center = (0,0), title='Strong selection',
     label = "Theory", show = False)
fig, ax = get_sympy_subplots(p)
ax.plot(Bs, fraction_communicative_strong, label="Simulation")
fig.legend()
fig.show()
