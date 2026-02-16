# SPDX-License-Identifier: GPL-3.0-or-later

@enum GameType chicken battle hero compromise deadlock dilemma staghunt assurance coordination peace harmony concord neutral allCommunicative allNoncommunicative disconnectedSynchronizedPopulations

@enum TieType lowTie midTie highTie doubleTie tripleTie basicTie zeroTie

const game_type_full_names = Dict(chicken => "Snowdrift", # Also called chicken
                                  battle => "Battle of the Sexes",
                                  hero => "Hero",
                                  compromise => "Compromise",
                                  deadlock => "Deadlock",
                                  dilemma => "Prisoner's Dilemma",
                                  staghunt => "Mutualism", # Also called Staghunt
                                  assurance => "Assurance",
                                  coordination => "Coordination",
                                  peace => "Peace",
                                  harmony => "Harmony",
                                  concord => "Concord",
                                  neutral => "Neutral",
				  )

# Source: doi:10.3390/g6040495
# In left-up convention (modified from right-up convention by swapping columns)
# = means exactly the same payoff matrix, \approx means equivalent up to swap_strategies!
const game_taxonomy = Dict(
			 [4 2;3 1] => (missing, concord),
			 [4 3;2 1] => (missing, harmony),
			 [4 3;1 2] => (missing, peace),
			 [4 2;1 3] => (missing, coordination),
			 [4 1;2 3] => (missing, assurance),
			 [4 1;3 2] => (missing, staghunt),
			 [3 1;4 2] => (missing, dilemma),
			 [2 1;4 3] => (missing, deadlock),
			 [1 2;4 3] => (missing, compromise),
			 [1 3;4 2] => (missing, hero),
			 [2 3;4 1] => (missing, battle),
			 [3 2;4 1] => (missing, chicken),
			 [2 3;4 2] => (lowTie, battle),
			 [2 2;4 3] => (lowTie, deadlock),
			 [3 2;4 2] => (lowTie, dilemma),
			 [4 2;2 3] => (lowTie, coordination),
			 [4 3;2 2] => (lowTie, harmony),
			 [4 2;3 2] => (lowTie, concord),
			 [3 3;4 1] => (midTie, battle),
			 [1 3;4 3] => (midTie, compromise),
			 [3 1;4 3] => (midTie, deadlock),
			 [4 1;3 3] => (midTie, staghunt),
			 [4 3;1 3] => (midTie, peace),
			 [4 3;3 1] => (midTie, harmony),
			 [1 4;4 2] => (highTie, hero),
			 [2 4;4 1] => (highTie, hero), # high battle \approx high hero
			 [4 2;4 1] => (highTie, concord), # = high chicken
			 [4 1;4 2] => (highTie, staghunt), # = high dilemma
			 [4 2;1 4] => (highTie, coordination),
			 [4 1;2 4] => (highTie, coordination), # high assurance \approx high coord
			 [4 4;1 2] => (highTie, peace),
			 [2 1;4 4] => (highTie, peace), # high deadlock \approx high peace
			 [4 4;2 1] => (highTie, harmony),
			 [1 2;4 4] => (highTie, harmony), # high compromise \approx high harmony
			 [4 2;4 2] => (doubleTie, staghunt), # = double dilemma; note: bruns2015 claims double dilemma = double dilemma, which is tautologically true, but (correctly) states double dilemma = double staghunt elsewhere
			 [4 2;2 4] => (doubleTie, coordination),
			 [2 4;4 2] => (doubleTie, coordination), # double hero \approx double coordination; note: this is the only equivalent pair that is not related by swap_strategies!;
			 # Instead, they are related by swapping only the columns or rows, but not both. Usually, doing this produces a non-symmetric game, but it this case, it is still symmetric
			 [4 4;2 2] => (doubleTie, harmony),
			 [2 2;4 4] => (doubleTie, harmony), # double compromise \approx double harmony
			 [4 4;1 4] => (tripleTie, deadlock),
			 [4 1;4 4] => (tripleTie, deadlock), # note: not explicitly given in bruns2015, but equivalent under swap_strategies!
			 [4 4;4 1] => (tripleTie, harmony),
			 [1 4;4 4] => (tripleTie, harmony), # note: not explicitly given in bruns2015, but equivalent under swap_strategies!
			 [3 3;4 3] => (basicTie, dilemma),
			 [4 3;3 3] => (basicTie, harmony),
			 [4 4;4 4] => (zeroTie, neutral),
			 )
