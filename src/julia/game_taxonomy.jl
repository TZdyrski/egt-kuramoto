@enum GameType harmony chicken battle hero compromise concord staghunt dilemma deadlock assurance coordination peace neutral

@enum TieType low mid high double triple basic zero

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

@enum StrategyParity all_communicative all_noncommunicative disconnected_synchronized_populations mixed
