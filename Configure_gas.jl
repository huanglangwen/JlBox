using OrdinaryDiffEq
const to = TimerOutput()
const mechfile="data/BCR_rxn.txt"
const initfile="data/BCR_pop.txt"
const simulation_time= 10.0 # seconds
const tspan=(0,simulation_time)