# plants-pollinators

This repository contains supplemental materials for the term project in MATH3620: Dynamical Models in Biology. 

## Abstract. 
Fishman et al 2010 develop a novel model for plant-pollinator population dynamics and explore it analytically, while also relating it to the Beddington-DeAngelis function (1975), a generic model for consumer-resource interactions. The proposed system of ODEs for plant and pollinator population densities exhbits bistable behavior for certain combinations of parameters. I investigate (1) how changes in the parameters affect the opportunity of the system to recover from a dramatic decrease in either population and (2) how sensitive the coexistence equilibrium's location is to plants' and pollinators' ability to "capitalize" on a successful interaction. 

## Citation:
Fishman, M. A., Hadany L. 2010. “Plant–Pollinator Population Dynamics.” Theoretical Population Biology 78 (4): 270–77. https://doi.org/10.1016/j.tpb.2010.08.002.

## In this repository:

- figure-2a-left.R code to reproduce the phase-plane portrait of the system that only has the extinction equilibrium
- figure-2a-right.R code to reproduce the phase-plane portrait of the system that exhibits bistability
- figure-3.R boilerplate code to reproduce the animations for increasing parameters to observe a change in coexistence equilibrium location and the area of the basin of attraction to the coexistence equilibrium
- figure-4a.R code to reproduce a plot of series of coexistence equilibrium locations obtained by incrementally varying parameters, started from the parameter set suggested by Fishman et al 2010. 
- figure-4a.R code to reproduce phase-plane portrait of the bistable system with force field and nullclines. 
- increasing-eta-animation.gif animation of how increase in parameter eta affects the location of the coexistence equilibrium and its basin of attraction 
- increasing-eta-animation.gif animation of how increase in parameter mu affects the location of the coexistence equilibrium and its basin of attraction 
