# Ryegrass-simulation

This repository includes all the materials used for the simulation in the paper:

> Effects of Different Strategies for Exploiting Genomic Selection in Perennial Ryegrass Breeding Programs


### Simulation outline
The simulation study consisted of the following main steps: 
- Simulation of ryegrass base population and initial ryegrass varieties.
- Simulation of conventional breeding and GS schemes in various scenarios.


### Simulation of the base population and initial varieties
The [QMSim](http://animalbiosciences.uoguelph.ca/~msargol/qmsim/) software (Sargolzaei and Schenkel 2009) was used to simulate a historical population. 
The paremeter file of the QMSim (QMSim.par) can be used for the creation of base population.

### Scenarios
Different scenarios were compared to a conventional breeding program. Simulated scenarios differed in the method of selection and structure of the breeding program. Two scenarios (Phen-Y12 and Phen) for phenotypic selection and three scenarios (GS-Y12, GS and GS-SP) were considered for genomic breeding schemes. The R codes used for the simulation are attached. 

To start for the simulation, first QMSim should run with the parameter file (QMSim.par) to create initial varieties. The R codes will use the output of QMSim to do simulate diffrent scenarios.
