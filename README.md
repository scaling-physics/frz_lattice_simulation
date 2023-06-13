# FRZ Lattice Simulation

This project is a lattice simulation for studying the behavior of interacting particles on a 2D hexagonal lattice. The simulation is implemented on C++ and provides a framework for analyzing the properties and dynamics of the lattice system. 

## Installation

To use this simulation, follow these steps:

1. Clone the repository: `git clone git@gitlab.gwdg.de:murray-group/frz_lattice_simulation.git`

2. Compile main.cpp. 

3. There are two important branches. The master branch has the version of the simulation without FrzB particles. The two_particle_model branch contains the updated model in th presence of FrzB particles. 

4. Execute binary file with the following inputs: `./frz_lattice_model $J $alpha $FrzB_num $beta $index`. J=J_AA and beta = J_AB represent the interaction between AA and AB particles, FrzB_num: Number of FrzB particles, index: simulation instance.


## License

This project is licensed under the [MIT License](LICENSE.md). See the `LICENSE.md` file for more details.


