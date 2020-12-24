# structure_dissipation_active_matter
A study into how the energy dissipated when nonconservative forces act on Brownian particles change their structure.

This repository contains all the code necessary to reproduce the relationship between the pair correlation function of an active fluid and the energy dissipation in the liquid, as explained in https://arxiv.org/abs/2012.10441.

First of all, the lammps source code modifications that are required to run Active Orhnstein Uhlenbeck particle simulations, along with the C++ file used to run simulations, can be found in the folder /simulations

Second, the files used to run the machine learning component can be found in the folder /machine_learning
