# structure_dissipation_active_matter
A study into how the energy dissipated when nonconservative forces act on Brownian particles change their structure.

This repository contains all the code necessary to reproduce the relationship between the pair correlation function of an active fluid and the energy dissipation in the liquid, as explained in https://arxiv.org/abs/2012.10441, and verify it using a novel machine learning network. The codes can also be used to derive other relationships between dissipation and structure that we uncovered in the past, such as in https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041026 and https://www.pnas.org/content/115/14/3569.short.

First of all, the lammps source code modifications that are required to run Active Orhnstein Uhlenbeck particle (AOUP) simulations, along with the C++ and bash files used to run simulations, can be found in the folder /simulations. These files enable one to run AOUP simulations at varying fractions of active particles, calculate rate of work and extract the pair correlation function.

Second, the files used to pre-process the simulation data and run the machine learning component can be found in the folder /machine_learning.
