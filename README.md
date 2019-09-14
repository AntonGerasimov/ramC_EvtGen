Program for spectral functions' calculation


# Prerequisites

The program uses the following package:

* ROOT
* EvtGen (including HepMC)
* ramC

These packages should be installed separately and the paths to them
should be added in the corresponding Find*.cmake files

# Instalation

This program is a CMake package, so usual steps are required for
installation:


    cd build/
    cmake ../
	make

As a result several programs will be generated, the main one is
**pion.exe**



# Usage

This program accepts one argument --- the number of the produced
pi-mesons

    ./pion.exe <n>

where <n> should be between 2 and 5. As the result the file
plot<n>.txt will be generated. If the argument of the program is
"all", all spectral functions will be calculated.

# LOG

* Sep 14, 2019: ramC is removed from the project (should be installed separately)
* Sep 14, 2019: Required number of pions can be specified in command line
* Sep 13, 2019: Project created
