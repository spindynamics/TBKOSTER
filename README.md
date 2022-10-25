# DyNaMol

DyNaMol is an open source software simulation package for the simulation of magnetic materials developed at the [French Alternative Energies and Atomic Energy Commission](http://www.cea.fr). DyNaMol is designed to provide a community standard tool for ab initio simulations of magnetic materials with high performance and an easy-to-use interface.

Using a variety of common simulation methods it can calculate the equilibrium and dynamic magnetic properties of a wide variety of magnetic materials and phenomena, including ferro-, ferri- and antiferro-magnets, nano objects, ultrafast spin dynamics, magnetic recording media, exchange bias, magnetic multilayer films and even complete devices.

DyNaMol is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/ or redistribute the software under the terms of the [CeCILL license](http://www.cecill.info) as circulated by CEA, CNRS and INRIA.

To account DyNaMol in your developments, please cite and refer to the following references for more details :
* [Ramon Cardias, Cyrille Barreteau, Pascal Thibaudeau and Chu Chun Fu, Phys. Rev. B **103**, 235436 (2021)](https://doi.org/10.1103/PhysRevB.103.235436) and [arXiv.2101.06121](https://arxiv.org/abs/2101.06121) 
* [Cyrille Barreteau, Daniel Spanjaard and Marie-Catherine Desjonquères, Comptes Rendus Physique **17** (3-4) 406-429 (2016)](https://www.sciencedirect.com/science/article/pii/S1631070515002601?via%3Dihub)


## Getting DyNaMol
DyNaMol is available from github. The code runs in serial on MacOS and linux. An OpenMP parallelism is deployed on k-points. A comprehensive overview of the software features and example input files are also provided.

## Capabilities
DyNaMol is designed to be highly flexible to deal with a wide variety of problems using a diverse set of simulation tools and methods. The capabilities of the code can be summarized broadly in terms of the simulation methods, standard problems, structural properties and features of the code, all of which can be combined to tackle almost any problem.

### Simulation methods
* -[x] Stochastic Landau-Lifshitz-Gilbert
* -[x] Direct Energy minimization of a spin system

### Standard calculations
* -[x] Tight-Binding total energy and atomic forces
* -[x] Atomic magnetic moments in non-colinear electronic configurations
* -[x] Magnetic effective fields for Ultrafast spin dynamics

### Structural properties
* -[x] Bulk-like systems
* -[x] Thin films
* -[x] Nanoparticles - spheres, cubes, truncated octahedra, cylinders
* -[x] Nanoparticle arrays
* -[x] Core-shell nanoparticles
* -[x] Multilayer thin films
* -[x] Interface roughness and intermixing
* -[x] Dilute magnetic systems
* -[x] full crystal structures
* -[x] User-defined atomic structures - for example from Molecular Dynamics simulations

### Magnetic properties
* -[x] Ferromagnets
* -[x] Antiferromagnets
* -[x] Ferri-magnets
* -[x] Spin glass
* -[x] User-defined Tight-Binding Hamiltonian from ab-initio Density Functional Theory (DFT) calculations

### Code features
* -[x] Modular Fortran code and shell script
* -[x] High performance code
* -[x] Parallelisation using OpenMP threads
* -[x] Output for visualisation and publication quality graphics
* -[x] Minimal dependence on external libraries for portability
* -[x] Freely available open source code

## Contributors

* [Pascal Thibaudeau](https://github.com/pthibaud)
* [Cyrille Barreteau](https://github.com/CyrilleBarreteau)
* [Mathieu César](https://github.com/MathieuCesar)
* [Ramon Cardias](https://github.com/ramoncardias)

## Installation

There are many ways to install DyNaMol.
The simplest way is to use Docker/Podman with the Dockerfile. 
```bash
docker build -t dynamol .
docker run -it localhost/dynamol:latest
```
The second way is to install the necessary packages locally.
The installation process relies heavily on [CMake](https://cmake.org). CMake is used to control the software compilation process using simple platform and compiler independent configuration files, and generate native makefiles and workspaces that can be used in the compiler environment of your choice.
Please read the corresponding manual. For further details, a UserManual will be created during the compilation process if a LaTeX distribution is found.  
