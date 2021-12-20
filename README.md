# Virtual_Bioreactors

#### Authors:
#### Daniel Harris — Brown University
#### Radu Cimpeanu — University of Warwick
#### Benny Smith — Brown University
#### Peyton Newman — Brown University

#### This repository hosts code written using [Basilisk C](http://basilisk.fr) to simulate a [Rocking Wave bioreactor](https://link.springer.com/article/10.1023/A:1008025016272). The goal is to support the Cultivated Meat industry by providing the ability to simulate bioreactors with a variety of setups in order to study shear stress and mixing processes. Thanks to [The Good Institute](https://gfi.org) for funding this project.
#### There are two main components of this work: fluid dynamical simulations of a [rocking ellipse](https://github.com/austinbennysmith/Virtual_Bioreactors/tree/main/Rocking_Ellipse) containing water and air, as well as an advection-diffusion solver to keep track of a tracer representing oxygen in the reactor. Originally, we attempted to come up with our own [advectiond-diffusion discretization](https://github.com/austinbennysmith/Virtual_Bioreactors/tree/main/advection_diffusion). It worked well for gentle parameters (i.e. water & oil), but broke down when parameters were adjusted to represent air. Therefore, the [latest version](https://github.com/austinbennysmith/Virtual_Bioreactors/tree/main/Rocking_Ellipse/Builtin_Solver/FullSized) of the rocking ellipse simulation just uses built-in [advection](http://basilisk.fr/src/advection.h) and [diffusion](http://basilisk.fr/src/diffusion.h) functions.

### **Directions for running code yourself:**
