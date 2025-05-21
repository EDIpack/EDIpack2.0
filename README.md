# EDIpack: a generic and interoperable Lanczos-based  exact diagonalization solver for Quantum Impurity problems 

[![TestSuite](https://img.shields.io/github/actions/workflow/status/edipack/EDIpack/PushWorkflow.yml?label=TestSuite&logo=Fortran&style=flat-square)](https://github.com/edipack/EDIpack/actions/workflows/PushWorkflow.yml) 
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://edipack.github.io/EDIpack/)
[![Anaconda-Server Badge](https://anaconda.org/edipack/edipack/badges/version.svg)](https://anaconda.org/edipack/edipack)

<!-- TO BE SETUP ASAP
[![Coverage]()]()
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/DMFT_ED)
-->

This is the latest version of [EDIpack](https://github.com/edipack/EDIpack): a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
 


### Install 
*EDIpack* is available in the form of a static Fortran library (`libedipack.a`) and the related Fortran module `EDIPACK`.
The release version includes additional modules to extend the software functionalities: i) an inequivalent impurities extension `Edipack2ineq`
and ii) a shared dynamical library `edipack_cbinding.so` implementing the Fortran-C interface. 

A standard installation from source is available through `CMake`, via the standard out-of-source method. 

An alternative approach is provided via `Anaconda`. 

Detailed information can be found at [edipack.github.io/EDIpack/installation](https://edipack.github.io/EDIpack/installation.html)



### Documentation
All the informations about the structure of the library and its use can be found in the documenation at [edipack.github.io/EDIpack/](https://edipack.github.io/EDIpack/)  


### Use
In [Quickstart](https://edipack.github.io/EDIpack/quickstart/02_dmft.html) we illustrate the use and the capabilities of EDIpack as a solver for Dynamical Mean-Field Theory calculation. 



### Authors
[Adriano Amaricci](https://github.com/aamaricci)  
[Lorenzo Crippa](https://github.com/lcrippa)  
[Samuele Giuli](https://github.com/SamueleGiuli)  
[Gabriele Bellomia](https://github.com/beddalumia)  
[Giacomo Mazza](https://github.com/GiacMazza)  
[Francesco Petocchi](mailto:francesco.petocchi@gmail.com)  
[Massimo Capone](mailto:capone@sissa.it)

### Reference
*A specific work is on the way to be published.*  
In the meantime you can find most of the information about EDIpack in the documetation page [edipack.github.io/EDIpack/](https://edipack.github.io/EDIpack/).
Useful information about the distributed memory parallel algorithms underlying the functioning of EDIpack are presented in [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261). 
Should you use EDIpack or any of the derived packages, please consider citing both papers.

### Issues
If you encounter bugs or difficulties, please [file an issue](https://github.com/edipack/EDIpack/issues/new/choose). For any other communication, please reach out any of the developers.          
