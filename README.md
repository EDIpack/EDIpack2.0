#  Lanczos Exact Diagonalization for Quantum Impurity problems 

A Lanczos based solver for Quantum Impurity problems with a special application to Dynamical Mean-Field Theory. 

This library provides a simple, yet generic, interface to the solution of multi-orbital quantum impurity models with a finite, discrete electronic bath. The solution is obtained in either for *normal*, *superconducting* (s-wave) or *Sz-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases. The code works at zero and finite, but low, temperature. 



### Dependencies

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

  


### Installation

Installation is  available using CMake. In the current v0.0.1 API are only provided in Fortran.   

Clone the repo:

`git clone https://github.com/aamaricci/lib_dmft_ed scifor`

And from the repository directory (`cd lib_dmft_ed`) make a standard out-of-source CMake compilation:

`mkdir build`
`cd build`
`cmake ..`     
`make`     
`make install`   
`make post-install`    

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/dmft_ed.pc`  
* environment module file `~/.modules.d/dmft_ed/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/dmft_ed/VERSION/PLAT/[GIT_BRANCH]>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright 2020- (C) Adriano Amaricci

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.