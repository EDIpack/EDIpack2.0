.. _edipack_install:

Install
#####################

|edipack2| is available in the form of a static Fortran library
(`libedipack2.a`) and the related Fortran module `EDIPACK2`.
The release version includes additional packages to extend the
software functionalities: an inequivalent impurities extension |edipack2ineq|
and two shared dynamical libraries implementing the
Fortran-C interface. 

A standard installation method is available through CMake.

An alternative approach is provided via Anaconda. 



Compiling from source
======================

Building
---------

We assume that `SciFortran` and `MPI` have been correctly installed
and are available in the system. See related documentation. Note that
the installation of |edipack2| closely follows the `SciFortran`
template.


Clone the repo:

.. code-block:: bash
		
   git clone https://github.com/edipack/EDIpack2 EDIpack2



Optionally define the fortran compiler:

.. code-block:: bash
		
   export FC=mpif90/gfortran/ifort


From the repository directory (`cd EDIpack2`) make a standard
out-of-source CMake compilation:

**GNU Make**

Using GNU `make` is the default CMake workflow, with widest version
support (CMake > 3.0). Note that parallel `make` execution is tested
and working.

.. code-block:: bash
		
   mkdir build 
   cd build  
   cmake .. 
   make -j



**Ninja**

Using `ninja` if a fortran-capable version of `ninja
<https://ninja-build.org>`_ is available in your system (and CMake can
take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

.. code-block:: bash
		
   mkdir build    
   cd build  
   cmake -GNinja ..  
   ninja

The `CMake` compilation can be customized using the following
additional variables:   

.. list-table:: CMake Options
   :widths: 30 20 50
   :header-rows: 1

   * - Option
     - Scope
     - Value
       
   * - :code:`-DPREFIX`
     - prefix directory  
     - ~/opt/EDIpack2/VERSION/PLATFORM/[GIT_BRANCH]
       
   * - :code:`-DUSE_MPI`
     - MPI support pre-compilation flag
     - yes/True OR no/False (default: True)

   * - :code:`-DVERBOSE`
     - Verbose CMake output 
     - yes/True OR no/False (default: True, superseded by :code:`make VERBOSE=yes/no`

   * - :code:`-DBUILD_TYPE`
     - Compilation flags
     - RELEASE/TESTING/DEBUG/AGGRESSIVE (Default RELEASE)

..
   TESTING:mild or no optimization,  DEBUG:relevant debugging options,  
.. warning::
   
   :code:`BUILD_TYPE=AGGRESSIVE`  includes many deep level debug options which might not compile on some systems or breakdown compilation at linking step.  


These steps build the whole set of libraries provided by |edipack2|. A
specific building for each library is also available specifying the
required target. For user convenience at the end CMake configuration
step the following recap is printed:

.. code-block:: bash

   *Build edipack2 [Default]:  
   $ make -j OR  $ make -j edipack2
   
   *Build edipack2_cbinding C-bindings: 
   $ make edipack2_cbinding
   
   *Build edipack2ineq Inequivalent Sites Extension: 
   $ make edipack2ineq
   
   *Build edipack2ineq_cbinding Inequivalent Sites Extension C-bindings: 
   $ make edipack2ineq_cbinding
   
   *Install: 
   $ make (edipack2/edipack2ineq/edipack2_cbinding/, default=all)   install
   
   *Uninstall: 
   $ make uninstall
   
   *Build documenation: 
   $ make doc
   
   *Build and Runtest: 
   $ make test



   
Installing
------------

System-wide installation is completed after the build step using either:

.. code-block:: bash

   make install

or

.. code-block:: bash
		
   ninja install

  
Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

*  A generated `environment module`_ , installed to`~/.modules.d/EDIpack2/<PLAT>`
  
* A generated `bash` script at `<PREFIX>/bin/configvars.sh`, to be sourced for permanent loading.

*  A generated `pkg-config`_ file to, installed to `~/.pkg-config.d/EDIpack2.pc`  

.. _environment module: https://github.com/cea-hpc/modules
.. _pkg-config: https://github.com/freedesktop/pkg-config


Uninstalling
--------------

Although CMake does not officially provide uninstall procedures in the
generated Make/Ninja files. Hence SciFortran supplies a homebrew
method to remove the generated files by calling (from the relevant
build folder):

.. code-block:: bash
		
   make uninstall

or

.. code-block:: bash
		
   ninja uninstall



Anaconda
======================

We provide Linux and MacOS packages for the `Anaconda <https://www.anaconda.com/>`_ 
distribution. To install the module, the virtual environment of choice should include
python 3.10 or later.

Once a command-line tool such as `conda <https://www.anaconda.com/>`_ or 
`mamba <https://mamba.readthedocs.io/en/latest/>`_ is installed, an environment 
using one of the available python version can be created, and then the EDIpack2.0 
package can be installed:

.. code-block:: shell

   conda create -n edipack
   conda activate edipack
   conda install -c conda-forge -c edipack edipack2


this installs a bundle of the `scifor` and `edipack2` libraries. In order to compile a
fortran program linking the libraries, we provide  `.pc` files which are readable via 
:code:`pkg-config`. If not present, the :code:`compilers` and :code:`pkg-config` conda
packages need to be installed

.. code-block:: shell

   conda install compilers
   conda install pkg-config
   
The inclusion and linking flag can then be obtained via 

.. code-block:: shell

   pkg-config --cflags edipack2 scifor
   pkg-config --libs   edipack2 scifor

