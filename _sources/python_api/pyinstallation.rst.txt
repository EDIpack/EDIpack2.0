.. _edipy_install:

EDIpy2 Install
==============


Anaconda
------------

We provide Linux and MacOS packages for the `Anaconda <https://www.anaconda.com/>`_ distribution. Packages are available for Python version 3.10 and later.
Once a command-line tool such as `conda <https://www.anaconda.com/>`_ or `mamba <https://mamba.readthedocs.io/en/latest/>`_ is installed, an environment using one of the available python version can be created, and then the EDIpack2.0 package can be installed:

.. code-block:: shell

   conda create -n edipack
   conda activate edipack
   conda install -c conda-forge -c edipack edipack2


the python module `edipy2` can then be directly imported.

Compile from source
---------------------

The python module `edipy2` requires:

* `SciFortran <https://github.com/scifortran/SciFortran>`_

* `EDIpack2 <https://github.com/edipack/EDIpack2.0>`_

libraries to be installed beforehand, please see related documentation
to install such Fortran libraries. Once both are set up, the python module can be installed from the EDIpack2 folder via

.. code-block:: shell

   pip install . --break-system-packages
   
The latter option may not be required in all cases, but it is in recent versions of Debian and OSX. Since no edipy2 package is provided by any distro, this will not create problems. If the user is using a virtual environment, the option is not necessary.





