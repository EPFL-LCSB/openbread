openbread
=========

|Build Status| |Codacy branch grade| |license|

A python 3 implementation of OPENFPM based brownian reaction dynamics.


Container-based install
-----------------------

We recommend to use this package inside of the provided DOCKER container. This requires an installation of Docker (R)
(Community edition sufficient) that is freely available |here|_

See |docker|_ page for the setup and use of the container.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/openbread/tree/master/docker

.. |here| replace:: ``here``
.. _here: https://www.docker.com/`


Source-based install
-----------------------


For a source based it is necessary to install |OPENFPM|_ with CXXFLAGS=-fPIC CFLAGS=-fPIC.
We recommend to use our fork of the openfpm package that is using openmpi 3.0.0. If  |OPENFPM|_ is compiled with
openmpi 2.X.X openbread cannot be compiled due to a bug impacting shared libraries.

.. code:: bash

    git clone https://github.com/weilandtd/openfpm_pdata.git
    cd openfpm_pdata
    ./install -s -i "/installation/directory" \
                 -c "--prefix=/dependency/directory  CXXFLAGS=-fPIC  CFLAGS=-fPIC"
    make install

Note that correct include and library paths need to be added. Be aware to specify your python version 3.X for the
PYTHON_INCLUDE environmental variable:

.. code:: bash

    cp /path/to/openfpm/ ~/example.mk
    source /openfpm_install_dir/openfpm_vars
    export PYTHON_INCLUDE=$PYTHON_ROOT/include/python3.Xm

The package can then be installed using the distutils:

.. code:: bash

    pip install -e /path/to/openbread

The source-based install of openbread was only tested in Linux based environments. Platform compatibility is achieved
using the container-based install.


Requirements
------------

This module was developed in Python 3.5, and it is recommended to run Python 3.5.
The module also was tested in Python 3.6.

This module requires |OPENFPM|_ to work properly.

.. |OPENFPM| replace:: ``OPENFPM``
.. _OPENFPM: http://openfpm.mpi-cbg.de/install_from_source#intro-wrapper

Further the following pip-python packages are required
    - sympy >= 1.1.
    - pytest
    - scipy
    - numpy
    - dill

Quick start
===========

.. code:: python

    from openbread.core import *


    # Construct species in the model

    A = Species(    name = 'A',
                    diffusion_constant = 100e-12,
                    collision_radius   = 2e-9,
                    mass = 1 )

    B = Species(    name = 'B',
                    diffusion_constant = 100e-12,
                    collision_radius   = 2e-9,
                    mass = 1 )

    C = Species(    name = 'C',
                    diffusion_constant = 50e-12,
                    collision_radius   = 4e-9,
                    mass = 2 )


    # Rate constant
    k = 1e12
    Keq = 50e-6

    # Setup particle simulation environemnt
    volume = 10e-18 # (0.1 mum)^3 in L

    medium = ParticleModel.Medium(  viscosity=0.7e-3, # Pa s
                                    temperatur=310.15)
    volume_fraction = 0.0

    crowding = ParticleModel.Crowding( volume_fraction = volume_fraction,
                                       mu = np.log(31.9),
                                       sigma = 0.825,
                                       max_size = 10e-3)

    particle_model = ParticleModel(medium,
                                   crowding,
                                   volume)

    particle_model.add_reaction(Reaction('A+B->C', {A:-1,B:-1,C:1},  k ))
    particle_model.add_reaction(Reaction('C->A+B', {A:1,B:1,C:-1},   k*Keq ))

    # Define initial conditions
    particle_model.initial_conditions['A'] = 50e-6
    particle_model.initial_conditions['B'] = 50e-6


    result = particle_model.simulate(   dt=1e-9,
                                        max_time=1e-5,
                                        log_step=10,
                                        random_seed=1,
                                        is_hardsphere=False,
                                        is_constant_state=False,
                                        t_equlibriate=0.0)



License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the
`LICENSE <https://github.com/EPFL-LCSB/openbread/blob/master/LICENSE.txt>`_ file for more details



.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/openbread/blob/master/LICENSE.txt
   
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/openbread.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/openbread
   
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/a914add308be46418e4e1ecc9d197317/master.svg
   :target: https://www.codacy.com/app/realLCSB/openbread
   
   

