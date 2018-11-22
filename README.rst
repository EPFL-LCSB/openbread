openbread
=========

A python 3 implementation of OPENFPM based brownian reaction dynamics.


Requirements
------------

**This module was developed in Python 3.5, and it is recommended to run Python 3.5

This module requires `OPENFPM http://openfpm.mpi-cbg.de/`_ to work properly.

Container-based install
-----------------------

We recommend to use this package inside of the provided DOCKER container. This requieres an installation of Docker (R)
(Community editions sufficient) that is freely available `here https://www.docker.com/`_

See |docker|_ page for the setup and use of the container.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/openbread/tree/master/docker

To

Source-based install
-----------------------

For a source based install a working installation of |OPENFPM|_ is required.
And the variables of the correct include and library paths are added using:

.. code:: bash
    source /openfpm_install_dir/openfpm_vars


.. |OPENFPM| replace:: ``openfpm``
.. _OPENFPM: http://openfpm.mpi-cbg.de/install_from_source#intro-wrapper


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
`LICENSE <https://github.com/EPFL-LCSB/openbread/blob/master/LICENSE.txt>`_ file for more details.
