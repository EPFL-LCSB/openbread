# -*- coding: utf-8 -*-
"""
.. module:: openbread
   :platform: Unix, Windows
   :synopsis: OPENFPM based brownian reaction dynamics

.. moduleauthor:: openbread team

[---------]

Copyright 2017 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

from openbread.core import *



def test_model_build():
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

    assert A.name in particle_model.species
    assert A in particle_model.species.values()


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

    # Test mass conservation
    assert result.species['A'][0]  + result.species['C'][0] == \
           result.species['A'][-1] + result.species['C'][-1]
