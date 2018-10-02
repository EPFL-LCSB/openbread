# -*- coding: utf-8 -*-
"""
.. module:: openbread
   :platform: Unix, Windows
   :synopsis: Python Multiscale Enzyme Simulator

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


# Construct species in the model

g3p = Species(  name = 'g3p',
                diffusion_constant = 940e-12,
                collision_radius   = 1.11e-9,
                mass = 0.186 )
g2p = Species(  name = 'g2p',
                diffusion_constant = 940e-12,
                collision_radius   = 1.11e-9,
                mass = 0.186 )
pgm = Species(  name = 'pgm',
                diffusion_constant = 84.8e-12 ,
                collision_radius   = 3.87e-9,
                mass = 61 )
EC_pgm = Species(  name = 'EC_pgm',
                   diffusion_constant = 84.8e-12,
                   collision_radius   = 3.87e-9,
                   mass = 61.186 )

species = [g3p, g2p, pgm, EC_pgm]

# Define microscopic reaction rate constants:
k1f = 1e5       # 1/Ms
k1b = 20.0      # 1/s
k2f = 10.0      # 1/s
k2b = 1e5       # 1/Ms


# Setup particle simulation environemnt
volume = 10e-18 # (0.1 mum)^3 in L

medium = ParticleModel.Medium(  viscosity=0.7e-3, # Pa s
                                temperatur=310.15)
volume_fraction = 0.1

crowding = ParticleModel.Crowding( volume_fraction = volume_fraction,
                                   mu = np.log(31.9),
                                   sigma = 0.825,
                                   max_size = 10e-3)

particle_model = ParticleModel(medium,
                               crowding,
                               volume)

particle_model.add_reaction(Reaction('E+S->ES', {g2p:-1,pgm:-1,EC_pgm:1},  k1f ))
particle_model.add_reaction(Reaction('ES->E+S', {g2p: 1,pgm: 1,EC_pgm:-1}, k1b ))
particle_model.add_reaction(Reaction('ES->E+P', {g3p: 1,pgm: 1,EC_pgm:-1}, k2f ))
particle_model.add_reaction(Reaction('E+P->ES', {g3p:-1,pgm:-1,EC_pgm:1},  k2b ))


# Define initial conditions
particle_model.initial_conditions['pgm'] = 50e-6
particle_model.initial_conditions['EC_pgm'] = 50e-6
particle_model.initial_conditions['g3p'] = 50e-6
particle_model.initial_conditions['g2p'] = 50e-6



result = particle_model.simulate(   dt=1e-9,
                                    max_time=1e-7,
                                    log_step=10,
                                    random_seed=1,
                                    is_hardsphere=True,
                                    is_constant_state=False,
                                    t_equlibriate=0.0)



