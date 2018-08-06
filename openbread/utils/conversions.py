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


from .constants import AVOGADRO_NUMBER
from numpy import pi

# TODO this function should be using a proper unit conversion!!!!

def match_reaction_rate(rate_constant, sum_of_species ):
   # TODO convert rate to IS units
   # Here we assime that the rate constant input was in L/mol/s
   # we convert it to m^3/mol/s
   rate_constant_standard = rate_constant/1000.0
   # Converserion to m^3/s
   rate_constant_standard = rate_constant_standard/AVOGADRO_NUMBER

   # TODO rescaling of the input space
   # ASSUME here that IS units are used
   sum_diffusion = sum_of_species.diffusion_constant
   sum_radius    = sum_of_species.collision_radius
   diffus_limited_rate_constant = 4.0*pi*sum_diffusion*sum_radius

   if diffus_limited_rate_constant <= rate_constant_standard:
      raise ValueError("The particle properties are inconsitent with the\
                        reaction rate properties k_difflim <= k_macro")

   nominator   = diffus_limited_rate_constant*rate_constant_standard
   denominator = diffus_limited_rate_constant-rate_constant_standard

   # TODO Scale to partile units system
   # Here we assume that particle simulations has units of
   # Da , mum , s
   # m^3/s -> 1e18 mum^3/s

   particle_rate_constant = nominator/denominator*1e18

   return particle_rate_constant
