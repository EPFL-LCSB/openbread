# -*- coding: utf-8 -*-
"""
.. module:: pymes
   :platform: Unix, Windows
   :synopsis: Python Multiscale Enzyme Simulator

.. moduleauthor:: pyMES team

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
from collections import namedtuple

class Species(object):
   """
   Class to defines the microscopic species in the particle simulation
   :return:
   """

   def __init__(self,name,diffusion_constant,collision_radius,mass):
      self.name               = name
      self.diffusion_constant = diffusion_constant
      self.collision_radius   = collision_radius
      self.mass               = mass


   def __add__(self,other_species):
      summed_name               = self.name + other_species.name
      summed_diffusion_constant = self.diffusion_constant + \
                                  other_species.diffusion_constant
      summed_collision_radius   = self.collision_radius + \
                                  other_species.collision_radius
      summed_mass               = self.mass + other_species.mass

      sum_of_species = Species(summed_name,
                               summed_diffusion_constant,
                               summed_collision_radius,
                               summed_mass)

      return sum_of_species
