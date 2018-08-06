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
from time import time

from skimpy.utils import iterable_to_tabdict

from .solutions import ParticleModelSolution
from ..utils.constants import BOLTZMANN_CONSTANT, AVOGADRO_NUMBER


class ParticleModel(object):
   """
   This class contains a multiscale particle model descibed by the reactions,
   species participating in the reactions
   """
   Crowding = namedtuple('Crowding',['volume_fraction',
                                     'mu',
                                     'sigma',
                                     'max_size'
                                     ] )

   Medium = namedtuple('Medium', ['viscosity', 'temperatur'])

   def __init__(self,name,reactions,species,medium,crowding,volume):
      self.name             = name
      self.reactions        = reactions
      self.species          = species
      self.volume           = volume

      # Assumes a quadratic box!
      # TODO PROPER SCALING  WE ASSUME THAT THIS INPUT VOLUME IS IN L
      #
      size = (volume)**(1.0/3.0) # in L = dm**3
      size = size*1e5 # in mum (particle simulation unit )

      self.simulation_space = {"box_size": size}

      self.medium = medium
      self.crowding = crowding

      # defoult initial conditsions with 0
      self.initial_conditions = iterable_to_tabdict([])

      for this_species in self.species.keys():
         if this_species != "CRWD":
             self.initial_conditions[this_species] = 0.0



   def simulate( self,
                 dt,
                 max_time,
                 log_step,
                 random_seed = 1,
                 is_hardsphere = True,
                 n_sample = 10,
                 is_constant_state = True,
                 t_equlibriate = 0):
      """ Simulate the particle model and return reults """

      simulation_id =  "pymes_call_" + str(int(time()*1e24))
      py_openfpm_init(simulation_id)


      initial_conditions_scaled = {}

      for this_species,this_concetration in self.initial_conditions.items():
         this_species_number = int(round(this_concetration \
                                       *self.volume    \
                                       *AVOGADRO_NUMBER))

         initial_conditions_scaled[this_species] = this_species_number

         this_species_id = self.species[this_species]["id"]

         if is_constant_state:
            # Add reactions to keep state constant
            self.reactions[this_species+"const"] = {"educts":[],
                                                    "products":[this_species_id],
                                                    "rates":[0,this_species_number]}


      crowding_pyhsbrd= {  "id":self.species["CRWD"]["id"],
                           "mu":self.crowding.mu,
                           "sigma":self.crowding.sigma,
                           "volume_fraction":self.crowding.volume_fraction,
                           "viscosity":self.medium.viscosity,
                           "kBT":BOLTZMANN_CONSTANT*self.medium.temperatur,
                           "max_radius":self.crowding.max_size}


      sim = HSBRDSimulation(self.species,
                            self.reactions,
                            initial_conditions_scaled,
                            crowding_pyhsbrd,
                            self.simulation_space,
                            random_seed,
                            is_hardsphere,
                            n_sample)

      # Depedneted on the simulation equlibriaten
      if not is_constant_state:
         sim.equilibrate(dt,t_equlibriate)

      results = sim.simulate(dt,max_time,log_step)

      particle_solution = ParticleModelSolution(results,
                                                self.reactions,
                                                self.species,
                                                self.volume,
                                                dt,
                                                log_step,
                                                is_hardsphere)
      # Todo Wrap results in a readable format



      return particle_solution
