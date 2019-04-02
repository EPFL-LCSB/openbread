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


from collections import OrderedDict

import numpy as np

from ..utils.constants import AVOGADRO_NUMBER
from ..utils.probabilities import calc_collision_success_probability

class ParticleModelSolution(object):
   """
   This class contains the solution of a particle simulation run
   result from C++ wrapped data
   """

   def __init__(self,results,reactions,species,volume,delta_t,n_log,is_hardsphere=True):

      # Save the raw output
      self._interface_output = results

      delta_t_log = delta_t*n_log

      species_dict = {}
      for this_species_name, this_species_dict in species.items():
         this_species_key = this_species_dict['id']
         species_dict[this_species_key] = this_species_name

      collison_dict = {}
      acceptance_dict = {}
      for this_reaction_name, reaction_dict in reactions.items():
         if len(reaction_dict['educts']) == 2:
            # Calcualte collision key

            a = reaction_dict['educts'][0]
            b = reaction_dict['educts'][1]
            if a < b:
               collison_key = (a+b)*(a+b+1)+b
            else:
               collison_key = (a+b)*(a+b+1)+a

            collison_dict[collison_key] = this_reaction_name

         if len(reaction_dict['educts']) == 1:
            acceptance_key = int((reaction_dict['educts'][0]+1)*1e6)
            counter = 0
            for product_id in reaction_dict['products']:
               acceptance_key += product_id*10**(counter*2)
               counter =+ 1

            acceptance_dict[acceptance_key] = this_reaction_name


      self.species = OrderedDict([])
      for this_key in results[0]["species"].keys():
         this_species_name = species_dict[this_key]
         self.species[this_species_name] = [result['species'][this_key] \
                                               for result in results ]

      self.mean_squared_disp = OrderedDict([])
      for this_key in results[0]["species"].keys():
         this_species_name = species_dict[this_key]
         self.mean_squared_disp[this_species_name] = [result['sq_disp'][this_key] \
                                                      / result['species'][this_key] \
                                                      if result['species'][this_key]  > 0 else 0.0 \
                                                      for result in results]

      self.collisions = OrderedDict([])
      for this_key,this_reaction_name in collison_dict.items():
         self.collisions[this_reaction_name] = []
         for result in results:
            try:
               self.collisions[this_reaction_name].append(result['collisions'][this_key])
            except IndexError:
               self.collisions[this_reaction_name].append(0.0)


      self.acceptance = OrderedDict([])
      for this_key,this_reaction_name in acceptance_dict.items():
        self.acceptance[this_reaction_name] = []
        for result in results:
           try:
             self.acceptance[this_reaction_name].append(
                result['acceptance'][this_key]/float(n_log))
           except IndexError:
             self.acceptance[this_reaction_name].append(0.0)


      self.time = [result['time'][0] for result in results ]

      # Effective macroscopic rate constans
      self.effective_rate_constants = OrderedDict([])
      self.error_effective_rate_constants = OrderedDict([])
      # Calculate the effective rate constants

      for this_reaction_name, reaction_dict in reactions.items():
         if len(reaction_dict['educts']) == 2:

            n_msrmts = round(len(self.collisions[this_reaction_name])*0.1)

            mean_collisions = np.mean(self.collisions[this_reaction_name][n_msrmts:])
            sample_size = len(self.collisions[this_reaction_name]) - n_msrmts

            first_species_id = reaction_dict['educts'][0]
            second_species_id = reaction_dict['educts'][1]

            num_first_educt = results[0]['species'][first_species_id]
            num_second_educt =  results[0]['species'][second_species_id]
            possible_collisions = float(num_first_educt*num_second_educt)
            total_possible_collisions = possible_collisions*sample_size

            first_species_name = species_dict[first_species_id]
            second_species_name = species_dict[second_species_id]
            sum_diffusion = species[first_species_name]['diff']  \
                           + species[second_species_name]['diff']
            sum_radius = species[first_species_name]['rad']  \
                         + species[second_species_name]['rad']

            micro_rate_constant = reaction_dict['rates'][0]

            avg_succes = calc_collision_success_probability(micro_rate_constant,
                                                            sum_diffusion,
                                                            sum_radius,
                                                            delta_t,
                                                            is_hardsphere)

            probablity_estimator = mean_collisions/possible_collisions

            scaling_factor = avg_succes/delta_t_log*volume*AVOGADRO_NUMBER

            this_effective_rate_constant = probablity_estimator*scaling_factor
            correction_factor = np.sqrt( probablity_estimator\
                                         *(1-probablity_estimator)
                                         /total_possible_collisions)

            # 95% confidence level for the Wald estimation of the
            # confidence interval for binomial distributions
            this_error_effective_rate_constants = 1.96*correction_factor


            # TODO Implemention using a proper unit converter
            # Scale to 1/Ms assumed input in mu^3/s
            self.effective_rate_constants[this_reaction_name] = \
                                             this_effective_rate_constant
            self.error_effective_rate_constants[this_reaction_name] = \
                                             this_error_effective_rate_constants

         if len(reaction_dict['educts']) == 1:
            micro_rate_constant = reaction_dict['rates'][0]
            n_msrmts = round(len(self.acceptance[this_reaction_name])*0.1)

            avg_succes =  np.mean(self.acceptance[this_reaction_name][n_msrmts:])
            this_effective_rate_constant = micro_rate_constant*avg_succes

            self.effective_rate_constants[this_reaction_name] = \
                                                   this_effective_rate_constant
