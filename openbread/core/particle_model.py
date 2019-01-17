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


from openbread import *
from collections import namedtuple, OrderedDict
from time import time

from .reactions import Reaction
from .solutions import ParticleModelSolution
from ..utils.tabdict import  iterable_to_tabdict, TabDict
from ..utils.constants import BOLTZMANN_CONSTANT, AVOGADRO_NUMBER
from ..utils.conversions import match_reaction_rate

class ParticleModel(object):
   """
   This class contains an OPENbread model described by elementary reactions

   """
   Crowding = namedtuple('Crowding',['volume_fraction',
                                     'mu',
                                     'sigma',
                                     'max_size'
                                     ] )

   Medium = namedtuple('Medium', ['viscosity', 'temperatur'])

   def __init__(self,medium ,crowding,volume,reactions = []):
      """

      :param medium: ParticleModel.Medium with viscosity in mum/s and temperature in K
      :param crowding: ParticleModel.Crowding with mu = log(mean kDa ), sigma in kDa sigma in mum
      :param volume: volume in L
      :param reactions: list of openbread.Reaction
      """
      self.reactions        = iterable_to_tabdict(reactions)

      self.volume           = volume

      # Assumes a quadratic box!
      size = (volume)**(1.0/3.0) # in L = dm**3
      size = size*1e5 # in mum (particle simulation unit )
      self.simulation_space = {"box_size": size}

      self.medium = medium
      self.crowding = crowding

      # defoult initial conditsions with 0
      self.initial_conditions = TabDict([])
      for this_species in self.species:
         if this_species != "CRWD":
             self.initial_conditions[this_species] = 0.0


   def add_reaction(self, reaction):
        if reaction.name in self.reactions.keys():
            raise ValueError("The reaction already exists!")
        else:
            self.reactions[reaction.name] = reaction

   @property
   def _reactions(self):
       
       return {rxn.name:{"products": [self._species[p.name]["id"] for p in rxn.products],
                         "educts": [self._species[p.name]["id"] for p in rxn.educts],
                         "rates": [match_reaction_rate(rxn.rate_constant, rxn.educts[0]+rxn.educts[1] ),] if len(rxn.educts)==2
                                  else [rxn.rate_constant,] if not rxn.constant_particle_rxn 
                                  else [0.,rxn.rate_constant] }
               for rxn in self.reactions.values()}

   @property
   def species(self):
       return iterable_to_tabdict(
           [sp for rxn in self.reactions.values()
            for sp in rxn.stoichiometry.keys()])

   @property
   def _species(self):
        dict = {e.name:{'id':i,
                       "diff":e.diffusion_constant*1e12, # Unit conversion from m^2/s to mum^2/s
                       "rad":e.collision_radius*1e6, # Unit conversion from m to mum
                       "mass":e.mass}
               for i,e in enumerate(self.species.values())}

        # Add crowding species dummy
        dict['CRWD'] = {'id': len(self.species)+1,
                        "diff": 0.0,
                        "rad": 0.0,
                        "mass": 0.0}

        return dict

   def simulate( self,
                 dt,
                 max_time,
                 log_step,
                 random_seed=1,
                 is_hardsphere=True,
                 n_sample=10,
                 is_constant_state=True,
                 is_reactive=True,
                 t_equlibriate=0):
      """
      Simulate the particle model and return reults

      """

      simulation_id =  "pymes_call_" + str(int(time()*1e24)*100+random_seed)
      py_openfpm_init(simulation_id)

      initial_conditions_scaled = {}

      for this_species,this_concetration in self.initial_conditions.items():
          this_species_number = int(round(this_concetration \
                                          *self.volume \
                                          *AVOGADRO_NUMBER))

          initial_conditions_scaled[this_species] = this_species_number

          this_species_obj = self.species[this_species]

          if is_constant_state:
              # Add reactions to keep state constant
              self.reactions[this_species+"const"] = Reaction(this_species+"const",
                                                              {this_species_obj:1},
                                                              this_species_number,
                                                              constant_particle_rxn=True)


      crowding_pyhsbrd= {  "id":self._species["CRWD"]["id"],
                           "mu":self.crowding.mu,
                           "sigma":self.crowding.sigma,
                           "volume_fraction":self.crowding.volume_fraction,
                           "viscosity":self.medium.viscosity,
                           "kBT":BOLTZMANN_CONSTANT*self.medium.temperatur,
                           "max_radius":self.crowding.max_size}


      sim = HSBRDSimulation(self._species,
                            self._reactions,
                            initial_conditions_scaled,
                            crowding_pyhsbrd,
                            self.simulation_space,
                            random_seed,
                            is_hardsphere,
                            n_sample,
                            is_constant_state)

      # Depedneted on the simulation equlibriaten
      if not is_constant_state:
         sim.equilibrate(dt,t_equlibriate)

      results = sim.simulate(dt, max_time, log_step, is_reactive)

      particle_solution = ParticleModelSolution(results,
                                                self._reactions,
                                                self._species,
                                                self.volume,
                                                dt,
                                                log_step,
                                                is_hardsphere)

      py_openfpm_finalize()

      return particle_solution
