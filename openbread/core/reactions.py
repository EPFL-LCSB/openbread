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

class Reaction:
    def __init__(self,name,stoichiometry, rate_constant, constant_particle_rxn = False):
        """
        Reaction object for bimolecular and monomolecular reactions

        Warning: Reactions can not have more than 2 educts and products!

        :param name: Name of the reaction
        :param stoichiometry: Dict of {openbrad.Species species : int stoichiometry}
        :param rate_constant: Double rate mass action rate constant in 1/Ms or 1/s

        """

        self.name = name

        self.stoichiometry = stoichiometry

        self.rate_constant = rate_constant

        self.constant_particle_rxn = constant_particle_rxn

    @property
    def educts(self):
        return [sp for sp, stoich in self.stoichiometry.items()
                          for i in range(abs(stoich)) if stoich < 0 ]

    @property
    def products(self):
        return [sp for sp, stoich in self.stoichiometry.items()
                          for i in range(abs(stoich)) if stoich > 0 ]

    def __repr__(self):
        reactants = ""
        products = ""

        # ToDo nice string repr

        return "{} --> {}".format(reactants,products)