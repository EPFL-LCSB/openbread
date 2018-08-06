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


#Example of how the python interface of HSBRD works
import openbread
import numpy as np


species = { "S":{"id":0,"diff":300.0,"rad":1e-3,"mass":0.2},
            "P":{"id":1,"diff":300.0,"rad":1e-3,"mass":0.2},
            "E":{"id":2,"diff":100.0,"rad":3e-3,"mass":60.0},
            "ES":{"id":3,"diff":100.0,"rad":3e-3,"mass":60.2}, }

k1f = 1.01685795e-04;
k1b = 4.49346405e+03;
k2f = 4.49346405e+03;
k2b = 1.22023123e-04;

reactions = { "E+S->ES": {"educts":[0,2], "products":[3], "rates":[k1f,]},
              "ES->E+S": {"educts":[3], "products":[0,2], "rates":[k1b,]},
              "P+E->ES": {"educts":[1,2], "products":[3], "rates":[k2f,]},
              "ES->P+E": {"educts":[3], "products":[1,2], "rates":[k2b,]},
             }

initial_conditions = {"S":1350,"P":934,"E":120,"ES":2,}

crowding = {"id":13.0,
            "mu":np.log(31.9),
            "sigma":0.825,
            "volume_fraction":0.1,
            "viscosity":1e-3,
            "kBT":1.38e-23*298.15,
            "max_radius":10e-3}

simulation_space = {"box_size":0.1}

openbread.py_openfpm_init()

sim = openbread.HSBRDSimulation(species,
                                reactions,
                                initial_conditions,
                                crowding,
                                simulation_space,
                                1,
                                True,
                                100)

results = sim.simulate(1e-9,1e-8,1)
