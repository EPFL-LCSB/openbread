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


from scipy.special import erfc
from numpy import pi,sqrt,exp


def calc_collision_success_probability(micro_rate_constant,
                                       sum_diffusion,
                                       sum_radius,
                                       delta_t,
                                       is_hardsphere = True):
   if is_hardsphere:
      v_eff = calc_effective_volume(sum_diffusion,sum_radius,delta_t)
   else:
      v_eff = 4.0/3.0*pi*sum_radius**3

   avg_succes_probability = 1.0-exp(-micro_rate_constant/v_eff*delta_t)

   return avg_succes_probability






def calc_effective_volume(diffusion,dist,delta_t):
   """ Normalization factor for Brownian dyanamics simulation """

   #Bi mol rxn scaling
   sig_2 =4.0*delta_t*diffusion
   sig = sqrt(sig_2)

   exp_4_r_sig = exp(-4.0*dist**2/sig_2)

   #Expresion
   A = (sig**3 - 2.0*dist**2*sig)*exp_4_r_sig
   B = 6.0*dist**2*sig - sig**3 + 4.0*sqrt(pi)*dist**3*erfc(2.0*dist/sig)

   effective_volume = 4.0*pi*(A + B)/12.0/sqrt(pi)
   
   return effective_volume
