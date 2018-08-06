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

import numpy as np


from scipy.special import erfc
from numpy import pi,sqrt,exp

class post_proccesing:

	def __init__(self,_rates, _diffusion , _radius ,_volume , _delta_t):
		self.rates = _rates
		self.diffusion = _diffusion
		self.radius = _radius
		self.delta_t = _delta_t
		self.volume = _volume

	def calculate_rates(self, _folder, _is_volume_exlusion,):

		is_volume_exlusion = _is_volume_exlusion;

		self.collision_data = np.loadtxt(_folder+"/collisions_0", skiprows =2 )
		self.species_data = np.loadtxt(_folder+"/species_number_0", skiprows =2 )
		self.acceptance_data = np.loadtxt(_folder+"/acceptance_0", skiprows =2 )


		D = self.diffusion[0] + self.diffusion[1]
		R = self.radius[0] + self.radius[1]

		if _is_volume_exlusion:
			norm = norm_bimol_rate(D,R,self.delta_t)
		else:
			norm = R**3/3.0


		# Maximum likly hood estimator for bimolecular rates
		r_1f = 	np.mean(self.collision_data[:,3]/self.species_data[:,1]/self.species_data[:,3])/self.delta_t*self.volume*6e23*(1.0-exp(-self.delta_t*self.rates[0]/(4.0*pi*norm)))
		r_2b = 	np.mean(self.collision_data[:,5]/self.species_data[:,1]/self.species_data[:,4])/self.delta_t*self.volume*6e23*(1.0-exp(-self.delta_t*self.rates[3]/(4.0*pi*norm)))

		# Monomolceular rates


		if _is_volume_exlusion:
			r_1b = np.mean(self.acceptance_data[:,1]/self.species_data[:,2])*self.rates[1]
			r_2f = np.mean(self.acceptance_data[:,2]/self.species_data[:,2])*self.rates[2]
		else:
			r_1b = self.rates[1]
			r_2f = self.rates[2]

		return r_1f, r_1b , r_2f , r_2b



# Normalization factor for Brownian dyanamics simulation
def norm_bimol_rate(diffusion,dist,delta_t):
		#Bi mol rxn scaling
		sig_2 =4.0*delta_t*diffusion
		sig = sqrt(sig_2)

		exp_4_r_sig = exp(-4.0*dist**2/sig_2)

		#Expresion
		A = (sig**3 - 2.0*dist**2*sig)*exp_4_r_sig
		B = 6.0*dist**2*sig - sig**3 + 4.0*sqrt(pi)*dist**3*erfc(2.0*dist/sig)

		return (A + B )/12.0/sqrt(pi)
