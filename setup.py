# -*- coding: utf-8 -*-
"""
.. module:: openbread
   :platform: Unix, Windows
   :synopsis: OPENFPM based brownian reaction dynamics

.. moduleauthor:: openbread team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
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

import subprocess
from setuptools import setup, find_packages
from distutils.core import Extension
from distutils.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def run(self):

        subprocess.call(['make', 'clean', '-C', 'openfpm_core'],)
        subprocess.call(['ls','openfpm_core'], )

        try:
            subprocess.check_output('make -C openfpm_core',
                                    stderr=subprocess.STDOUT,
                                    shell=True)
        except subprocess.CalledProcessError as inst:

            raise RuntimeError("CMD: {} failed with error {}".format(inst.cmd,inst.output))

        _build_ext.run(self)


version_tag = '0.0.1'


setup(name='openbread',
      version=version_tag,
      author='OPENBREAD Team ',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/OPENBREAD/',
      download_url='https://github.com/EPFL-LCSB/OPENBREAD/archive/'+version_tag+'.tar.gz',
      install_requires=['sympy >= 1.1.',
                        'pytest',
                        'scipy',
                        'numpy',
                        'dill'],

      packages = find_packages(),
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='OPRNFPM based Brownian Reaction Dyncamics',
      keywords=['OPENFPM',
                'Brownian Reaction Dyncamics',
                'agent-based','kinetic','models'],

      license='Apache2',

      # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry'
            'Environment :: Console',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: Apache Software License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],

      cmdclass={'build_ext': build_ext},
     )
