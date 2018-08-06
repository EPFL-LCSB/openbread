""" OPENBREAD '

.. moduleauthor:: pyMES team


"""
import subprocess
from setuptools import setup, find_packages
from distutils.core import Extension
from distutils.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def run(self):
        subprocess.call(['make', 'clean', '-C', 'openfpm_core'])
        subprocess.call(['make', '-C', 'openfpm_core'])
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
                        'dill',
                        'salib',],

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
