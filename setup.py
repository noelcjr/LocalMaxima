#!/usr/bin/env python

"""
Created on Fri Jun 24 12:10:40 2016

@author: noel

"""

# References:
# https://docs.python.org/2/distutils/setupscript.html
# https://docs.python.org/2/distutils/introduction.html#a-simple-example
# https://codeyarns.com/2014/05/30/how-to-install-and-uninstall-python-package-from-source/

import setuptools

setuptools.setup(
      name='Local Maxima',
      version='1.0.0',
      description='Python Distribution for Protein Modeling and Design',
      author='Noel Carrascal, PhD',
      author_email='noel.carrascal@entropymaxima.com',
      url='https://github.com/EntropyMaxima/LocalMaxima',
      packages=['LocalMaxima','LocalMaxima.scripts',
                'LocalMaxima.subscripts','LocalMaxima.sequence',
                'LocalMaxima.structure',
                'LocalMaxima.energy',
                'LocalMaxima.utilities','LocalMaxima.ThirdPartySoftware.charmm',
                'LocalMaxima.ThirdPartySoftware'],
      # package_data={'em': ['params/*.str','params/*.pdb','params/charmm27.ff/*']},
      license='GNU GENERAL PUBLIC LICENSE',
      long_description='A program for the classification and manipulation of protein structures.',
      platforms='Tested on a Ubuntu biolinux 16.04 LTS machine with python 3.7.3. Should work with Windows or Macs but not tested.',
      scripts=['LocalMaxima/scripts/strctr',
               'LocalMaxima/scripts/modifr',
               'LocalMaxima/scripts/prepr',
               'LocalMaxima/scripts/mmgbsa',
               'LocalMaxima/scripts/aguan',
               'LocalMaxima/scripts/ntraj',
               'LocalMaxima/scripts/nplot'],
      install_requires=[
        'biopython',
        'pandas',
      ]
)
