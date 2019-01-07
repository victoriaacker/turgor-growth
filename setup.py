# -*- coding: latin-1 -*-
"""

    setup
    ~~~~~
    
    Setup script for installation.
    
    See README.md for installing procedure.

    :license: CeCILL-C, see LICENSE for details.
    
"""

import ez_setup
import pkg_resources
import sys
from setuptools import setup, find_packages

import turgorgrowth

ez_setup.use_setuptools()

if sys.version_info < (2, 7):
    print('ERROR: Turgor-Growth requires at least Python 2.7 to run.')
    sys.exit(1)

if sys.version_info >= (3, 0):
    print('WARNING: Turgor-Growth has not been tested with Python 3.')

pkg_resources.require('numpy>=1.15.1', 'pandas>=0.23.4', 'scipy>=1.1.0')

setup(
    name="Turgor-Growth",
    version=turgorgrowth.__version__,
    packages=find_packages(),
    include_package_data=True,
    author="R.Barillot and T.de Swaef",
    author_email="romain.barillot@inra.fr, tom.deswaef@ilvo.vlaanderen.be",
    description="A turgor-driven model of leaf growth ",
    long_description="A turgor-driven model of leaf growth linking water and C dynamics. This project is adapted from Coussement et al. (2018).",
    license="CeCILL-C",
    keywords="functional-structural plant model, wheat, ode, system integration, scipy",
    url="https://sourcesup.renater.fr/turgor-growth/",
    download_url="",
)
