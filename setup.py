#!/usr/bin/env python2
import os
from setuptools import setup, find_packages
    
#HERE = os.path.abspath(os.path.split(os.path.realpath(__file__))[0])

CLASSIFIERS= [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License (GPL)",
    "Natural Language :: English",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
    ]


MOD_NAME = "eggnogmapper"
VERSION = '1.0'
LONG_DESCRIPTION="""
Fast functional annotation of novel sequences using eggNOG orthology assignments.
"""

try:
    _s = setup(
        include_package_data = False,

        name = MOD_NAME,
        version = VERSION,
        packages = ['eggnogmapper'],

        # Project uses reStructuredText, so ensure that the docutils get
        # installed or upgraded on the target machine
        install_requires = [
            ],
        package_data = {

        },
        data_files = [],
        scripts=['download_eggnog_data.py', 'emapper.py'],

        # metadata for upload to PyPI
        author = "Jaime Huerta-Cepas",
        author_email = "jhcepas@gmail.com",
        maintainer = "Jaime Huerta-Cepas",
        maintainer_email = "huerta@embl.de",
        platforms = "OS Independent",
        license = "GPLv3",
        description = "Fast functional annotation of novel sequences using eggNOG orthology assignments.",
        long_description = LONG_DESCRIPTION,
        classifiers = CLASSIFIERS,
        provides = [MOD_NAME],
        keywords = "functional annotation, orthology, eggNOG",
        url = "http://eggnogdb.embl.de",
    )

except:
    print("\033[91m - Errors found! - \033[0m")
    raise
