#!/usr/bin/env python3
import os, sys
from setuptools import setup, find_packages

if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported')
    
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
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
    ]


MOD_NAME = "eggnog-mapper"
VERSION = '2.0.5'
LONG_DESCRIPTION="""
Fast functional annotation of novel sequences using eggNOG orthology assignments.
"""

try:
    _s = setup(
        include_package_data = False,

        name = MOD_NAME,
        version = VERSION,
        packages = find_packages(),

        # Project uses reStructuredText, so ensure that the docutils get
        # installed or upgraded on the target machine
        install_requires = [
            ],
        package_data = {

        },
        data_files = [],
        scripts=['download_eggnog_data.py', 'emapper.py', 'hmm_mapper.py', 'hmm_server.py', 'hmm_worker.py'],

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
        keywords = "functional annotation, orthology, eggNOG",
        url = "http://eggnogdb.embl.de",
        python_requires='>=3.7',
    )

except:
    print("\033[91m - Errors found! - \033[0m")
    raise
