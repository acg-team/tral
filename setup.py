import glob
import os
from distutils.core import setup

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

HOME=os.path.expanduser('~')

setup(
    name='tral',
    version='0.3.0',
    author='Elke Schaper',
    author_email='elke.schaper@isb-sib.ch',
    packages=['tral', 'tral.hmm', 'tral.hmm.test', 'tral.repeat', 'tral.repeat.test', 'tral.repeat_list', 'tral.repeat_list.test', 'tral.sequence', 'tral.sequence.test'],
    #packages=find_packages(exclude=['tests*']),
    scripts=[os.path.join('tral', 'scripts', i) for i in ['create_and_annotate_sequence_pickles.py', 'create_hmm_pickles.py', 'detect_tandem_repeats_in_sequence.py', 'example_pipeline.py', 'split_sequence_file.py']],
    url='http://pypi.python.org/pypi/tral/',
    license='LICENSE.txt',
    description='Detect and evaluate tandem repeats in genomic sequence data.',
    long_description=read('README.markdown'),
    #include_package_data=True, # If you want files mentioned in MANIFEST.in also to be installed...
    classifiers = [
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        'Natural Language :: English',
        "Topic :: Software Development",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
        ],
    install_requires=[
        "argparse >= 1.2.1",
        "Biopython >= 1.64",
        "configobj >= 5.0.6",
        "docutils >= 0.11",
        "numpy >= 1.6.1",
        "pytest >= 2.5.2",
        "scipy >=0.12.0",
        "setuptools >= 5.1",
        "Sphinx >= 1.2.2",
    ],
    data_files=[(os.path.join(HOME, ".tral"), glob.glob("tral/data/*.ini")),
                (os.path.join(HOME, ".tral", "data", "hhrepid"), glob.glob("tral/data/hhrepid/*")),
                (os.path.join(HOME, ".tral", "data", "pValue"), []),
                (os.path.join(HOME, ".tral", "data", "substitution_rate_matrices"), glob.glob("tral/data/substitution_rate_matrices/*"))],
    package_data={'tral': ['data/*.ini', 'data/paml/*', 'data/hhrepid/*']},
    package_dir={'tral': 'tral'},
)
