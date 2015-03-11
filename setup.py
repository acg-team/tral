import glob
import os
import shutil
import sys

try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

# Set the home variable with user argument:
# Prettier solutions might be possible: http://stackoverflow.com/questions/677577/distutils-how-to-pass-a-user-defined-parameter-to-setup-py
try:
    i = sys.argv.index("--home")
    HOME = sys.argv[i + 1]
    del sys.argv[i+1]
    del sys.argv[i]
    if not os.path.exists(HOME):
        raise ValueError('The argument supplied in --home is not a valid path: {}'.format(HOME))
except:
    HOME=os.path.expanduser('~')


setup(
    name='tral',
    version='0.3.0',
    author='Elke Schaper',
    author_email='elke.schaper@isb-sib.ch',
    packages=['tral', 'tral.test', 'tral.hmm', 'tral.hmm.test', 'tral.repeat', 'tral.repeat.test', 'tral.repeat_list', 'tral.repeat_list.test', 'tral.sequence', 'tral.sequence.test'],
    #packages=find_packages(exclude=['tests*']),
    scripts=[os.path.join('tral', 'examples', i) for i in ['create_and_annotate_sequence_pickles.py', 'create_hmm_pickles.py', 'detect_tandem_repeats_in_sequence.py', 'example_pipeline.py', 'split_sequence_file.py']],
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
        "Biopython >= 1.64",
        "configobj >= 5.0.6",
        "docutils >= 0.11",
        "numpy >= 1.6.1",
        "pytest >= 2.5.2",
        "scipy >=0.12.0",
        "setuptools >= 5.1",
        "Sphinx >= 1.2.2",
    ],
    package_data={'tral': ['data/*.ini', 'data/hhrepid/*', 'data/substitution_rate_matrices/*']},
    package_dir={'tral': 'tral'},
)


TRAL = os.path.join(HOME, ".tral")
if os.path.exists(TRAL):
    print("Beware: The TRAL configuration directory {} already exists. The template configuration and datafiles are not copied to the already existing directory".format(TRAL))
else:
    shutil.copytree("tral/data", TRAL)
    print("The TRAL configuration files and data files are now located in {}".format(TRAL))
