import os
from setuptools import setup, find_packages


def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), "r") as f:
        return f.read()


this_directory = os.path.abspath(os.path.dirname(__file__))
long_description = read(this_directory, "README.md")


SCRIPTS1 = [os.path.join("tral", "examples", i) for i in ["example_workflow_MBE2014.py"]]
SCRIPTS2 = [os.path.join("tral", "examples", "workflow", i) for i in ["tandem_repeat_annotation_scripts.py",
                                                                      "tandem_repeat_annotation_workflow.py"]]

# Load the version number from tral/__init__.py
__version__ = "Undefined"
for line in open('tral/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(
    name="tral",
    version=__version__,
    author="Elke Schaper",
    author_email="elke.schaper@isb-sib.ch",
    packages=find_packages(),
    scripts=SCRIPTS1 + SCRIPTS2,
    url="http://pypi.python.org/pypi/tral/",
    license="LICENSE.txt",
    description="Detect and evaluate tandem repeats in genomic sequence data.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Natural Language :: English",
        "Topic :: Software Development",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython >= 1.64",
        "configobj >= 5.0.6",
        "numpy >= 1.6.1",
        "scipy >=0.12.0",
    ],
    setup_requires=["pytest-runner"],
    tests_require=[
        "pytest >= 2.5.2",
    ],
    # Install with e.g. `python setup.py install tral[docs]`
    extras_require={
        'docs': [
            "docutils >= 0.11",
            "pypandoc >= 0.9.6",
            "Sphinx >= 1.2.2",
        ],
        'develop': [
            "flake8 >= 3.6",
            "flake8-colors",
            "tox >= 3.5",
            "pytest >= 2.5.2",
        ],
        'workflow': [
            "pyfaidx==0.4.7.1",  # Used by the workflow example for fasta indexing
        ]
    },
    # package_data: None-module files, which should still be distributed are mentioned here:
    package_data={"tral": ["tral_configuration/*.ini",
                           "tral_configuration/data/*",
                           "examples/*.py", "examples/data/*",
                           "examples/workflow/*.py", "examples/workflow/*.tsv",
                           "examples/workflow/*.ini", "examples/workflow/*.hmm",
                           "examples/workflow/*.fasta",
                           "examples/workflow/split_sequence_data/*.fasta"]},
    package_dir={"tral": "tral"},
)
