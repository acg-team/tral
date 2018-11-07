import os
import shutil

try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command


def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), "r") as f:
        return f.read()


# Set the home variable with user argument:
class InstallCommand(Command):
    description = "Installs TRAL"
    user_options = [
        ('home=', None, 'Directory containing the `.tral` data directory (default to user home)'),
    ]

    def initialize_options(self):
        self.home = None

    def finalize_options(self):
        if not os.path.exists(self.home):
            raise ValueError("The argument supplied in --home is not a valid path: {}".format(self.home))

    def run(self):
        datadir = os.path.join(self.home, ".tral")
        if os.path.exists(datadir):
            print("The TRAL configuration directory {} already exists. The "
                  "template configuration and datafiles are not copied to the "
                  "already existing directory at this step.".format(datadir))
        else:
            shutil.copytree("tral/tral_configuration", datadir)
            print("The TRAL configuration files and data files are now located in {}".format(datadir))


SCRIPTS1 = [os.path.join("tral", "examples", i) for i in ["example_workflow_MBE2014.py"]]
SCRIPTS2 = [os.path.join("tral", "examples", "workflow", i) for i in ["tandem_repeat_annotation_scripts.py",
                                                                      "tandem_repeat_annotation_workflow.py"]]
packages = ["tral", "tral.test", "tral.hmm", "tral.hmm.test", "tral.repeat", "tral.repeat.test",
            "tral.repeat_list", "tral.repeat_list.test", "tral.sequence", "tral.sequence.test"]

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
    packages=packages,
    scripts=SCRIPTS1 + SCRIPTS2,
    url="http://pypi.python.org/pypi/tral/",
    license="LICENSE.txt",
    description="Detect and evaluate tandem repeats in genomic sequence data.",
    long_description=read("README.rst"),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
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
            "tox >= 3.5",
        ],
        'workflow': [
            "pyfaidx==0.4.7.1",  # Used by the workflow example for fasta indexing
        ]
    },
    # package_data: None-module files, which should still be distributed are mentioned here:
    package_data={"tral": ["tral_configuration/*.ini",
                           "tral_configuration/data/hhrepid/*",
                           "tral_configuration/data/substitution_rate_matrices/*",
                           "examples/*.py", "examples/data/*",
                           "examples/workflow/*.py", "examples/workflow/*.tsv",
                           "examples/workflow/*.ini", "examples/workflow/*.hmm",
                           "examples/workflow/*.fasta",
                           "examples/workflow/split_sequence_data/*.fasta"]},
    package_dir={"tral": "tral"},
    cmdclass={
        'install': InstallCommand,
    }
)
