from distutils.core import setup

setup(
    name='tandemrepeats',
    version='0.1.0',
    author='Elke Schaper',
    author_email='elke@inf.ethz.ch',
    packages=['hmm', 'hmm.test', 'repeat', 'repeat.test', 'repeat_list', 'repeat_list.test', 'sequence', 'sequence.test'],
    scripts=[],
    url='http://pypi.python.org/pypi/tandemrepeats/',
    license='LICENSE.txt',
    description='Detect and evaluate tandem repeats in genomic sequence data.',
    long_description=open('README.txt').read(),
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        ],
    install_requires=[
        "argparse >= 1.2.1",
        "Biopython >= 1.64",
        "docutils >= 0.11",
        "numpy >= 1.6.1",
        "pytest >= 2.5.2",
        "scipy >=0.12.0",
        "setuptools >= 5.1",
        "Sphinx >= 1.2.2",
        "XSTREAM >= 1.72",
        "T-REKS >= 1.3",
        "HHrepID >= 1.1.0",
        "TRUST >= 1.0",
    ],
    package_data={'tandemrepeats': ['data/*']}
    package_dir={'tandemrepeats': 'src/tandemrepeats'}
    scripts=['scripts/detect_tandem_repeats_in_sequence', 'scripts/example_pipeline']
)
