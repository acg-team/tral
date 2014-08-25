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
        "numpy >= 1.6.1",
    ],
    package_data={'tandemrepeats': ['data/*']}
    package_dir={'tandemrepeats': 'src/tandemrepeats'}
    scripts=['scripts/detect_tandem_repeats_in_sequence', 'scripts/example_pipeline']
)
