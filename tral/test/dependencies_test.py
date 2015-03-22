import pytest
import sys

@pytest.mark.no_external_software_required
def test_imports():
    import argparse
    import Bio
    import configobj
    import docutils
    import numpy
    import pytest
    import scipy
    import setuptools
    import sphinx
    import scipy.linalg,scipy.special,scipy.stats

@pytest.mark.no_external_software_required
def test_python_version():

    # Check if Python3 is used
    assert sys.version_info[0] == 3

