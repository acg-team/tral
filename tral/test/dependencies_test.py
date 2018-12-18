import pytest
import sys


@pytest.mark.no_external_software_required
def test_imports():
    import argparse  # noqa: F401
    import Bio  # noqa: F401
    import configobj  # noqa: F401
    import docutils  # noqa: F401
    import numpy  # noqa: F401
    import pytest  # noqa: F401
    import scipy  # noqa: F401
    import setuptools  # noqa: F401
    import sphinx  # noqa: F401
    import scipy.linalg  # noqa: F401
    import scipy.special  # noqa: F401
    import scipy.stats  # noqa: F401


@pytest.mark.no_external_software_required
def test_python_version():

    # Check if Python3 is used
    assert sys.version_info[0] == 3
