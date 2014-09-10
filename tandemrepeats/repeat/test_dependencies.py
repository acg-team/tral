import argparse
import Biopython
import configobj
import docutils
import numpy
import pytest
import scipy
import setuptools
import Sphinx
import scipy.linalg,scipy.special,scipy.stats

try:
    a = {i:2 for i in range(10)}
except:
    print("Please use Python 3 to test this code.")

print("SUCCESS! WELL DONE!")