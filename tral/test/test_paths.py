import os
import tempfile
import pytest
from tral import paths
from urllib.request import URLError


@pytest.fixture
def config_dir():
    with tempfile.TemporaryDirectory(prefix="tral-") as config_dir:
        yield config_dir

    # config_dir automatically cleaned up after yield


DISABLED_URL = "http://0.0.0.0/disabled/"  # sentinel non-working URL


def test_package_config(config_dir):
    # Should get copied from package
    ini_file = paths.config_file('config.ini', config_dir=config_dir, config_url=DISABLED_URL)

    assert ini_file == os.path.join(config_dir, 'config.ini')
    assert os.path.exists(ini_file)

    # Not found in package and can't be downloaded
    with pytest.raises(URLError):
        ini_file = paths.config_file('foo', config_dir=config_dir, config_url=DISABLED_URL)


def test_download_config(config_dir):
    small_file = "pvalue/AA/PhyloScore_AA_Empirical_MolecularClock_SameRateOnBranches_10000.descriptor"

    # not in package
    with pytest.raises(URLError):
        downloaded = paths.config_file(small_file, config_dir=config_dir, config_url=DISABLED_URL)

    # can be downloaded
    downloaded = paths.config_file(small_file, config_dir=config_dir)

    assert downloaded == os.path.join(config_dir, small_file)
    assert os.path.exists(downloaded)

    # Check file contents
    expected = """Phyloscore on random data - AA based on 3mers
sequenceType: AA
nTotalRuns: 10000
batchID: PhyloScore_AA_Empirical_MolecularClock_SameRateOnBranches_10000
CODEML: Empirical, molecular clock assumption, star tree, constant branch lengths
SamplesPerRun: 10000"""

    with open(downloaded, 'r') as f:
        assert f.read() == expected
