import pytest

def pytest_runtest_setup(item):
    if 'notfixed' in item.keywords:
        pytest.skip("Skipping tests that are not fixed yet.")
