"""
Unit and regression test for the pydrachem package.
"""

# Import package, test suite, and other packages as needed
import pydrachem
import pytest
import sys

def test_pydrachem_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pydrachem" in sys.modules
