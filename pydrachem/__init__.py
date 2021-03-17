"""
PydraChem
Python package to streamline figure generation for systems
"""

# Add imports here
from .pydrachem import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
