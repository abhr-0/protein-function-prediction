"""Sanity test verifying that the package skeleton is correctly configured."""

import src
import src.data
import src.features
import src.models
import src.utils


def test_package_imports():
    """Verify all submodules can be imported and version is defined."""
    assert src.__version__ == "0.1.0"
    assert src.data is not None
    assert src.features is not None
    assert src.models is not None
    assert src.utils is not None
