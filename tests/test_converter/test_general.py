import pytest
import torch

from smi2xyz import Converter


def test_unsupported_framework():
    """Test for unsupported framework."""
    Converter.framework = "unsupported_framework"
    with pytest.raises(NotImplementedError):
        Converter.smiles_to_xyz("CC(=O)OC1=CC=CC=C1C(=O)O")


def test_supported_framework():
    """Test for ssupported framework."""
    assert all(frm in Converter.options for frm in ["rdkit", "xtb"])
