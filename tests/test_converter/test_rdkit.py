import pytest
import torch

from smi2xyz import Converter

from .references import refs2

smiles = [
    # Benzene
    "c1ccccc1",
    # Aspirin
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    # Ciprofloxacin
    "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
    # Vancomycin
    "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O",
]


@pytest.mark.parametrize("smile", smiles)
def test_xtb_conversion(smile: str):
    """Test conversion using the rdkit framework."""
    atol = 1e-3

    Converter.framework = "rdkit"
    xyz = Converter.smiles_to_xyz(smile)
    assert pytest.approx(refs2[smile], abs=atol) == xyz
