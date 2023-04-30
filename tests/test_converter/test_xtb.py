import pytest
import torch

from smi2xyz import Converter

from .references import refs

smiles = [
    # Benzene
    "c1ccccc1",
    # Aspirin
    # "CC(=O)OC1=CC=CC=C1C(=O)O",
]


@pytest.mark.parametrize("smile", smiles)
def test_xtb_conversion(smile: str):
    """Test conversion using the xtb framework."""
    atol = 1e1
    # NOTE: highly unstable for larger molecules (atol > 1e+1)

    Converter.framework = "xtb"
    xyz = Converter.smiles_to_xyz(smile)

    assert isinstance(xyz, torch.Tensor)
    assert pytest.approx(refs[smile], abs=atol) == xyz
