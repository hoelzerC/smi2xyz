"""
Run tests for IO of `xyz` files.
"""
from pathlib import Path
import pytest
import tempfile

from smi2xyz import XYZ_Handler

from .references import refs

sample_path = Path(__file__).resolve().parents[0] / "samples"
samples = ["h2o.xyz", "nh3.xyz", "vancoh2.xyz"]


@pytest.mark.parametrize("sample", samples)
def test_read(sample: str) -> None:
    """Test reading of xyz files."""
    atol = 1e-4

    data, chrg = XYZ_Handler.read_xyz(sample_path / sample)

    assert pytest.approx(refs[sample]["data"], abs=atol) == data
    assert pytest.approx(refs[sample]["chrg"], abs=0.0) == chrg


@pytest.mark.parametrize("sample", samples)
def test_write(sample: str) -> None:
    """Test writing of xyz files via cyclic read/write."""
    atol = 1e-4

    data, chrg = XYZ_Handler.read_xyz(sample_path / sample)

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        XYZ_Handler.write_xyz(temp_file.name, data)
        written_data, written_chrg = XYZ_Handler.read_xyz(temp_file.name)

    assert pytest.approx(data, abs=atol) == written_data
    assert pytest.approx(chrg, abs=0.0) == written_chrg
