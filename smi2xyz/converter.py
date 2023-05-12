import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import tempfile
import torch

from .constants import AA2AU
from .xyz_handler import XYZ_Handler

Tensor = torch.tensor


class Converter:
    """Simple converter capable of creating 3D `.xyz` files from SMILES `.smi` files."""

    options = ["rdkit", "xtb"]
    """Supported frameworks."""

    framework = "rdkit"
    """Default framework used for conversion."""

    xtb_path = "xtb"
    """Path to installed `xtb` binary. Only required for xtb framework."""

    obabel_path = "obabel"
    """Path to installed `obabel` binary. Only required for xtb framework."""

    @classmethod
    def smiles_to_xyz(cls, smiles: str):

        if cls.framework == "rdkit":
            return cls._rdkit(smiles)
        elif cls.framework == "xtb":
            return cls._xtb(smiles)
        else:
            raise NotImplementedError(
                "Unsupported framework for SMILES to xyz conversion."
            )

    @classmethod
    def _rdkit(cls, smiles: str) -> Tensor:
        """Converts a SMILES string to a 3D structure using RDKit.

        Parameters
        ----------
        smiles : str
            A SMILES string to be converted.

        Returns
        -------
        Tensor
            The molecule information in a 2D tensor of shape (num_atoms, 4) 
            where each row contains the atomic number, x, y, and z coordinates
            of an atom in the molecule.
        """

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        # generate molecule object with 3D coordinates
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        return torch.tensor(
            [
                (
                    mol.GetAtomWithIdx(i).GetAtomicNum(),
                    mol.GetConformer().GetAtomPosition(i).x * AA2AU,
                    mol.GetConformer().GetAtomPosition(i).y * AA2AU,
                    mol.GetConformer().GetAtomPosition(i).z * AA2AU,
                )
                for i in range(mol.GetNumAtoms())
            ]
        )

    @classmethod
    def _xtb(cls, smiles: str) -> Tensor:
        """Converts a SMILES string to a 3D structure using xtb (and obabel).
           Essentially, obabel is used to convert SMILES to 2D and xtb to 
           convert to 3D.

        Parameters
        ----------
        smiles : str
            A SMILES string to be converted.

        Returns
        -------
        Tensor
            The molecule information in a 2D tensor of shape (num_atoms, 4) 
            where each row contains the atomic number, x, y, and z coordinates
            of an atom in the molecule.

        Raises
        ------
        IOError
            Raise error if calculation failed.
        """

        def run(cmd: str, cwd=None):
            """Helper function."""
            subprocess.run(
                cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd
            )

        # NOTE: `xtb` and `obabel` required on the system
        run(f"{cls.xtb_path} -h")
        run(f"{cls.obabel_path} -h")  # will raise a FileNotFoundError

        with tempfile.TemporaryDirectory() as tmpdir:

            outSMI = tmpdir + "/smiles.smi"
            out2D = tmpdir + "/2D.sdf"
            out3D = tmpdir + "/3D.xyz"

            with open(outSMI, "w",encoding="UTF-8") as file:
                file.write(smiles)

            # convert SMILES to 2D using obabel
            obabel_cmd = f"{cls.obabel_path} {outSMI} -O {out2D} --gen2d -h"
            run(obabel_cmd, cwd=tmpdir)

            # convert to 3D using xtb
            xtb_cmd = f"{cls.xtb_path} {out2D} --gfn 2 --sp"
            run(xtb_cmd, cwd=tmpdir)

            # for gfn2 usage
            if os.path.isfile(os.path.join(tmpdir, "xtbtopo.sdf")):
                file_opt = "xtbtopo.sdf"
            # for gfnff usage
            elif os.path.isfile(os.path.join(tmpdir, "gfnff_convert.sdf")):
                file_opt = "gfnff_convert.sdf"
            else:
                raise IOError(f"No file found for {outSMI}. Exit.")

            # conduct geometry optimisation on 3D-structures
            opt_cmd = f"{cls.xtb_path} {file_opt} --gfn2 --opt"
            run(opt_cmd, cwd=tmpdir)

            # convert 3D-structures to .xyz files
            obabel_cmd = f"{cls.obabel_path} xtbopt.sdf -O {out3D}"
            run(obabel_cmd, cwd=tmpdir)

            # check for convergence
            if not os.path.isfile(os.path.join(tmpdir, ".xtboptok")):
                raise IOError(f"Optimization for {outSMI} failed. Exit.")

            # read output
            xyz, chrg = XYZ_Handler.read_xyz(self=None, fp=out3D)
        return xyz
