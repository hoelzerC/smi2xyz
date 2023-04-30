import torch

Tensor = torch.tensor

from .constants import ATOMIC_NUMBER, PSE, AA2AU


class XYZ_Handler:
    """
        A class for reading and writing xyz files.
        For the xyt format specification see https://en.wikipedia.org/wiki/XYZ_file_format

        NOTE: All output is given in a.u.
    """

    def read_xyz(fp: str) -> tuple[Tensor, Tensor]:
        """Reads an xyz file and returns a tuple of two tensors.

        Parameters
        ----------
        fp : str
            The file path of the xyz file to be read.

        Returns
        -------
        tuple[Tensor, Tensor]
            A tuple containing two tensors. The first tensor is a
            2D tensor of shape (num_atoms, 4) where each row contains the atomic number,
            x, y, and z coordinates of an atom in the molecule. The second tensor is a
            1D tensor of length 1 containing the total charge of the molecule. If no charge
            is specified in the xyz file, the charge tensor will contain a value of 0.
        """
        data, charge = [], 0

        with open(fp, "r") as file:
            for line_number, line in enumerate(file):
                if line_number == 0:
                    num_atoms = int(line)
                elif line_number == 1:
                    if "charge=" in line:
                        charge = int(line.split("=")[1])
                else:
                    ele, x, y, z = line.split()
                    data.append(
                        [
                            ATOMIC_NUMBER.get(ele),
                            float(x) * AA2AU,
                            float(y) * AA2AU,
                            float(z) * AA2AU,
                        ]
                    )
        return torch.tensor(data), torch.tensor([charge])

    def write_xyz(fp: str, data: Tensor, comment: str = None) -> None:
        """Writes a Tensor to an xyz file.

        Parameters
        ----------
        fp : str
            The file path to write the data to.
        data : Tensor
            The input Tensor of shape (N, 4), where N is the number of atoms in the molecule.
            The columns of the Tensor should be in the order: atomic number, x coordinate, y 
            coordinate, and z coordinate.
        comment : str, optional
            A comment to include in the header of the xyz file, by default None
        """

        with open(fp, "w") as f:
            # write the number of atoms as the first line
            f.write(f"{data.shape[0]}\n")

            # write a comment line as the second line
            f.write(f"{comment}\n")

            # write the atomic symbols and positions
            for d in data:
                f.write(
                    "{} {:.6f} {:.6f} {:.6f}\n".format(
                        PSE[d[0].item()], d[1] / AA2AU, d[2] / AA2AU, d[3] / AA2AU
                    )
                )