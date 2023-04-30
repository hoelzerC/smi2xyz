# SMILES2XYZ

SMILES2XYZ is a Python library that allows SMILES strings to be converted into 3D structures (XYZ files). These can then be used in various applications, e.g. for molecular representations or simulations.

# Overview
```
smi2xyz/
│
├── smi2xyz/
│   ├── constants.py
│   ├── converter.py
│   └── xyz_handler.py
│
├── tests/
│   ├── test_converter.py
│   └── test_xyz_handler.py
│
├── .gitignore
├── environment.yaml
├── LICENSE
└── README.md
```

## Setup

You can create a conda environment specific for this project with the following command:

```conda env create -f environment.yml```


## Example Usage

The usage is very simple. Here is an example:

```
from smi2xyz import Converter

smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
xyz = Converter.smiles_to_xyz(smiles)
```

## LICENSE

SMILES2XYZ is licensed under the GNU GPLv3  license. For more information, see the LICENSE file.
