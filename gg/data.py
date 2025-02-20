"""Common Adsorates"""

from ase import Atoms

# Dictionary of common adsorbates with their molecular structures
adsorbates = {
    "H2": Atoms(
        "H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]], cell=[10, 10, 10]  # H-H
    ),
    "O2": Atoms(
        "O2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.21]], cell=[10, 10, 10]  # O=O
    ),
    "OH": Atoms(
        "OH",
        positions=[[0.0, 0.0, 0.0], [0.15, 0.15, 0.80]],
        cell=[10, 10, 10],  # O  # H
    ),
    "OOH": Atoms(
        "OOH",
        positions=[[0.0, 0.0, 0.0], [0.46, 1.20, 0.63], [1.44, 1.06, 0.65]],
        cell=[10, 10, 10],
    ),
    "H2O": Atoms(
        "OH2",
        positions=[
            [0.0, 0.0, 0.0],  # O at origin
            [0.0, 0.757, 0.587],  # H1: positioned to create proper H-O-H angle
            [0.0, -0.757, 0.587],  # H2: symmetric position
        ],
        cell=[10, 10, 10],
    ),
    "NO": Atoms(
        "NO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.15]], cell=[10, 10, 10]  # N  # O
    ),
    "NH3": Atoms(
        "NH3",
        positions=[
            [0.0, 0.0, 0.0],  # N
            [0.0, 0.944, 0.384],  # H1
            [0.817, -0.472, 0.384],  # H2
            [-0.817, -0.472, 0.384],
        ],  # H3
        cell=[10, 10, 10],
    ),
    "NH2": Atoms(
        "NH2",
        positions=[
            [0.0, 0.0, 0.0],  # N at origin
            [0.9, 0.0, 0.15],  # H1
            [-0.9, 0.0, 0.15],  # H2
        ],
        cell=[10, 10, 10],
    ),
    "NH": Atoms(
        "NH",
        positions=[
            [0.0, 0.0, 0.0],  # N at origin
            [0.0, 0.0, 1.02],  # H at correct bond length
        ],
        cell=[10, 10, 10],
    ),
    "CO": Atoms(
        "CO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.13]], cell=[10, 10, 10]  # C-O
    ),
    "CH3": Atoms(
        "CH3",
        positions=[
            [0.0, 0.0, 0.0],  # C
            [0.0, 0.944, 0.384],  # H1
            [0.817, -0.472, 0.384],  # H2
            [-0.817, -0.472, 0.384],
        ],  # H3
        cell=[10, 10, 10],
    ),
    "CH2": Atoms(
        "CH2",
        positions=[
            [0.0, 0.0, 0.0],  # C at origin
            [0.84, 0.0, 0.1],  # H1
            [-0.84, 0.0, 0.1],  # H2
        ],
        cell=[10, 10, 10],
    ),
    "CH": Atoms(
        "CH",
        positions=[
            [0.0, 0.0, 0.0],  # C at origin
            [0.1, 0.1, 1.1],  # H at correct bond length
        ],
        cell=[10, 10, 10],
    ),
    "HCOO": Atoms(
        "OCHO",
        positions=[
            [0.0, 0.0, 0.0],  # O at center
            [0.77, 0.40, 0.90],  # C
            [0.53, 0.93, 1.82],  # H
            [2.10, 0.16, 0.84],  # O
        ],
        cell=[10, 10, 10],
    ),
}
