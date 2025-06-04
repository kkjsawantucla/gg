"""Common Adsorates"""

from ase import Atoms
import numpy as np

# These values are approximate and can change in systems.
# Please input values that make sense for you.
MISS = 0
max_coord_arr = np.array(
    [
        MISS,  # X
        2,  # H
        MISS,  # He
        MISS,  # Li
        MISS,  # Be
        MISS,  # B
        4,  # C
        4,  # N
        4,  # O
        4,  # F
        MISS,  # Ne
        4,  # Na
        MISS,  # Mg
        6,  # Al
        6,  # Si
        6,  # P
        6,  # S
        MISS,  # Cl
        MISS,  # Ar
        MISS,  # K
        MISS,  # Ca
        MISS,  # Sc
        12,  # Ti
        12,  # V
        12,  # Cr
        12,  # Mn
        12,  # Fe
        12,  # Co
        12,  # Ni
        12,  # Cu
        12,  # Zn
        6,  # Ga
        6,  # Ge
        6,  # As
        6,  # Se
        MISS,  # Br
        MISS,  # Kr
        MISS,  # Rb
        MISS,  # Sr
        MISS,  # Y
        MISS,  # Zr
        MISS,  # Nb
        4,  # Mo
        MISS,  # Tc
        12,  # Ru
        12,  # Rh
        12,  # Pd
        12,  # Ag
        MISS,  # Cd
        MISS,  # In
        MISS,  # Sn
        MISS,  # Sb
        MISS,  # Te
        MISS,  # I
        MISS,  # Xe
        MISS,  # Cs
        MISS,  # Ba
        MISS,  # La
        MISS,  # Ce
        MISS,  # Pr
        MISS,  # Nd
        MISS,  # Pm
        MISS,  # Sm
        MISS,  # Eu
        MISS,  # Gd
        MISS,  # Tb
        MISS,  # Dy
        MISS,  # Ho
        MISS,  # Er
        MISS,  # Tm
        MISS,  # Yb
        MISS,  # Lu
        MISS,  # Hf
        MISS,  # Ta
        MISS,  # W
        12,  # Re
        12,  # Os
        12,  # Ir
        12,  # Pt
        12,  # Au
        MISS,  # Hg
        MISS,  # Tl
        6,  # Pb
        6,  # Bi
        MISS,  # Po
        MISS,  # At
        MISS,  # Rn
        MISS,  # Fr
        MISS,  # Ra
        MISS,  # Ac
        MISS,  # Th
        MISS,  # Pa
        MISS,  # U
        MISS,  # Np
        MISS,  # Pu
        MISS,  # Am
        MISS,  # Cm
        MISS,  # Bk
        MISS,  # Cf
        MISS,  # Es
        MISS,  # Fm
        MISS,  # Md
        MISS,  # No
        MISS,  # Lr
        MISS,  # Rf
        MISS,  # Db
        MISS,  # Sg
        MISS,  # Bh
        MISS,  # Hs
        MISS,  # Mt
        MISS,  # Ds
        MISS,  # Rg
        MISS,  # Cn
        MISS,  # Nh
        MISS,  # Fl
        MISS,  # Mc
        MISS,  # Lv
        MISS,  # Ts
        MISS,  # Og
    ]
)


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
    "CO2": Atoms(
        "CO2",
        positions=[[0.0, 0.0, 0.0], [-1.1, 0.15, 0.50], [1.1, -0.15, 0.50]],
        cell=[10, 10, 10],  # C-O
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

adsorbate_names = {
    "hydrogen": adsorbates["H2"],
    "oxygen": adsorbates["O2"],
    "water": adsorbates["H2O"],
    "hydroxyl": adsorbates["OH"],
    "hydroperoxyl": adsorbates["OOH"],
    "carbon monoxide": adsorbates["CO"],
    "carbon dioxide": adsorbates["CO2"],
    "nitric oxide": adsorbates["NO"],
    "ammonia": adsorbates["NH3"],
    "amino": adsorbates["NH2"],
    "imidogen": adsorbates["NH2"],
    "methyl": adsorbates["CH3"],
    "methylene": adsorbates["CH2"],
    "methylidyne": adsorbates["CH"],
    "formate": adsorbates["HCOO"],
}
