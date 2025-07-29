import os
import yaml
from ase import Atoms
from gg.phase_diag import get_entries_from_folders


class DummyModifier:
    def __init__(self, weight=0.5):
        self.weight = weight

    def get_ind_to_remove_list(self, atoms):
        return [[0, 1]]


def test_get_entries_with_vib(tmp_path):
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
    contcar = tmp_path / 'CONTCAR'
    oszicar = tmp_path / 'OSZICAR'
    atoms.write(contcar, format='vasp')
    with open(oszicar, 'w') as f:
        f.write('E0= -10.0\n')

    mu = {'chemical_potential': {'H': 0.0, 'O': 0.0}}
    mu_file = tmp_path / 'input.yaml'
    with open(mu_file, 'w') as f:
        yaml.safe_dump(mu, f)

    modifier = DummyModifier(weight=0.5)
    entries = get_entries_from_folders(
        'H',
        'O',
        base_folders=[str(tmp_path)],
        mu_path=str(mu_file),
        vib_corrections={'dummy': modifier},
    )
    assert len(entries) == 1
    area = 100.0
    assert abs(entries[0].energy - (-9.5 / area)) < 1e-6
    assert abs(entries[0].n1 - (2.0 / area)) < 1e-6
    assert abs(entries[0].n2) < 1e-8
