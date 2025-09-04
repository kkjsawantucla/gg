import os
from gg.gcbh import Gcbh
from ase import Atoms
from ase.calculators.emt import EMT


def test_optimize_handles_exception(tmp_path):
    class BoomOptimizer:
        def __init__(self, atoms, logfile=None):
            self.nsteps = 0
        def run(self, fmax=None, steps=None):  # pragma: no cover - raises intentionally
            raise RuntimeError("boom")

    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = EMT()

    config = tmp_path / "config.yaml"
    config.write_text("chemical_potential:\n  H: 0\ninitialize: false\n")

    oldcwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        bh = Gcbh(atoms, config_file=str(config), optimizer=BoomOptimizer)
        _, en = bh.optimize(atoms)
    finally:
        os.chdir(oldcwd)

    assert bh.c["opt_on"] == -1
    assert en < -100000
