import os
from ase.io import read
from ase.optimize import BFGS
from nequip.ase import NequIPCalculator


# Function to run ground-state relaxation
def run_gs(calculator, poscar_path):
    """
    Args:
        calculator (ase.calc):
        poscar_path (str):
    """
    try:
        atoms = read(poscar_path)
        atoms.calc = calculator
        dyn = BFGS(atoms, logfile=os.path.join(os.path.dirname(poscar_path), "out.log"))
        dyn.run(fmax=0.01, steps=500)

        energy = atoms.get_potential_energy()
        print(f"Energy for {poscar_path}: {energy} eV")

        contcar_path = os.path.join(os.path.dirname(poscar_path), "CONTCAR_mlp")
        atoms.write(contcar_path, format="vasp")
    except Exception as e:
        print(f"Error processing {poscar_path}: {e}")


def find_poscar_files(root_dir):
    """
    Args:
        root_dir (str):

    Returns:
        (list[str]):
    """
    poscar_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "POSCAR" in filenames and "out.log" not in filenames:
            poscar_files.append(os.path.join(dirpath, "POSCAR"))
    return poscar_files


if __name__ == "__main__":
    CALCPATH = "./deployed_model.pth"
    atom_dict = {"Si": "Si", "H": "H", "O": "O", "Al": "Al"}
    calc = NequIPCalculator.from_deployed_model(
        CALCPATH, species_to_type_name=atom_dict, device="cuda"
    )

    # Get current directory and search for POSCAR files
    og_directory = os.getcwd()
    poscarfiles = find_poscar_files(og_directory)

    # Process each POSCAR file
    for poscar in poscarfiles:
        run_gs(calc, poscar)
