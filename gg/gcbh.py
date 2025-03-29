""" Main grand canonical basin hopping class """

import os
import sys
import subprocess
import shutil
import pickle
from time import strftime, localtime
import numpy as np
import yaml
from ase import Atoms
from ase import units
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.io.formats import UnknownFileTypeError
from ase.optimize.optimize import Dynamics
from ase.optimize import BFGS
from ase.neighborlist import NeighborList, natural_cutoffs
from gg.reference import get_ref_coeff
from gg.utils import (
    NoReasonableStructureFound,
    get_area,
    extract_lowest_energy_from_oszicar,
)
from gg.utils_graph import atoms_to_graph, is_unique_graph
from gg.logo import logo


__author__ = "Kaustubh Sawant, Geng Sun"


def get_current_time():
    """
    Returns:
        str: Current time
    """
    time_label = strftime("%d-%b-%Y %H:%M:%S", localtime())
    return time_label


class Gcbh(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(
        self,
        atoms: Atoms,
        logfile: str = "gcbh.log",
        trajectory: str = "gcbh.traj",
        config_file: str = None,
        restart: bool = False,
        optimizer=BFGS,
    ):
        """

        Args:
            atoms (ase.Atoms): The Atoms object to operate on.

            logfile (str, optional): If *logfile* is a string, a file will be opened.
            Defaults to "gcbh.log".

            trajectory (str, optional): Pickle file used to store the trajectory of atomic movement.
            Defaults to "gcbh.traj".

            config_file (str[.yaml file], optional): Input file to gcbh. Defaults to None.
            restart (bool, optional):

            optimizer: ase optimizer for geometric relaxation.
            Defaults to ase.optimize.BFGS
        """
        if not isinstance(atoms, Atoms):
            raise RuntimeError("The input atoms is not as Atoms object")

        if atoms.calc is None:
            raise RuntimeError("The atoms instance has no calculator")
        else:
            calc = atoms.calc

        # Intitalize by setting up the parent Dynamics Class
        super().__init__(atoms=atoms, logfile=logfile, trajectory=None)
        self.logfile.write(logo())
        self.logfile.flush()
        self.atoms.center()

        # Read Config File if it exists
        self.c = {
            "temp": 1500,
            "max_temp": None,
            "min_temp": None,
            "stop_steps": 40,
            "stop_opt": 500,
            "vasp_opt": False,
            "chemical_potential": None,
            "max_history": 25,
            "max_bond": 2,
            "max_bond_ratio": 0,
            "check_graphs": True,
            "area": False,
            "fmax": 0.05,
            "vib_correction": False,
            "initialize": True,
        }
        self.optimizer = optimizer
        if config_file:
            self.set_config(config_file)

        self.logtxt(f'Current Temp : {self.c["temp"]}')
        # Setup Temperature
        if self.c["max_temp"] is None:
            self.c["max_temp"] = 1.0 / ((1.0 / self.c["temp"]) / 1.5)
        else:
            self.c["max_temp"] = max([self.c["max_temp"], self.c["temp"]])
        self.logtxt(f'Setting Max Temp to {self.c["max_temp"]}')

        if self.c["min_temp"] is None:
            self.c["min_temp"] = 1.0 / ((1.0 / self.c["temp"]) * 1.5)
        else:
            self.c["min_temp"] = min([self.c["min_temp"], self.c["temp"]])
        self.logtxt(f'Setting Min Temp to {self.c["min_temp"]}')

        # Some file names and folders are hardcoded
        self.status_file = "current_status.pkl"
        self.opt_folder = "opt_folder"

        if os.path.exists("local_minima.traj"):
            self.lm_trajectory = Trajectory("local_minima.traj", "a", atoms)
        else:
            self.lm_trajectory = Trajectory("local_minima.traj", "w", atoms)

        if os.path.exists(trajectory):
            self.traj = Trajectory(trajectory, "a", atoms)
        else:
            self.traj = Trajectory(trajectory, "w", atoms)

        self.structure_modifiers = {}  # Setup empty class to add structure modifiers
        self.vib_correction = {}  # Setup empty class to add specific corrections
        self.c["acc_hist"] = []
        # used for adjusting the temperature of Metropolis algorithm
        # a series of 0 and 1, 0 stands for not accepted, 1 stands for accepted

        # Print the chemical potential for different elements
        if self.c["chemical_potential"]:
            self.mu = self.c["chemical_potential"]
            for k, v in self.mu.items():
                self.logtxt(f"Chemical potential of {k} is {v}")
        else:
            raise RuntimeError("No chemical potential was provided")

        self.c["energy"] = None
        self.c["fe"] = None
        self.c["fe_min"] = None
        self.c["no_impr_step"] = 0
        self.c["nsteps"] = 0
        self.c["mod_weights"] = {}

        # negative value indicates no ongoing structure optimization
        self.c["opt_on"] = -1

        # Graphing
        self.c["graphs"] = []

        if restart:
            if os.path.exists(self.status_file):
                self.update_from_file(self.status_file)
                if any(
                    arg is None
                    for arg in [self.c["energy"], self.c["fe"], self.c["fe_min"]]
                ):
                    self.initialize()

                elif self.c["opt_on"] > 0:
                    subdir = os.path.join(
                        os.getcwd(), self.opt_folder, f'opt_{self.c["opt_on"]}'
                    )
                    if os.path.isdir(subdir):
                        print(f'Restarting from {self.c["opt_on"]}')
            else:
                self.initialize()
            try:
                self.atoms = read("local_minima.traj")
                self.atoms.calc = calc
            except UnknownFileTypeError as e:
                print(f"Cannot read local_minima.traj due to {e}")
        else:
            self.initialize()

    def set_config(self, config_file: str):
        """_
        Args:
            config_file (dict): Dictionary of key inputs
        """
        with open(config_file, "r", encoding="utf-8") as f:
            input_config = yaml.safe_load(f)
        self.c.update(input_config)

    def todict(self) -> dict:
        description = {
            "type": "optimization",
            "optimizer": self.__class__.__name__,
        }
        return description

    def dump(self, filename: str):
        """dump dictionary of variables to store"""
        for key, value in self.structure_modifiers.items():
            del value.atoms
            self.c["mod_weights"][key] = value.weight
        with open(filename, "wb") as file:
            pickle.dump(self.c, file, protocol=pickle.HIGHEST_PROTOCOL)
        return

    def update_from_file(self, filename: str):
        """Update variable dictionary from file"""
        with open(filename, "rb") as file:
            c_old = pickle.load(file)
        self.c.update(c_old)
        if self.c["mod_weights"]:
            for key, value in self.c["mod_weights"].items():
                self.logtxt(f"{key} weight = {value:.2f}")
        return

    def logtxt(
        self,
        msg="",
    ):
        """Dump txt into logfile
        Args:
            msg (str, optional): The message to dump. Defaults to "".
        """
        real_message = f"{msg.strip()} \n"
        self.logfile.write(real_message)
        self.logfile.flush()
        return

    def add_modifier(self, instance, name: str):
        """
        Args:
            instance (Modifier): Instance of a ParentModifier or a child
            name (str): Name for the instance/modifier
        """
        if name in self.structure_modifiers:
            raise RuntimeError(f"Structure modifier {name} exists already!\n")
        self.structure_modifiers[name] = instance
        if name in self.c["mod_weights"]:
            self.structure_modifiers[name].weight = self.c["mod_weights"][name]
        return

    def select_modifier(self) -> str:
        """
        Returns:
            str: random modifier name
        """
        modifier_names = list(self.structure_modifiers.keys())

        modifier_weights = np.asarray(
            [self.structure_modifiers[key].weight for key in modifier_names]
        )
        modifier_weights = modifier_weights / modifier_weights.sum()
        return np.random.choice(modifier_names, p=modifier_weights)

    def add_vib_correction(self, instance, name: str):
        """
        Args:
            instance (Modifier): Instance of a ParentModifier or a child
            name (str): Name for the instance/modifier
        """
        if name in self.vib_correction:
            raise RuntimeError(f"Correction: {name} exists already!\n")
        self.vib_correction[name] = instance
        return

    def update_modifier_weights(self, name: str, action: str):
        """
        Args:
            name (str): _description_
            action (str): _description_
        """
        if name not in self.structure_modifiers:
            raise RuntimeError(f"operator name {name} not recognized")
        action = action.split()[0][0]

        if action not in ["i", "d", "r"]:
            raise RuntimeError("action must be 'increase', 'decrease' or 'rest'")

        elif action == "r":
            og_weight = self.structure_modifiers[name].og_weight
            self.structure_modifiers[name].weight = og_weight
            self.logtxt(
                f"Modifier {name} weight is reset to original weight : {og_weight:.2f}"
            )

        elif action == "i":
            og_weight = self.structure_modifiers[name].og_weight
            self.structure_modifiers[name].weight = min(
                [
                    og_weight * 2,
                    self.structure_modifiers[name].weight * 1.05,
                ]
            )

        elif action == "d":
            og_weight = self.structure_modifiers[name].og_weight
            self.structure_modifiers[name].weight = max(
                [
                    og_weight / 2.0,
                    self.structure_modifiers[name].weight / 1.05,
                ]
            )

        for key, value in self.structure_modifiers.items():
            self.logtxt(f"{key} weight = {value.weight:.2f}")
        return

    def get_ref_potential(self, atoms: Atoms):
        """
        Args:
            atoms (ase.Atoms):
        Returns:
            float: total ref value to substract
        """
        if self.mu:
            formula = atoms.get_chemical_formula()
            ref_sum = 0
            to_print = f"{formula} ="
            ref_coeff = get_ref_coeff(self.mu, formula)
            for key, value in self.mu.items():
                ref_sum += ref_coeff[key] * value
                to_print += f" {ref_coeff[key]:+.2f} {key} "
            self.logtxt(to_print)
            return ref_sum
        else:
            return 0

    def get_vib_correction(self, atoms: Atoms):
        """
        Args:
            atoms (ase.Atoms):
        Returns:
            float: total correction to add
        """
        if self.c["vib_correction"]:
            if len(self.vib_correction) > 0:
                corr_sum = 0
                for key, value in self.vib_correction.items():
                    n = value.get_n(atoms)
                    corr_value = value.weight
                    corr_sum += n * corr_value
                    self.logtxt(f"Adding correction {key}:{n}*{corr_value}")
                return corr_sum
            else:
                return 0
        else:
            return 0

    def adjust_temp(self, int_accept: int):
        """Start adjusting temperature beyond max_history
        Args:
            int_accept (int): 0 or 1
        """
        # adjust the temperatures
        self.c["acc_hist"].append(int_accept)
        if len(self.c["acc_hist"]) > self.c["max_history"]:
            self.c["acc_hist"].pop(0)
            _balance = sum(self.c["acc_hist"]) / float(self.c["max_history"])
            if _balance > 2.0 * (1 - _balance):
                self.c["temp"] = self.c["temp"] / 1.03
            elif _balance < 0.5 * (1 - _balance):
                self.c["temp"] = self.c["temp"] * 1.03

        if self.c["temp"] < self.c["min_temp"]:
            self.c["temp"] = self.c["min_temp"]
        elif self.c["temp"] > self.c["max_temp"]:
            self.c["temp"] = self.c["max_temp"]
        self.logtxt(f'Current Temp: {self.c["temp"]}')
        return

    def move(self, name: str) -> Atoms:
        """Move atoms by a random modifier."""
        atoms = self.atoms
        self.logtxt(f"Modifier '{name}' formula {atoms.get_chemical_formula()}")
        atoms = self.structure_modifiers[name].get_modified_atoms(atoms)
        atoms.center()
        return atoms

    def initialize(self):
        """Initialize Atoms"""
        if not self.c["initialize"]:
            self.logtxt(
                "Warning!!! Skipping initialization, I hope you know what you are doing"
            )
            self.c["fe"] = float("inf")
            self.c["energy"] = 0
            self.c["nsteps"] += 1
            self.append_graph(self.atoms)
            self.dump(self.status_file)
            self.traj.write(self.atoms, energy=0)
            self.lm_trajectory.write(self.atoms, energy=0, accept=1)
            return

        self.c["opt_on"] = 0
        self.atoms, en = self.optimize(self.atoms)
        self.dump(self.status_file)
        self.c["energy"] = en
        ref = self.get_ref_potential(self.atoms)
        self.c["fe"] = self.c["energy"] - ref
        if self.c["vib_correction"]:
            self.c["fe"] += self.get_vib_correction(self.atoms)
        self.c["fe_min"] = self.c["fe"]
        self.c["opt_on"] = -1
        self.dump(self.status_file)
        self.c["nsteps"] += 1
        formula = self.atoms.get_chemical_formula()
        en = self.c["energy"]
        self.append_graph(self.atoms)
        self.traj.write(self.atoms, energy=en)
        self.lm_trajectory.write(self.atoms, energy=en, accept=1)
        if self.c["area"]:
            en = en / get_area(self.atoms)
            self.c["fe"] = self.c["fe"] / get_area(self.atoms)
        self.logtxt(
            f'Atoms: {formula} E(initial): {en:.2f} F(initial) {self.c["fe"]:.2f}'
        )

    def run(self, steps: int = 4000, maximum_trial: int = 30):
        """
        Args:
            steps (int): Number of steps to run
            Defaults to 4000

            maximum_trial (int): Number of failed modification steps before terminating.
            Defaults to 30
        """
        self.logtxt("Intitial Weights:")
        for key, value in self.structure_modifiers.items():
            self.logtxt(f"{key} weight = {value.weight:.2f}")

        for step in range(steps):
            if self.c["no_impr_step"] >= self.c["stop_steps"]:
                self.logtxt(
                    f'The best solution has not improved after {self.c["no_impr_step"]} steps'
                )
                break
            self.logtxt(
                f"{step}-------------------------------------------------------"
            )

            self.logtxt(
                f'{get_current_time()}:  Starting Basin-Hopping Step {self.c["nsteps"]}'
            )

            for trials in range(maximum_trial):
                modifier_name = self.select_modifier()
                try:
                    newatoms = self.move(modifier_name)
                except (
                    NoReasonableStructureFound
                ) as emsg:  # emsg stands for error message
                    self.logtxt(
                        f"{modifier_name} did not find a good structure because {emsg} {type(emsg)}"
                    )
                else:
                    if self.append_graph(newatoms):
                        self.c["opt_on"] = self.c["nsteps"]
                        self.logtxt(
                            f"One structure found with modifier {modifier_name}"
                        )
                        self.dump(self.status_file)
                        converged_atoms, en = self.optimize(newatoms)
                        self.traj.write(converged_atoms)
                        self.append_graph(converged_atoms)
                        if self.c["opt_on"] == -1 or en < -100000:
                            self.c["opt_on"] = -1
                            self.c["nsteps"] += 1
                            continue
                        self.logtxt(f"Optimization Done with E = {en:.2f}")
                        self.accepting_new_structures(
                            converged_atoms, modifier_name, en
                        )
                        self.c["opt_on"] = -1  # switch off the optimization status
                        self.dump(self.status_file)
                        self.c["nsteps"] += 1
                        break
                    else:
                        self.logtxt("Atoms object visited previously")
                        continue
            else:
                self.logtxt(
                    f"Program does not find a good structure after {trials+1} tests"
                )
                raise RuntimeError(
                    f"Program does not find a good structure after {trials+1} tests"
                )

    def accepting_new_structures(self, newatoms: Atoms, modifier_name: str, en):
        """This function takes care of all the accepting algorithms. I.E metropolis algorithms
        newatoms is the newly optimized structure
        """
        assert newatoms is not None  # Energy_new
        fn = en - self.get_ref_potential(newatoms)  # Free_energy_new
        if self.c["vib_correction"]:
            fn += self.get_vib_correction(newatoms)

        if fn < self.c["fe"]:
            accept = True
            modifier_weight_action = "increase"

        # Check Probability for acceptance
        elif np.random.uniform() < np.exp(
            -(fn - self.c["fe"]) / self.c["temp"] / units.kB
        ):
            accept = True
            modifier_weight_action = "decrease"
        else:
            accept = False
            modifier_weight_action = "decrease"

        self.update_modifier_weights(name=modifier_name, action=modifier_weight_action)

        if accept:
            int_accept = 1
            self.logtxt(
                f'Accepted, F(old)={self.c["fe"]:.2f} F(new)={fn:.2f} at step {self.c["nsteps"]}'
            )
            self.atoms = newatoms
            self.c["energy"] = en
            self.c["fe"] = fn
            self.lm_trajectory.write(newatoms, energy=en, accept=1)

        else:
            int_accept = 0
            self.logtxt(
                f'Rejected, F(old)={self.c["fe"]:.2f} F(new)={fn:.2f} at step {self.c["nsteps"]}'
            )

        self.adjust_temp(int_accept)

        # update the best result for this basin-hopping
        if self.c["fe"] < self.c["fe_min"]:
            self.c["fe_min"] = self.c["fe"]
            self.c["no_impr_step"] = 0
        else:
            self.c["no_impr_step"] += 1

        self.dump(self.status_file)
        self.logtxt("-------------------------------------------------------")

    def optimize(self, atoms: Atoms):
        """Optimize atoms"""
        optimizer = self.optimizer
        if atoms.calc is None:
            raise RuntimeError("The atoms object has no calculator")

        self.logtxt(
            f'{get_current_time()}: Begin structure optimization routine at step {self.c["nsteps"]}'
        )
        opt_dir = self.opt_folder
        topdir = os.getcwd()
        subdir = os.path.join(topdir, opt_dir, f'opt_{self.c["nsteps"]}')
        if not os.path.isdir(subdir):
            os.makedirs(subdir)
        os.chdir(subdir)
        if self.c["vasp_opt"]:
            en = atoms.get_potential_energy()
        else:
            atoms.write("POSCAR", format="vasp")
            dyn = optimizer(atoms, trajectory="opt.traj", logfile="opt.log")
            dyn.run(fmax=self.c["fmax"], steps=self.c["stop_opt"])
            if dyn.nsteps == self.c["stop_opt"]:
                self.logtxt(
                    f'Opt is incomplete in {self.c["stop_opt"]} steps, ignore gcbh {self.c["nsteps"]} step'
                )
                self.c["opt_on"] = -1
            else:
                atoms.write("CONTCAR", format="vasp")
                self.logtxt(
                    f'{get_current_time()}: Structure optimization completed for {self.c["nsteps"]}'
                )
            en = atoms.get_potential_energy()
        os.chdir(topdir)
        return atoms, en

    def append_graph(self, atoms, unique_method = "fullgraph"):
        """Append the graph to list if its unique
        Args:
            atoms (_type_): _description_
        """
        if self.c["check_graphs"]:
            nl = NeighborList(
                natural_cutoffs(atoms), self_interaction=False, bothways=True
            )
            nl.update(atoms)
            new_g = atoms_to_graph(
                atoms,
                nl,
                max_bond=self.c["max_bond"],
                max_bond_ratio=self.c["max_bond_ratio"],
            )
            if self.c["graphs"]:
                if is_unique_graph(new_g, self.c["graphs"], comp_type = unique_method):
                    self.logtxt(
                        f"Add graph step:{self.c['nsteps']} and graph loc:{len(self.c['graphs'])}"
                    )
                    self.c["graphs"].append(new_g)
                    return True
                else:
                    return False
            else:
                self.c["graphs"].append(new_g)
                self.logtxt("Appending first graph")
        else:
            return True


class GcbhFlexOpt(Gcbh):
    """
    optimizer_file (str): Path to file that will run in opt_n folder
    copied_files (str): Files to move into opt_n folder to help in optimization
    """

    def __init__(
        self,
        atoms: Atoms,
        logfile: str = "gcbh.log",
        trajectory: str = "gcbh.traj",
        config_file: str = None,
        restart: bool = False,
        optimizer_file: str = "optimize.sh",
        copied_files: list = None,
    ):
        """

        Args:
            atoms (ase.Atoms): The Atoms object to operate on.

            logfile (str, optional): If *logfile* is a string, a file will be opened.
            Defaults to "gcbh.log".

            trajectory (str, optional): Pickle file used to store the trajectory of atomic movement.
            Defaults to "gcbh.traj".

            config_file (str[.yaml file], optional): Input file to gcbh. Defaults to None.
            restart (bool, optional):

            optimizer (str):
        """

        self.opt_file = optimizer_file
        self.copied_files = copied_files if copied_files is not None else ["opt.py"]
        super().__init__(atoms, logfile, trajectory, config_file, restart)

    def optimize(self, atoms: Atoms):
        """Optimize atoms"""

        self.logtxt(
            f'{get_current_time()}: Begin structure optimization routine at step {self.c["nsteps"]}'
        )
        script = self.opt_file
        copied_files = self.copied_files
        opt_dir = self.opt_folder
        topdir = os.getcwd()
        subdir = os.path.join(topdir, opt_dir, f'opt_{self.c["nsteps"]}')
        if not os.path.isdir(subdir):
            os.makedirs(subdir)

        if script not in copied_files:
            copied_files.append(script)
        for file in copied_files:
            assert os.path.isfile(file)
            shutil.copy(os.path.join(topdir, file), os.path.join(subdir, file))
        atoms.write(os.path.join(subdir, "input.traj"))

        try:
            os.chdir(subdir)
            opt_job = subprocess.Popen(["bash", script], cwd=subdir)
            opt_job.wait()
            if opt_job.returncode < 0:
                sys.stderr.write(
                    f"optimization does not terminate properly at {subdir}"
                )
                sys.exit(1)
        except Exception as esc:
            raise RuntimeError(
                f"some error encountered at folder {subdir} during optimizations"
            ) from esc
        else:
            fn = os.path.join(subdir, "opt.traj")
            assert os.path.isfile(fn)
            atoms = read(fn)
        finally:
            os.chdir(topdir)
        en = atoms.get_potential_energy()
        return atoms, en

    def generate_all(self, traj=None):
        """
        Generate all possible structures for all modifiers by looping over each
        modifier, generating modified atoms, filtering unique structures, and
        saving each unique structure as a POSCAR in a new subfolder with the
        required files copied.
        """
        self.logtxt(
            "Starting GcbhAll run: generating all unique structures for all modifiers."
        )
        if not traj:
            traj = [self.atoms]
        if not isinstance(traj, list):
            mod_structures = [traj]
        index = 0
        for a in traj:
            for mod_name, mod_instance in self.structure_modifiers.items():
                atoms = a.copy()
                self.logtxt(f"Processing modifier: {mod_name}")
                try:
                    # Generate modified structures.
                    mod_instance.unique = True
                    mod_instance.print_movie = True
                    mod_structures = mod_instance.get_modified_atoms(atoms)
                except (
                    NoReasonableStructureFound
                ) as emsg:  # emsg stands for error message
                    self.logtxt(
                        f"{mod_name} did not find a good structure because {emsg} {type(emsg)}"
                    )
                    continue

                if not isinstance(mod_structures, list):
                    mod_structures = [mod_structures]

                # Process each generated structure.
                for struct in mod_structures:
                    # Check uniqueness by using the already available append_graph method.
                    if self.c["check_graphs"]:
                        # append_graph returns True if the structure is unique.
                        is_unique = self.append_graph(struct,unique_method=mod_instance.unique_method)
                    else:
                        is_unique = True

                    if is_unique:
                        # Increment the step counter to generate a unique folder name.
                        subfolder = os.path.join(
                            os.getcwd(),
                            self.opt_folder,
                            f"opt_{self.c['nsteps']}",
                            f"opt_{index}",
                        )

                        if not os.path.isdir(subfolder):
                            os.makedirs(subfolder)

                        # Write the POSCAR file in VASP format.
                        poscar_file = os.path.join(subfolder, "POSCAR")
                        struct.write(poscar_file, format="vasp")

                        # Copy any additional required files into the new subfolder.
                        for file in self.copied_files:
                            if os.path.isfile(file):
                                shutil.copy(file, subfolder)
                            else:
                                self.logtxt(f"File {file} not found for copying.")

                        self.logtxt(
                            f"Unique structure from modifier '{mod_name}' saved in {subfolder}"
                        )
                    else:
                        self.logtxt(
                            f"Duplicate structure from modifier '{mod_name}' encountered; skipping."
                        )
                    index += 1
        self.c["nsteps"] += 1

    def update_lowest_energy(self):
        """
        Loops over all subdirectories in self.opt_folder.
        """

        best_fe = self.c["fe"]
        best_atoms = self.atoms

        opt_folder_path = os.path.join(os.getcwd(), self.opt_folder)
        if not os.path.isdir(opt_folder_path):
            self.logtxt(f"Optimization folder '{self.opt_folder}' does not exist.")
            return

        # Loop over every subdirectory in the optimization folder.
        for subdir in os.listdir(opt_folder_path):
            subdir_path = os.path.join(opt_folder_path, subdir)
            if os.path.isdir(subdir_path):
                contcar_path = os.path.join(subdir_path, "CONTCAR")
                oszicar_path = os.path.join(subdir_path, "OSZICAR")
                if os.path.isfile(contcar_path) and os.path.isfile(oszicar_path):
                    en = extract_lowest_energy_from_oszicar(oszicar_path)
                    if en is not None:
                        self.logtxt(f"Found energy {en:.2f} in {subdir_path}")
                        atoms = read(contcar_path, format="vasp")
                        fn = en - self.get_ref_potential(atoms)
                        if self.c["vib_correction"]:
                            fn += self.get_vib_correction(atoms)
                        self.logtxt(f"F at {subdir_path} is {fn}")
                        if fn < best_fe:
                            best_fe = fn
                            best_atoms = atoms
                            self.logtxt(f"Accepted; F(new)={fn} at {subdir_path}")
                        else:
                            self.logtxt(f"Rejected; F(new)={fn} at {subdir_path}")
        self.atoms = best_atoms
        self.c["fe"] = best_fe
