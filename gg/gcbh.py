""" Main grand canonical basin hopping class """

import os
import pickle
from time import strftime, localtime
import numpy as np
import yaml
from ase import Atoms
from ase import units
from ase.optimize.optimize import Dynamics
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.neighborlist import NeighborList, natural_cutoffs
from gg.reference import get_ref_coeff
from gg.utils import NoReasonableStructureFound
from gg.utils_graph import atoms_to_graph, is_unique_graph


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
    ):
        """_summary_

        Args:
            atoms (ase.Atoms): The Atoms object to operate on.
            logfile (str, optional): Pickle file used to store the trajectory of atomic movement.
            Defaults to "gcbh.log".
            trajectory (str, optional):  If *logfile* is a string, a file will be opened.
            Defaults to "gcbh.traj".
            config_file (str[.yaml file], optional): Input file to gcbh. Defaults to None.
            restart (bool, optional):
        """
        # Intitalize by setting up the parent Dynamics Class
        super().__init__(atoms, logfile, trajectory)
        self.logfile.write("Begin GCBH Graph \n")
        self.logfile.flush()

        # Read Config File if it exists
        self.c = {
            "temp": 1500,
            "max_temp": None,
            "min_temp": None,
            "stop_steps": 400,
            "chemical_potential": None,
            "max_history": 25,
            "max_bond": 2.5,
            "max_bond_ratio": 0,
        }

        if config_file:
            self.set_config(config_file)

        # Setup Temperature
        if self.c["max_temp"] is None:
            self.c["max_temp"] = 1.0 / ((1.0 / self.c["temp"]) / 1.5)
        else:
            self.c["max_temp"] = max([self.c["max_temp"], self.c["temp"]])

        if self.c["min_temp"] is None:
            self.c["min_temp"] = 1.0 / ((1.0 / self.c["temp"]) * 1.5)
        else:
            self.c["min_temp"] = min([self.c["min_temp"], self.c["temp"]])

        # Some file names and folders are hardcoded
        self.current_atoms_name = "CONTCAR"
        self.status_file = "current_status.pkl"
        self.opt_folder = "opt_folder"
        self.lm_trajectory = Trajectory("local_minima.traj", "a", atoms)

        self.structure_modifiers = {}  # Setup empty class to add structure modifiers
        self.c["acc_hist"] = (
            []
        )  # used for adjusting the temperature of Metropolis algorithm
        # a series of 0 and 1, 0 stands for not accepted, 1 stands for accepted

        # Print the chemical potential for different elements
        if self.c["chemical_potential"]:
            self.mu = self.c["chemical_potential"]
            for k, v in self.mu.items():
                self.logtxt(f"Chemical potential of {k} is {v}")
        else:
            self.mu = None

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

    def todict(self):
        return self.__dict__

    def dump(self, filename: str):
        """dump dictionary of variables to store"""
        for key, value in self.structure_modifiers.items():
            value.atoms = None
            self.c["mod_weights"][key] = value.weight
        with open(filename, "wb") as file:
            pickle.dump(self.c, file, protocol=pickle.HIGHEST_PROTOCOL)

    def update_from_file(self, filename: str):
        """Update variable dictionary from file"""
        with open(filename, "rb") as file:
            c_old = pickle.load(file)
        self.c.update(c_old)
        if self.c["mod_weights"]:
            for key, value in self.c["mod_weights"].items():
                self.logtxt(f"{key} weight = {value:.2f}")

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

    def select_modifier(self):
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
            atoms (ase.Atoms): _description_
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

    def move(self, name: str):
        """Move atoms by a random modifier."""
        atoms = self.atoms
        self.logtxt(f"Modifier '{name}' formula {atoms.get_chemical_formula()}")
        atoms = self.structure_modifiers[name].get_modified_atoms(atoms)
        atoms.center()
        return atoms

    def initialize(self):
        """Initialize Atoms"""
        self.c["opt_on"] = 0
        self.atoms = self.optimize(self.atoms)
        self.dump(self.status_file)
        self.c["energy"] = self.atoms.get_potential_energy()
        ref = self.get_ref_potential(self.atoms)
        self.c["fe"] = self.c["energy"] - ref
        self.c["fe_min"] = self.c["fe"]
        self.c["opt_on"] = -1
        self.dump(self.status_file)
        self.c["nsteps"] += 1
        formula = self.atoms.get_chemical_formula()
        en = self.c["energy"]
        self.append_graph(self.atoms)
        self.logtxt(
            f'Atoms: {formula} E(initial): {en:.2f} F(initial) {self.c["fe"]:.2f}'
        )

    def run(self, steps: int = 4000, maximum_trial: int = 30):
        """Hop the basins for defined number of steps."""
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
                    if not isinstance(emsg, str):
                        emsg = "Unknown"
                    self.logtxt(
                        f"{modifier_name} did not find a good structure because of {emsg}"
                    )
                else:
                    if self.append_graph(newatoms):
                        self.c["opt_on"] = self.c["nsteps"]
                        self.logtxt(
                            f"One structure found with modifier {modifier_name}"
                        )
                        self.dump(self.status_file)
                        converged_atoms = self.optimize(newatoms)
                        _ = self.append_graph(converged_atoms)
                        en = converged_atoms.get_potential_energy()
                        self.logtxt(f"Optimization Done with E = {en:.2f}")
                        self.accepting_new_structures(converged_atoms, modifier_name)
                        self.c["opt_on"] = -1  # switch off the optimization status
                        self.dump(self.status_file)
                        self.c["nsteps"] += 1
                        break
                    else:
                        self.logtxt("Atoms object visited previously")
                        continue
            else:
                raise RuntimeError(
                    f"Program does not find a good structure after {trials+1} tests"
                )

    def accepting_new_structures(self, newatoms: Atoms, modifier_name: str):
        """This function takes care of all the accepting algorithms. I.E metropolis algorithms
        newatoms is the newly optimized structure
        """
        assert newatoms is not None
        en = newatoms.get_potential_energy()  # Energy_new
        fn = en - self.get_ref_potential(newatoms)  # Free_energy_new

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
            self.logtxt(f'Accepted, F(old)={self.c["fe"]:.2f} F(new)={fn:.2f}')
            self.atoms = newatoms
            self.c["energy"] = en
            self.c["fe"] = fn
            self.lm_trajectory.write(newatoms, accept=1)

        else:
            int_accept = 0
            self.logtxt(f'Rejected, F(old)={self.c["fe"]:.2f} F(new)={fn:.2f}')
            self.lm_trajectory.write(newatoms, accept=0)

        self.adjust_temp(int_accept)

        # update the best result for this basin-hopping
        if self.c["fe"] < self.c["fe_min"]:
            self.c["fe_min"] = self.c["fe"]
            self.c["no_impr_step"] = 0
        else:
            self.c["no_impr_step"] += 1

        self.dump(self.status_file)
        self.logtxt("-------------------------------------------------------")

    def optimize(self, atoms: Atoms, optimizer=BFGS, fmax: float = 0.05):
        """Optimize atoms"""
        if atoms.get_calculator() is None:
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
        atoms.write("POSCAR", format="vasp")
        dyn = optimizer(atoms, trajectory="opt.traj", logfile="opt_log")
        dyn.run(fmax=fmax)
        atoms.write("CONTCAR", format="vasp")
        self.logtxt(
            f'{get_current_time()}: Structure optimization completed for {self.c["nsteps"]}'
        )
        os.chdir(topdir)
        return atoms

    def append_graph(self, atoms):
        """Append the graph to list if its unique
        Args:
            atoms (_type_): _description_
        """
        nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
        nl.update(atoms)
        new_g = atoms_to_graph(
            atoms,
            nl,
            max_bond=self.c["max_bond"],
            max_bond_ratio=self.c["max_bond_ratio"],
        )
        if self.c["graphs"]:
            if is_unique_graph(new_g, self.c["graphs"]):
                self.c["graphs"].append(new_g)
                return True
            else:
                return False
        else:
            self.c["graphs"].append(new_g)
