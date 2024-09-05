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
from gg.reference import get_ref_coeff
from gg.utils import NoReasonableStructureFound


__author__ = "Kaustubh Sawant, Geng Sun"


def get_current_time():
    """
    Returns:
        str: Current time
    """
    time_label = strftime("%d-%b-%Y %H:%M:%S", localtime())
    return time_label


class Gcbh(Dynamics):
    """ Basin hopping algorithm.

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
        self.config = {
            "temp": 1500,
            "max_temp": None,
            "min_temp": None,
            "stop_steps": 400,
            "chemical_potential": None,
            "max_history": 25,
        }

        if config_file:
            self.set_config(config_file)

        # Setup Temperature
        self.t = self.config["temp"]
        if self.config["max_temp"] is None:
            self.config["max_temp"] = 1.0 / ((1.0 / self.t) / 1.5)
        else:
            self.config["max_temp"] = max([self.config["max_temp"], self.t])

        if self.config["min_temp"] is None:
            self.config["min_temp"] = 1.0 / ((1.0 / self.t) * 1.5)
        else:
            self.config["min_temp"] = min([self.config["min_temp"], self.t])

        # Some file names and folders are hardcoded
        self.current_atoms_name = "CONTCAR"
        self.status_file = "current_status.pkl"
        self.opt_folder = "opt_folder"
        self.lm_trajectory = Trajectory("local_minima.traj", "a", atoms)

        self.structure_modifiers = {}  # Setup empty class to add structure modifiers
        self.accept_history = []  # used for adjusting the temperature of Metropolis algorithm
        # a series of 0 and 1, 0 stands for not accepted, 1 stands for accepted

        # Print the chemical potential for different elements
        if self.config["chemical_potential"]:
            self.mu = self.config["chemical_potential"]
            for k, v in self.mu.items():
                self.logtxt(f"Chemical potential of {k} is {v}")
        else:
            self.mu = None

        self.energy = None
        self.free_energy = None
        self.free_energy_min = None
        self.no_improvement_step = 0
        self.nsteps = 0

        # negative value indicates no ongoing structure optimization
        self.on_optimization = -1

        if restart:
            if os.path.exists(self.status_file):
                self.update_from_file(self.status_file)
        else:
            self.initialize()

    def set_config(self, config_file: str):
        """_
        Args:
            config_file (dict): Dictionary of key inputs
        """
        with open(config_file, "r", encoding="utf-8") as f:
            input_config = yaml.safe_load(f)
        self.config.update(input_config)

    def todict(self):
        d = self.__dict__.copy()
        keys_to_remove = [
            "logfile",
            "_lazy_cache",
            "optimizable",
        ]
        for key in keys_to_remove:
            d.pop(key, None)
        return d

    def dump(self, filename: str):
        """dump dictionary of variables to store"""
        if os.path.exists(filename):
            with open(filename, "rb") as file:
                new_state = pickle.load(file)
            for k, _ in new_state.items():
                new_state[k] = self.__dict__[k]
            with open(filename, "wb") as file:
                pickle.dump(new_state, file)
        else:
            with open(filename, "wb") as file:
                state = self.get_state()
                pickle.dump(state, file)

    def get_state(self):
        """get dictionary of variables to store"""
        state = self.__dict__.copy()
        # Remove the wierd file handle
        keys_to_remove = [
            "atoms",
            "logfile",
            "trajectory",
            "observers",
            "_lazy_cache",
            "optimizable",
            "lm_trajectory",
        ]
        for key in keys_to_remove:
            state.pop(key, None)
        return state

    def update_from_file(self, filename: str):
        """Update variable dictionary from file"""
        with open(filename, "rb") as file:
            new_state = pickle.load(file)
        self.__dict__.update(new_state)

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
        """ Start adjusting temperature beyond max_history
        Args:
            int_accept (int): 0 or 1
        """
        # adjust the temperatures
        self.accept_history.append(int_accept)
        if len(self.accept_history) > self.config["max_history"]:
            self.accept_history.pop(0)
            _balance = sum(self.accept_history) / float(self.config["max_history"])
            if _balance > 2.0 * (1 - _balance):
                self.t = self.t / 1.03
            elif _balance < 0.5 * (1 - _balance):
                self.t = self.t * 1.03

        if self.t < self.config["min_temp"]:
            self.t = self.config["min_temp"]
        elif self.t > self.config["max_temp"]:
            self.t = self.config["max_temp"]

    def move(self, name: str):
        """Move atoms by a random modifier."""
        atoms = self.atoms
        self.logtxt(f"Modifier '{name}' formula {atoms.get_chemical_formula()}")
        atoms = self.structure_modifiers[name].get_modified_atoms(atoms)
        atoms.center()
        return atoms

    def initialize(self):
        """ Initialize Atoms """
        self.on_optimization = 0
        self.atoms = self.optimize(self.atoms)
        self.dump(self.status_file)
        self.energy = self.atoms.get_potential_energy()
        ref = self.get_ref_potential(self.atoms)
        self.free_energy = self.energy - ref
        self.free_energy_min = self.free_energy
        self.on_optimization = -1
        self.dump(self.status_file)
        self.nsteps += 1
        formula = self.atoms.get_chemical_formula()
        self.logtxt(
            f"Atoms: {formula} E(initial): {self.energy:.2f} F(initial) {self.free_energy:.2f}"
        )

    def run(self, steps: int = 4000, maximum_trial: int = 30):
        """Hop the basins for defined number of steps."""
        self.logtxt("Intitial Weights:")
        for key, value in self.structure_modifiers.items():
            self.logtxt(f"{key} weight = {value.weight:.2f}")

        for step in range(steps):
            if self.no_improvement_step >= self.config["stop_steps"]:
                self.logtxt(
                    f"The best solution has not improved after {self.no_improvement_step} steps"
                )
                break
            self.logtxt(
                f"{step}-------------------------------------------------------"
            )

            self.logtxt(
                f"{get_current_time()}:  Starting Basin-Hopping Step {self.nsteps}"
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
                    self.on_optimization = self.nsteps
                    self.logtxt(f"One structure found with modifier {modifier_name}")
                    self.dump(self.status_file)
                    converged_atoms = self.optimize(newatoms)
                    en = converged_atoms.get_potential_energy()
                    self.logtxt(
                        f"{get_current_time()}: Optimization Done with E = {en:.2f}"
                    )
                    self.accepting_new_structures(converged_atoms, modifier_name)
                    self.on_optimization = -1  # switch off the optimization status
                    self.dump(self.status_file)
                    self.nsteps += 1
                    break
            else:
                raise RuntimeError(
                    f"Program does not find a good structure after {trials} tests"
                )

    def accepting_new_structures(self, newatoms: Atoms, modifier_name: str):
        """ This function takes care of all the accepting algorithms. I.E metropolis algorithms
        newatoms is the newly optimized structure
        """
        assert newatoms is not None
        en = newatoms.get_potential_energy()  # Energy_new
        fn = en - self.get_ref_potential(newatoms)  # Free_energy_new

        if fn < self.free_energy:
            accept = True
            modifier_weight_action = "increase"

        # Check Probability for acceptance
        elif np.random.uniform() < np.exp(-(fn - self.free_energy) / self.t / units.kB):
            accept = True
            modifier_weight_action = "decrease"
        else:
            accept = False
            modifier_weight_action = "decrease"

        self.update_modifier_weights(name=modifier_name, action=modifier_weight_action)

        if accept:
            int_accept = 1
            self.logtxt(f"Accepted, F(old)={self.free_energy:.2f} F(new)={fn:.2f}")
            self.atoms = newatoms
            self.energy = en
            self.free_energy = fn

        else:
            int_accept = 0
            self.logtxt(f"Accepted, F(old)={self.free_energy} F(new)={fn}")

        self.adjust_temp(int_accept)
        if accept:
            self.lm_trajectory.write(self.atoms, accept=1)
        else:
            self.lm_trajectory.write(self.atoms, accept=0)

        # update the best result for this basin-hopping
        if self.free_energy < self.free_energy_min:
            self.free_energy_min = self.free_energy
            self.no_improvement_step = 0
        else:
            self.no_improvement_step += 1

        self.dump(self.status_file)
        self.logtxt("-------------------------------------------------------")

    def optimize(self, atoms: Atoms, optimizer=BFGS, fmax: float = 0.05):
        """ Optimize atoms"""
        if atoms.get_calculator() is None:
            raise RuntimeError("The atoms object has no calculator")

        self.logtxt(
            f"{get_current_time()}: Begin structure optimization subroutine at step {self.nsteps}"
        )
        opt_dir = self.opt_folder
        topdir = os.getcwd()
        subdir = os.path.join(topdir, opt_dir, f"opt_{self.nsteps}")
        if not os.path.isdir(subdir):
            os.makedirs(subdir)
        os.chdir(subdir)
        atoms.write("POSCAR", format="vasp")
        dyn = optimizer(atoms, trajectory="opt.traj", logfile="opt_log")
        dyn.run(fmax=fmax)
        atoms.write("CONTCAR", format="vasp")
        self.logtxt(
            f"{get_current_time()}: Structure optimization completed for {self.nsteps}"
        )
        os.chdir(topdir)
        return atoms
