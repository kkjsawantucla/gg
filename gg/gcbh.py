"""Importing Modules"""

import subprocess
import os
import sys
import shutil
import json
from time import strftime, localtime
import numpy as np
import yaml
from ase import units
from ase.optimize.optimize import Dynamics
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import PropertyNotImplementedError
from gg.reference import get_ref_coeff

__author__ = "Kaustubh Sawant, Geng Sun"

def get_current_time():
    """
    Returns:
        str: time
    """
    time_label = strftime("%d-%b-%Y %H:%M:%S", localtime())
    return time_label

class GrandCanonicalBasinHopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(
        self,
        atoms,
        logfile="grandcanonical.log",
        trajectory="grandcanonical.traj",
        config_file=None,
    ):
        """Parameters:

        atoms: Atoms object
            The Atoms object to operate on.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.
        """
        #Intitalize by setting up the parent Dynamics Class
        super().__init__(atoms, logfile, trajectory)
        self.logfile.write("Begin GCBH Graph \n")
        self.logfile.flush()

        #Read Config File if it exists
        self.config = {
            "temp": 1500,
            "max_temp": None,
            "min_temp": None,
            "stop_steps": 400,
            "restart": False,
            "chemical_potential": None,
            "bash_script": "optimize.sh",
            "files_to_copied": None,
            "max_history":25,
        }
        if config_file:
            self.set_config(config_file)

       
        #Setup Temperature
        self.t = self.config["temp"]
        if self.config["max_temp"] is None:
            self.config["max_temp"] = 1.0 / ((1.0 / self.t) / 1.5)
        else:
            self.config["max_temp"] = max([self.config["max_temp"], self.t])

        if self.config["min_temp"] is None:
            self.config["min_temp"] = 1.0 / ((1.0 / self.t) * 1.5)
        else:
            self.config["min_temp"] = min([self.config["min_temp"], self.t])

        #Some file names and folders are hardcoded
        self.current_atoms = "Current_atoms.traj"
        self.status_file = "Current_Status.json"
        self.opt_folder = "opt_folder"
        self.structure_modifiers = {} #Setup empty class to add structure modifiers
        self.accept_history = [] # this is used for adjusting the temperature of Metropolis algorithm
        # a series of 0 and 1, 0 stands for not accpeted, 1 stands for accepted

        # Print the chemical potential for different elements
        if self.config["chemical_potential"]:
            self.mu = self.config["chemical_potential"]
            for k, v in self.mu.items():
                self.dumplog(f"Chemical potential of {k} is {v}")
        else:
            self.mu = None

        self.energy = None
        self.free_energy = None
        self.free_energy_min = None
        self.no_improvement_step = 0
        
        # negative value indicates no on-going structure optimization,
        # otherwise it will be the on-going optimization
        self.on_optimization = -1
        self.initialize()

    def set_config(self, config_file):
        """_
        Args:
            config_file (dict): Dictionary of key inputs
        """
        with open(config_file, "r", encoding="utf-8") as f:
            input_config = yaml.safe_load(f)
        self.config.update(input_config)

    def todict(self):
        d = {}
        return d

    def dumplog(
        self,
        msg="",
    ):
        """Dump into logfile
        Args:
            msg (str, optional): The message to dump. Defaults to "".
        """
        real_message = f"{msg.strip()} + \n"
        self.logfile.write(real_message)
        self.logfile.flush()

    def initialize(self):
        self.on_optimization = 0
        self.nsteps = 0
        self.optimize(self.atoms)
        self.save_current_status()
        self.energy = self.atoms.get_potential_energy()
        ref = self.get_ref_potential(self.atoms)
        self.free_energy = self.energy - ref
        self.free_energy_min = self.free_energy
        self.no_improvement_step = 0
        self.on_optimization = -1
        self.save_current_status()
        self.nsteps += 1

    def save_current_status(self):
        # save current atoms
        t = self.atoms.copy()
        t.info = self.atoms.info.copy()
        e = self.atoms.get_potential_energy()
        f = self.atoms.get_forces()
        spc = SinglePointCalculator(t, energy=e, forces=f)
        t.set_calculator(spc)
        write(self.current_atoms, t)
        accept_digits = ""
        for ii in self.accept_history:
            accept_digits += str(ii)
            accept_digits += ","
        accept_digits = accept_digits[:-1]

        # save the current status of the basin hopping
        info = {
            "nsteps": self.nsteps,
            "no_improvement_step": self.no_improvement_step,
            "Temperature": self.t,
            "free_energy_min": self.free_energy_min,
            "history": accept_digits,
            "on_optimization": self.on_optimization,
        }
        with open(self.status_file, "w", encoding="utf-8") as fp:
            json.dump(info, fp, sort_keys=True, indent=4, separators=(",", ": "))

    def add_modifier(self, instance, name):
        """
        Args:
            instance (Modifier): Instance of a ParentModifier or a child
            name (str): Name for the instance/modifier
        """
        if name in self.structure_modifiers:
            raise RuntimeError(f"Structure modifier {name} exists already!\n")
        self.structure_modifiers[name] = [instance]
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

    def update_modifier_weights(self, name, action):
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
            self.dumplog(
                f"Modifier {name} weight is reset to original weight : {og_weight}"
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
        return

    def move(self, name):
        """Move atoms by a random step."""
        atoms = self.atoms.copy()
        self.dumplog(
            f"{get_current_time()} : Modifier '{name}' formula {atoms.get_chemical_formula()}"
        )
        atoms = self.structure_modifiers[name].get_modified_atoms(atoms)
        atoms.center()
        return atoms

    def log_status(self):
        """Print Current Status Log"""
        time_label = get_current_time()
        natoms = self.atoms.get_number_of_atoms()
        formula = self.atoms.get_chemical_formula()
        self.dumplog(
            "%20s%6s (natoms=%3d, %8s) Steps:%8d E=%15.8f F=%15.8f"
            % (
                time_label,
                "GCBH",
                natoms,
                formula,
                self.nsteps - 1,
                self.energy,
                self.free_energy,
            )
        )
        for key, value in self.structure_modifiers:
            self.dumplog(f"modifier:{key} with weight:{value.weight}")
            self.dumplog(f"Current temperature is {self.t}")

    def run(self, steps=4000, maximum_trial=30):
        """Hop the basins for defined number of steps."""
        for step in range(steps):
            if self.no_improvement_step >= self.config["stop_steps"]:
                self.dumplog(
                    f"The best solution has not improved after {self.no_improvement_step} steps"
                )
                break
            self.dumplog(f"{step}-------------------------------------------------------")

            self.dumplog(f"{get_current_time()}:  Starting Basin-Hopping Step {self.nsteps}")

            for number_of_trials in range(maximum_trial):
                modifier_name = self.select_modifier()
                new_atoms = self.move(modifier_name)
                self.on_optimization = self.nsteps
                self.dumplog("One structure found, begin to optimize this structure\n")

                self.save_current_status()
                self.optimize(inatoms=new_atoms)

                self.accepting_new_structures(
                    newatoms=new_atoms, move_action=modifier_name
                )
                self.on_optimization = -1

                self.save_current_status()
                self.nsteps += 1

    def accepting_new_structures(self, newatoms=None, move_action=None):
        """This function takes care of all the accepting algorithm. I.E metropolis algorithms
        newatoms is the newly optimized structure
        move_action is action (modifier name) to  produce the initial structure for newatoms;
        If move_action is specified, its weights will be adjusted according to the acception or rejection; otherwise,
        the weights are not altered"""

        assert newatoms is not None

        En = newatoms.get_potential_energy()  # Energy_new
        Fn = En - self.get_ref_potential(newatoms)  # Free_energy_new

        accept = False
        modifier_weight_action = "decrease"
        if Fn < self.free_energy:
            accept = True
            modifier_weight_action = "increase"
        elif np.random.uniform() < np.exp(-(Fn - self.free_energy) / self.t / units.kB):
            accept = True

        if move_action is not None:
            self.update_modifier_weights(
                name=move_action, action=modifier_weight_action
            )

        if accept:
            _int_accept = 1
            self.dumplog("Accepted, F(old)=%.3f F(new)=%.3f\n" % (self.free_energy, Fn))
            self.update_self_atoms(newatoms)
            self.energy = En
            self.free_energy = Fn
            # if move_action is not None:
            #     self.update_modifier_weights(name=move_action, action='increase')
        else:
            _int_accept = 0
            self.dumplog("Rejected, F(old)=%.3f F(new)=%.3f\n" % (self.free_energy, Fn))
            # if move_action is not None:
            #     self.update_modifier_weights(name=move_action, action='decrease')

        # if accept and self.lm_trajectory is not None:
        #     self.lm_trajectory.write(self.atoms)
        if accept:
            self.lm_trajectory.write(self.atoms, accept=1)
        else:
            self.lm_trajectory.write(self.atoms, accept=0)

        # adjust the temperatures
        self.accept_history.append(_int_accept)
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

        # update the best result for this basin-hopping
        if self.free_energy < self.free_energy_min:
            self.free_energy_min = self.free_energy
            self.no_improvement_step = 0
        else:
            self.no_improvement_step += 1

        # self.log_status()
        self.save_current_status()
        self.log_status()
        self.dumplog("-------------------------------------------------------")

    def optimize(self, inatoms=None, restart=False):
        self.dumplog(
            "{}: begin structure optimization subroutine at step {}".format(
                get_current_time(), self.nsteps
            )
        )
        atoms = inatoms.copy()
        opt_dir = self.opt_folder
        steps = self.nsteps
        script = self.config["bash_script"]
        copied_files = self.copied_files[:]
        topdir = os.getcwd()
        subdir = os.path.join(topdir, opt_dir, "opt_%05d" % steps)
        if restart:
            assert os.path.isdir(subdir)
        else:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

            if script not in copied_files:
                copied_files.append(script)
            for fn in copied_files:
                assert os.path.isfile(fn)
                shutil.copy(os.path.join(topdir, fn), os.path.join(subdir, fn))
            write(os.path.join(subdir, "input.traj"), atoms)
        try:
            os.chdir(subdir)
            opt_job = subprocess.Popen(["bash", script], cwd=subdir)
            opt_job.wait()
            if opt_job.returncode < 0:
                sys.stderr.write(
                    "optimization does not terminate properly at {}".format(subdir)
                )
                sys.exit(1)
        except:
            raise RuntimeError(
                "some error encountered at folder {} during optimizations".format(
                    subdir
                )
            )
        else:
            fn = os.path.join(subdir, "optimized.traj")
            assert os.path.isfile(fn)
            optimized_atoms = read(fn)
        finally:
            os.chdir(topdir)

        e = optimized_atoms.get_potential_energy()
        f = optimized_atoms.get_forces()

        cell = optimized_atoms.get_cell()
        pbc = optimized_atoms.get_pbc()
        inatoms.set_constraint()
        del inatoms[range(inatoms.get_number_of_atoms())]
        inatoms.extend(optimized_atoms)
        inatoms.set_pbc(pbc)
        inatoms.set_cell(cell)
        inatoms.set_constraint(optimized_atoms.constraints)
        spc = SinglePointCalculator(inatoms, energy=e, forces=f)
        inatoms.set_calculator(spc)
        self.dumplog("{}: Optimization Done\n".format(get_current_time()))

    def get_ref_potential(self, atoms):
        """
        Args:
            atoms (ase.Atoms): _description_
        Returns:
            float: total ref value to substract
        """
        if self.mu:
            formula = atoms.get_chemical_formula()
            ref_sum = 0
            to_print = f"{formula} = "
            ref_coeff = get_ref_coeff(self.mu, formula)
            for key,value in self.mu.items():
                ref_sum += ref_coeff[key]*value
                to_print += f"{ref_coeff[key]} {value} +"
            self.dumplog(to_print)
            return ref_sum
        else:
            return 0

    def update_self_atoms(self, a):
        """
        This function will keep the original reference of self.atoms,
        but refresh it with new structures.

        You have to keep the reference of self.atoms, otherwise, self.call_observers will not work.
        :param a: ase.atoms.Atoms object.
        :return: None
        """
        self.atoms.set_constraint()
        del self.atoms[range(self.atoms.get_number_of_atoms())]
        cell = a.get_cell()
        pbc = a.get_pbc()
        self.atoms.extend(a.copy())
        self.atoms.set_pbc(pbc)
        self.atoms.set_cell(cell)
        self.atoms.set_constraint(a.constraints)
        try:
            e = a.get_potential_energy()
            f = a.get_forces()
        except PropertyNotImplementedError:
            self.dumplog("Warnning : self.atoms no energy !!!!")
        else:
            spc = SinglePointCalculator(self.atoms, forces=f, energy=e)
            self.atoms.set_calculator(spc)
