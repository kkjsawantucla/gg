from ase.optimize.optimize import Dynamics
from ase import units
from ase.io import read, write
import json
from time import strftime, localtime
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import PropertyNotImplementedError
from ase.io.trajectory import Trajectory
import subprocess
import os
import sys
import shutil
from distutils.version import LooseVersion

import numpy as np
assert LooseVersion(np.version.version) > LooseVersion("1.7.0")

class PropertyNotImplementedError(NotImplementedError):
    pass


class NoReasonableStructureFound(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class UnreasonableStructureFound(Exception):
    pass


class FragmentedStructure(Exception):
    pass


def get_current_time():
    time_label = strftime("%d-%b-%Y %H:%M:%S", localtime())
    return time_label


class GrandCanonicalBasinHopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
                 temperature=1500.0,
                 maximum_temp=None,
                 minimum_temp=None,
                 stop_steps=400,
                 logfile='grandcanonical.log',
                 trajectory='grandcanonical.traj',
                 local_minima_trajectory='local_minima.traj',
                 adjust_cm=False,
                 restart=False,
                 chemical_potential=None,
                 bash_script="optimize.sh",
                 files_to_copied=None,
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
        self.T = temperature
        if maximum_temp is None:
            self.max_T = 1.0/((1.0/self.T)/1.5)
        else:
            self.max_T = max([maximum_temp, self.T])
        if minimum_temp is None:
            self.min_T = 1.0/((1.0/self.T)*1.5)
        else:
            self.min_T = min([minimum_temp, self.T])
        self.stop_steps = stop_steps
        self.restart = restart
        self.bash_script = bash_script
        self.copied_files = files_to_copied

        # some file names and folders are hardcoded
        self.fn_current_atoms = "Current_atoms.traj"
        self.fn_status_file = "Current_Status.json"
        self.opt_folder = "opt_folder"

        self.structure_modifiers={}

        self.adjust_cm = adjust_cm

        self.lm_trajectory = local_minima_trajectory
        if isinstance(local_minima_trajectory, str):
            self.lm_trajectory = Trajectory(local_minima_trajectory,'a', atoms)

        # Dynamics.__init__ simply set
        Dynamics.__init__(self, atoms, logfile, trajectory)

        # print the program logo at the beginning of the output file
        self.logfile.write("%s\n" % programlogo)
        self.logfile.flush()

        # setup the chemical potential for different elements
        self.mu={}
        if chemical_potential is not None and os.path.isfile(chemical_potential):
            with open(chemical_potential, "r") as fp:
                for i, istr in enumerate(fp):
                    if istr.strip() == "":
                        continue
                    k, v=istr.split()
                    self.mu[k]=float(v)
        else:
            raise RuntimeError("chemical potential file %s is not found" % chemical_potential)
        for k, v in self.mu.items():
            self.dumplog("Chemical potential of %s is %.3f" % (k, v))

        # try to read previous result
        if self.restart:
            if (not os.path.isfile(self.fn_status_file)) or (not os.path.isfile(self.fn_current_atoms)):
                self.dumplog("%s or %s no found, start from scratch\n"
                                   % (self.fn_current_atoms,self.fn_status_file))
                self.restart=False
            elif os.path.getsize(self.fn_current_atoms) == 0:
                self.dumplog("{} is empty, set self.restart=False".format(self.fn_current_atoms))
                self.restart = False
            else:
                try:
                    atoms=read(self.fn_current_atoms)
                    atoms.get_potential_energy()
                except PropertyNotImplementedError:
                    self.dumplog("No energy found in {}, set self.restart=False".format(self.fn_current_atoms))
                    self.restart = False
                except RuntimeError as e:
                    self.dumplog("Error when read {}, set self.restart=False".format(e))
                    self.restart = False

        self.energy = None
        self.free_energy = None
        self.energy_min = None
        self.free_energy_min = None
        self.no_improvement_step = 0
        # negative value indicates no on-going structure optimization, otherwise it will be the  on-going optimization
        self.on_optimization = -1

        # this is used for adjusting the temperature of Metropolis algorithm
        self.accept_history = [] # a series of 0 and 1, 0 stands for not accpeted, 1 stands for accepted
        self.max_history = 25 # max length of self.accept_history is 25

        if not self.restart:
            self.initialize()
        else:
            self.reload_previous_results()

    def todict(self):
        d = {}
        return d

    def dumplog(self, msg="", level=1, highlight=None):
        if level < 1:
            level = 1
        real_message = " " * level + msg.strip() +"\n"
        if highlight is None:
            self.logfile.write(real_message)
        else:
            bars=highlight * (len(real_message)-1) + "\n"
            self.logfile.write(bars)
            self.logfile.write(real_message)
            self.logfile.write(bars)
        self.logfile.flush()

    def initialize(self):
        self.on_optimization = 0
        self.nsteps = 0
        self.optimize(self.atoms)
        self.save_current_status()
        self.energy = self.atoms.get_potential_energy()
        ref = self.get_ref_potential(self.atoms)
        self.free_energy = self.energy - ref
        self.energy_min = self.energy
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
        spc=SinglePointCalculator(t,energy=e,forces=f)
        t.set_calculator(spc)
        write(self.fn_current_atoms, t)

        accept_digits = ""
        for ii in self.accept_history:
            accept_digits += str(ii)
            accept_digits += ","
        accept_digits = accept_digits[:-1]

        # save the current status of the basin hopping
        info = {"nsteps": self.nsteps,
                "no_improvement_step": self.no_improvement_step,
                'Temperature': self.T,
                "free_energy_min": self.free_energy_min,
                "energy_min": self.energy_min,
                'history': accept_digits,
                'on_optimization': self.on_optimization}
        with open(self.fn_status_file, "w") as fp:
            json.dump(info, fp, sort_keys=True, indent=4, separators=(',', ': '))

    def reload_previous_results(self):
        with open(self.fn_status_file) as fp:
            info = json.load(fp)
            for k, v in info.items():
                if hasattr(v, 'keys'):
                    # if v is also a dictionary, which is used for recording the weights of operators; but they are not
                    # saved in the current version
                    self.dumplog("Read in {}".format(k))
                    for sub_k, sub_v in v.items():
                        self.dumplog("{0}={1}".format(sub_k, sub_v), level=4)
                else:
                    self.dumplog("Read previous result {0} ={1}".format(k, v))
            tl = get_current_time()
            self.dumplog("### %s: Previous Status Read in Successfullly ###\n" % tl)
            self.nsteps = info['nsteps']
            self.no_improvement_step = info['no_improvement_step']
            self.free_energy_min = info['free_energy_min']
            self.energy_min = info['energy_min']
            # Temperature and history is collected
            # since some previous version does not have this two terms, we have to query about the existence.
            if "Temperature" in info.keys():
                self.dumplog("Previous temperature is read\n")
                self.T = info['Temperature']
            if 'history' in info.keys():
                for ii in info['history'].split(","):
                    if ii.isdigit():
                        self.accept_history.append(int(ii))
            if 'on_optimization' in info.keys():
                self.on_optimization = info['on_optimization']

        previous_atoms = read(self.fn_current_atoms)
        self.update_self_atoms(previous_atoms)
        # get the self.energy and self.free_energy
        self.energy = self.atoms.get_potential_energy()
        ref = self.get_ref_potential(self.atoms)
        self.free_energy = self.energy - ref
        self.dumplog("self.atoms read successfully")

        # try to relocate previous optimization result
        if self.on_optimization > -1:
            opt_folder = os.path.join(os.getcwd(), self.opt_folder, "opt_%05d" % self.on_optimization)
            assert os.path.isdir(opt_folder)
            self.nsteps = self.on_optimization
            a = previous_atoms.copy()
            self.save_current_status()
            self.optimize(inatoms=a)
            self.accepting_new_structures(newatoms=a)
            self.on_optimization = -1
            self.save_current_status()
            self.nsteps += 1
        else:
            self.dumplog("Start new optimization from current atoms")

    def add_modifier(self, instance, name="mutation"):
        if name in self.structure_modifiers.keys():
            raise RuntimeError("structure modifier %s exists already!\n" % name)
        self.structure_modifiers[name] = [instance]

    def select_modifier(self):
        operator_names = self.structure_modifiers.keys()
        if not isinstance(operator_names, list):
            operator_names = list(operator_names)
        operator_weights = np.asarray([self.structure_modifiers[key].weight for key in operator_names])
        operator_weights = operator_weights/operator_weights.sum()
        return np.random.choice(operator_names, p=operator_weights)

    def update_modifier_weights(self, name='mutation', action='increase'):
        if name not in self.structure_modifiers.keys():
            raise RuntimeError("operator name %s not recognized" % name)
        if action not in ["increase", "decrease", "reset"]:
            raise RuntimeError("action must be 'increase','decrease' or 'rest'")
        elif action == "reset":
            self.structure_modifiers[name].weight = self.structure_modifiers[name].og_weight     
            self.dumplog("All the modifier weights are reset as 1.0\n")
        elif action == 'increase':
            self.structure_modifiers[name][-1].weight = min([self.structure_modifiers[name].og_weight*2, self.structure_modifiers[name].weight*1.05])
        else:
            self.structure_modifiers[name][-1] = max([self.structure_modifiers[name].og_weight/2.0, self.structure_modifiers[name].weight/1.05])

    def move(self, name='mutation'):
        """Move atoms by a random step."""
        atoms = self.atoms.copy()
        self.dumplog("%s : Starting operator '%s' (formula %s) \n" %
                           (get_current_time(), name, atoms.get_chemical_formula()))
        atoms = self.structure_modifiers[name].get_modified_atoms(atoms)
        atoms.center()
        self.dumplog("%s : End operator (formula %s) \n" % (get_current_time(),atoms.get_chemical_formula()))
        return atoms

    def log_status(self):
        time_label = get_current_time()
        natoms = self.atoms.get_number_of_atoms()
        formula = self.atoms.get_chemical_formula()
        self.dumplog("%20s%6s (natoms=%3d, %8s) Steps:%8d E=%15.8f F=%15.8f \n" %
                               (time_label, "GCBH", natoms, formula, self.nsteps - 1,
                                self.energy, self.free_energy))
        for key in self.structure_modifiers.keys():
            self.dumplog("modifier %s (weight %3.2f)    " % (key, self.structure_modifiers[key].weight))
        self.dumplog("Current Temperature is %.2f" % self.T)

    def run(self, maximum_steps=4000, maximum_trial=30):
        """Hop the basins for defined number of steps."""
        while self.nsteps < maximum_steps:
            if self.no_improvement_step >= self.stop_steps:
                self.dumplog("The best solution has not "
                             "improved after {} steps\n".format(self.no_improvement_step), highlight="#")
            self.dumplog("-------------------------------------------------------")
            time_label = get_current_time()
            self.dumplog("%s:  Starting Basin-Hopping Step %05d\n" % (time_label, self.nsteps))

            for number_of_trials in range(maximum_trial):
                modifier_name = self.select_modifier()
                try:
                    new_atoms = self.move(modifier_name=modifier_name)
                except NoReasonableStructureFound as emsg:  # emsg stands for error message
                    if not isinstance(emsg, str):
                        emsg = "Unknown"
                    self.dumplog("%s did not find a good structure because of %s" % (modifier_name, emsg))
                else:
                    self.on_optimization = self.nsteps
                    self.dumplog("One structure found, begin to optimize this structure\n")

                    self.save_current_status()  
                    self.optimize(inatoms=new_atoms)
                   
                    self.accepting_new_structures(newatoms=new_atoms, move_action=modifier_name)
                    self.on_optimization = -1
                    
                    self.save_current_status()
                    self.nsteps += 1
                    break
            else:
                raise RuntimeError("Program does not find a good structure after {} tests".format(maximum_trial))

    def accepting_new_structures(self, newatoms=None, move_action=None):
        """This function takes care of all the accepting algorithm. I.E metropolis algorithms
        newatoms is the newly optimized structure
        move_action is action (modifier name) to  produce the initial structure for newatoms;
        If move_action is specified, its weights will be adjusted according to the acception or rejection; otherwise,
        the weights are not altered"""

        assert newatoms is not None

        En = newatoms.get_potential_energy()  # Energy_new
        Fn = En-self.get_ref_potential(newatoms) # Free_energy_new

        accept = False
        modifier_weight_action = 'decrease'
        if Fn < self.free_energy:
            accept = True
            modifier_weight_action = 'increase'
        elif np.random.uniform() < np.exp(-(Fn-self.free_energy)/self.T/units.kB):
            accept = True

        if move_action is not None:
            self.update_modifier_weights(name=move_action, action=modifier_weight_action)

        if accept:
            _int_accept=1
            self.dumplog("Accepted, F(old)=%.3f F(new)=%.3f\n" % (self.free_energy, Fn))
            self.update_self_atoms(newatoms)
            self.energy = En
            self.free_energy = Fn
            # if move_action is not None:
            #     self.update_modifier_weights(name=move_action, action='increase')
        else:
            _int_accept=0
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
        if len(self.accept_history) > self.max_history:
            self.accept_history.pop(0)
            _balance = sum(self.accept_history)/float(self.max_history)
            if _balance > 2.0* (1-_balance):
                self.T = self.T/1.03
            elif _balance < 0.5* (1-_balance):
                self.T = self.T*1.03

        if self.T < self.min_T:
            self.T = self.min_T
        elif self.T > self.max_T:
            self.T = self.max_T

        # update the best result for this basin-hopping
        if self.free_energy < self.free_energy_min:
            self.free_energy_min = self.free_energy
            self.no_improvement_step = 0
        else:
            self.no_improvement_step += 1

        # self.energy is not used for updating no_improvement_step
        if self.energy < self.energy_min:
            self.energy_min = self.energy

        # self.log_status()
        self.save_current_status()
        self.log_status()
        self.dumplog("-------------------------------------------------------")

    def optimize(self, inatoms=None, restart=False):
        self.dumplog("{}: begin structure optimization subroutine at step {}".format(get_current_time(), self.nsteps))
        atoms = inatoms.copy()
        opt_dir = self.opt_folder
        steps = self.nsteps
        script = self.bash_script
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
            opt_job = subprocess.Popen(['bash', script], cwd=subdir)
            opt_job.wait()
            if opt_job.returncode < 0:
                sys.stderr.write("optimization does not terminate properly at {}".format(subdir))
                sys.exit(1)
        except:
            raise RuntimeError("some error encountered at folder {} during optimizations".format(subdir))
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

    def get_ref_potential(self,atoms):
        """
        calculate the chemical potential of atoms
        :param atoms:
        :return:
        """
        ref=0.0
        for i,si in enumerate(atoms.get_chemical_symbols()):
            if si not in self.mu.keys():
                raise RuntimeError("I did not find the chemical potential for element %s" % si)
            else:
                ref += self.mu.get(si)
        return ref

    def update_self_atoms(self, a):
        """
        This function will keep the original reference of self.atoms, but refresh it with new structures.
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
            e=a.get_potential_energy()
            f=a.get_forces()
        except PropertyNotImplementedError:
            self.dumplog("Warnning : self.atoms no energy !!!!")
        else:
            spc=SinglePointCalculator(self.atoms,forces=f, energy=e)
            self.atoms.set_calculator(spc)
 
