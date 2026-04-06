Modifiers
=========

.. contents::
    :local:

The modifiers form the building block of the code. They determine how the atoms are modified during each basin hopping step. 
The code provides basic modifiers as building blocks for more complex modifiers.

Many modifiers support ``unique_method`` and ``unique_depth`` for controlling how duplicate structures are filtered.

Add Monodentate
---------------

The modifier can add a monodentate adsorbate, or moiety at specific sites on the parent atoms object.

.. code-block:: python

  from gg.modifiers import Add
  from gg.predefined_sites import FlexibleSites
  from ase.build import fcc111

  atoms = fcc111("Pt", size=(3, 3 , 4), vacuum=10.0)

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  add_OH = Add(FS, ads="OH", surf_coord=[1,2,3], ads_id=["O"], surf_sym=["Pt"], print_movie=True, unique=True)

  modified_atoms = add_OH.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.add.Add
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Add Bidentate
---------------

The modifier can add a bidentate adsorbate, or moiety at specific sites on the parent atoms object.

.. code-block:: python

  from gg.modifiers import AddBi
  from gg.predefined_sites import FlexibleSites
  from ase.build import fcc111

  atoms = fcc111("Pt", size=(3, 3 , 4), vacuum=10.0)

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  add_formate = AddBi(FS, ads="HCOO", surf_coord=[1,2,3], ads_id=["O"], surf_sym=["Pt"], print_movie=True, unique=True)

  modified_atoms = add_formate.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.add.AddBi
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Remove Adsorbate
----------------

Remove an adsorbate from the surface

.. code-block:: python

  from gg.modifiers import Remove
  from gg.predefined_sites import FlexibleSites

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface  
  remove_OH = Remove(FS, to_del="OH", print_movie=True)

  atoms = read("POSCAR_with_OH") #The atoms object that has OHs to be removed
  modified_atoms = remove_OH.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.modifiers.Remove
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Replace Atoms/Molecules
-----------------------

Remove an adsorbate from the surface

.. code-block:: python

  from gg.modifiers import Replace
  from gg.predefined_sites import FlexibleSites
  from ase.io import read

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  remove_OH = Replace(FS, to_del="OH", with_replace="NO", print_movie=True)

  atoms = read("POSCAR_with_OH") #The atoms object that has OHs to be removed
  modified_atoms = remove_OH.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.modifiers.Replace
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Swap Atoms
-----------------------

Swap two atoms on the surface

.. code-block:: python

  from gg.modifiers import Swap
  from gg.predefined_sites import FlexibleSites
  from ase.io import read

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  swap = Swap(FS, swap_sym=["Pt","Au"], print_movie=True)

  atoms = read("POSCAR_PtAu") #The atoms object with Pt and Au that can be swapped
  modified_atoms = swap.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.modifiers.Swap
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Cluster Rotate
-----------------------

Rotate a cluster of atoms

Make sure to tag the cluster atoms using atom.tag = -1

.. code-block:: python

  from gg.modifiers import ClusterRotate
  from gg.predefined_sites import FlexibleSites
  from ase.io import read

  atoms = read("POSCAR_Pt_TiO2")
  for a in atoms:
    if a.symbol=='Pt':
      a.tag=-1
  FS = FlexibleSites(tag=-1)
  rotate = ClusterRotate(FS,contact_error=0.2,rotate_vector=(0,0,1))
  modified_atoms = rotate.get_modified_atoms(atoms)

.. autoclass:: gg.modifiers.cluster.ClusterRotate
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Cluster Translate
-----------------------

Translate a cluster of atoms

Make sure to tag the cluster atoms using atom.tag = -1

.. code-block:: python

  from gg.modifiers import ClusterTranslate
  from gg.predefined_sites import FlexibleSites
  from ase.io import read

  atoms = read("POSCAR_Pt_TiO2")
  for a in atoms:
    if a.symbol=='Pt':
      a.tag=-1
  FS = FlexibleSites(tag=-1)
  trans = ClusterTranslate(FS,contact_error=0.2)
  modified_atoms = trans.get_modified_atoms(atoms)

.. autoclass:: gg.modifiers.cluster.ClusterTranslate
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:


Modifier Adder
-----------------------

Combine multiple modifiers into a single sequential modifier. The following example shows how two add modifiers can be combined to create a new dissociatively H2O adsorption modifier.

.. code-block:: python

  from gg.modifiers import Add, ModifierAdder
  from gg.predefined_sites import FlexibleSites
  from ase.build import fcc111

  atoms = fcc111("Pt", size=(3, 3 , 4), vacuum=10.0)

  FS = FlexibleSites(max_bond_ratio=1.2,com=0.5)
  add_H = Add(FS, ads="H", surf_coord=[1], ads_id=["H"], surf_sym=["Pt","O"],print_movie=True)
  add_OH = Add(FS, ads="OH", surf_coord=[1], ads_id=["O"], surf_sym=["Pt"],print_movie=True)
  add_H2O = ModifierAdder([add_OH, add_H], print_movie=True, unique=True)

  modified_atoms = add_H2O.get_modified_atoms(atoms)

.. autoclass:: gg.modifiers.modifiers.ModifierAdder
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:

Creating New Modifiers with ParentModifier
------------------------------------------

You can implement custom modifiers by inheriting from ``ParentModifier`` and defining
``get_modified_atoms``.

Use this pattern:

.. code-block:: python

  from ase import Atoms
  from gg.modifiers.modifiers import ParentModifier
  from gg.utils import NoReasonableStructureFound

  class RaiseTopLayer(ParentModifier):
      """Example custom modifier that shifts selected atoms in +z."""

      def __init__(self, z_shift: float = 0.1, weight: float = 1.0):
          super().__init__(weight)
          self.z_shift = z_shift

      def get_modified_atoms(self, atoms: Atoms) -> Atoms:
          # Always set self.atoms first so string/Atoms inputs are normalized
          self.atoms = atoms

          # Modify a copied Atoms object via self.atoms
          top_z = max(a.position[2] for a in self.atoms)
          for atom in self.atoms:
              if atom.position[2] > top_z - 0.5:
                  atom.position[2] += self.z_shift

          # Raise NoReasonableStructureFound when the move is invalid
          # raise NoReasonableStructureFound("No valid perturbation found")

          return self.atoms

Key points when implementing a new modifier:

- Call ``super().__init__(weight)`` in your custom ``__init__``.
- Set ``self.atoms = atoms`` at the top of ``get_modified_atoms``.
  ``ParentModifier`` handles both file paths and ``ase.Atoms`` inputs.
- Return an ``ase.Atoms`` object (or a list of them only for movie-like workflows).
- Raise ``NoReasonableStructureFound`` when your attempted move cannot produce a valid structure.

For modifiers that can return many candidates (``print_movie=True`` workflows), implement uniqueness filtering the same way as built-in modifiers:

.. code-block:: python

  from gg.utils_graph import get_unique_atoms

  if self.unique:
      return get_unique_atoms(
          movie,
          max_bond=self.max_bond,
          max_bond_ratio=self.max_bond_ratio,
          unique_method=self.unique_method,
          depth=self.unique_depth,
      )
  else:
      return movie

This pattern ensures your custom modifier honors ``unique_method`` and ``unique_depth`` consistently.

.. autoclass:: gg.modifiers.modifiers.ParentModifier
  :members: get_modified_atoms, atoms
  :undoc-members:
  :show-inheritance:

