Modifiers
=========

The modifiers form the building block of the code. They determine how the atoms are modified during each basin hopping step. 
The code provides basic modifiers as building blocks for more complex modifiers.

Add Monodentate
---------------

The modifier can add a monodentate adsorbate, or moiety at specific sites on the parent atoms object.

.. code-block:: python

  from gg.modifiers import Add
  from gg.sites import FlexibleSites

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  adsorbate_OH = read("OH.POSCAR") #adsorbate to be added
  add_OH = Add(FS, ads=adsorbate_OH, surf_coord=[1,2,3], ads_id=["O"], surf_sym=["Pt"],print_movie=True)

  atoms = read('POSCAR') #The atoms object that will adsorb
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
  from gg.sites import FlexibleSites

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  adsorbate_formate = read("OCHO.POSCAR") #adsorbate to be added (formate)
  add_formate = AddBi(FS, ads=adsorbate_OH, surf_coord=[1,2,3], ads_id=["O"], surf_sym=["Pt"], print_movie=True)

  atoms = read('POSCAR') #The atoms object that will adsorb
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
  
  adsorbate_OH = read("OH.POSCAR") #adsorbate to be removed
  remove_OH = Remove(FS, to_del=adsorbate_OH, print_movie=True)

  atoms = read('POSCAR_with_OH') #The atoms object that has OHs to be removed
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
  from gg.sites import FlexibleSites

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  adsorbate_OH = read("OH.POSCAR") #adsorbate to be removed
  adsorbate_NO = read("NO.POSCAR") #adsorbate to be replaced with
  remove_OH = Replace(FS, to_del=adsorbate_OH, with_replace= adsorbate_NO, print_movie=True)

  atoms = read('POSCAR_with_OH') #The atoms object that has OHs to be removed
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
  from gg.sites import FlexibleSites

  FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
  swap = Swap(FS, swap_sym=["Pt","Au"], print_movie=True)

  atoms = read('POSCAR_PtAu') #The atoms object with Pt and Au that can be swapped
  modified_atoms = swap.get_modified_atoms(atoms) #Atoms with the adsorbate

.. autoclass:: gg.modifiers.modifiers.Swap
  :members: get_modified_atoms
  :undoc-members:
  :show-inheritance:
