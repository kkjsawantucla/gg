Examples
========

1. ZnO/Cu Cluster Using MACE
--------------------------------

The only file you need to initialize is `input.yaml <https://github.com/kkjsawantucla/gg/blob/main/examples/Zn_Cu_cluster/input.yaml>`_, which contains information about chemical potential and temperature.

    .. code-block:: yaml

        temp: 550
        stop_steps: 2000
        stop_opt: 500 #Max steps for geometry optimization
        
        chemical_potential:
          Cu: -4.10
          ZnO: -9.10
          H: -3.43
          O: -7.51 
        
        check_graphs: False #Check graphs can remember unique structures to avoid redundancy
        vib_correction: False 
        max_bond_ratio: 1.2
        initialize: True


First, we need to set up the atoms. We are using a model developed by `Kempen et. al. <https://www.nature.com/articles/s41524-024-01507-z>`_ and the interatomic potential developed by the `MACE team <https://github.com/ACEsuit/mace/tree/main?tab=readme-ov-file#pretrained-foundation-models>`_
    .. code-block:: python

        from ase.io import read
        atoms = read('POSCAR')
        
        from mace.calculators import mace_mp
        calc = mace_mp(model="medium-mpa-0",default_dtype="float64",device='cuda')
        atoms.calc = calc


Define how to choose the surface atoms for modifications
    .. code-block:: python

        #Define surface site class
        from gg.predefined_sites import FlexibleSites

        for a in atoms:
            if a.symbol in ['Zn','O','H']:
                a.tag=-1
        FS = FlexibleSites(tag=-1) #Define adsorbate/cluster
        FS2 = FlexibleSites(constraints=True,com=0.75) #Get Surface atoms based on z co-ordinate

Define possible surface modifications allowed during basin hopping
    .. code-block:: python

        # Build Modifiers for the system
        from gg.modifiers import Add, Remove, Replace, ModifierAdder, ClusterRotate, ClusterTranslate, 

        #Add, Remove, Swap O
        addO = Add(FS2, "O", surf_coord=[2,3], ads_id = ["O"], surf_sym = ["Cu","Zn"], unique_method = "O")
        remO = Remove(FS, "O", max_bond_ratio = 1.2, unique_method = "O")
        swapO = ModifierAdder([addO,remO],unique_method = "O")

        #Add, Remove, Swap H
        addH = Add(FS2, "H", surf_coord=[1,2,3], ads_id = ["H"], surf_sym = ["Cu","O"], unique_method = "H")
        remH = Remove(FS, "H", max_bond_ratio = 1.2, unique_method = "H")
        swapH = ModifierAdder([addH,remH],unique_method = "H")

        #Rotate, Translate Clusters
        clstr_rot = ClusterRotate(FS,contact_error=0.25,rotate_vector=(0,0,1),nmovie=2)
        clstr_trans = ClusterTranslate(FS,contact_error=0.2,nmovie=2)

        #Add ZnOH
        adsZnOH = read('POSCAR_ZnOH')
        addZnOH = Add(FS2, adsZnOH, surf_coord=[2,3], ads_id = ["Zn"], surf_sym = ["Cu"], unique_method = "H")

        #Add, Remove OH
        addOH = Add(FS2, "OH", surf_coord=[2,3], ads_id = ["O"], surf_sym = ["Cu","Zn"], unique_method = "H")
        remOH = Remove(FS2, "OH", max_bond_ratio = 1.2, unique_method = "H")

        #Remove H2O
        remH2O = Remove(FS, "H2O", max_bond_ratio = 1.2, unique_method = "H")

        #Move, Remove Zn
        repl_Cu = Replace(FS2,to_del="Cu",with_replace="Zn", unique_method = "Zn")
        remZn = Remove(FS, "Zn", max_bond_ratio = 1.2, unique_method = "Zn")
        repl_Zn = Replace(FS2,to_del="Zn",with_replace="Cu", unique_method = "Zn")

Initialize the GCBH
    .. code-block:: python

        from gg.gcbh import Gcbh
        G = Gcbh(atoms,config_file='input.yaml')
        G.add_modifier(addO,'Add O')
        G.add_modifier(remO,'Remove O')
        G.add_modifier(swapO,'Swap O')
        G.add_modifier(addH,'Add H')
        G.add_modifier(remH,'Remove H')
        G.add_modifier(swapH,'Swap H')
        G.add_modifier(clstr_trans,'Cluster Translate')
        G.add_modifier(clstr_rot,'Cluster Rotate')
        G.add_modifier(repl_Cu,'Replace Cu with Zn')
        G.add_modifier(repl_Zn,'Replace Zn with Cu')
        G.add_modifier(remZn,"Remove Zn")
        G.add_modifier(remOH,"Remove OH")
        G.add_modifier(addZnOH,"Add ZnOH")
        G.add_modifier(addOH,"Add OH")
        G.add_modifier(remH2O,"Rem H2O")
        G.add_delete_gas(gas_species=["H2"])

Sometimes, the simulation can generate gas-phase species, which can skew results
    .. code-block:: python

        G.add_delete_gas(gas_species=["H2"])

Finally, run the code
    .. code-block:: python

        G.run(steps=1000)

This should generate the following files and folders:

- **local_minima.traj** : Trajectory file of accepted structures.
- **gcbh.log** : Log of the run.
- **gcbh.traj** : Trajectory file of all structures.
- **current_status.pkl** : current status of the run, useful in restarting.
- **opt_folder** : Folder containing individual geometry optimization steps.

 - opt_00
 - opt_01
 - ...


2. Adding H2O to ASA Surface
-----------------------------

The addition of dissociative water on complex aluminosilicate surfaces can be achieved with just a few lines of code.

Defining the surface
    .. code-block:: python

        from gg.sites import RuleSites, get_com_sites, get_surface_sites_by_coordination

        #Define maximum co-odrination each species can have
        max_coord = {"Al": 6, "Si": 4, "O": 4, "H": 1}

        ss = RuleSites(
            index_parsers=[
                lambda atoms: get_com_sites(atoms, fraction=0.50, direction="above"),
                lambda atoms: get_surface_sites_by_coordination(
                    atoms, max_coord, max_bond=2,
                ),
            ],
            combine_rules="intersection",
        )

Defining the addH2O modifier
    .. code-block:: python

        from gg.modifiers import Add, ModifierAdder

        addH = Add(ss,"H",1,ads_id="H",surf_sym=["O"],print_movie=True,unique=True,unique_method="H")
        addOH = Add(ss,"OH",1,ads_id="O",surf_sym=["Al", "Si"],print_movie=True,unique=True,unique_method="H")
        addH2O = ModifierAdder([addOH, addH],print_movie=True,unique=True,max_bond=2,unique_method="H")

Reading POSCAR and applying add modifier
    .. code-block:: python

        from ase.io import read

        atoms = read("POSCAR")
        add_atoms = addH2O.get_modified_atoms(atoms)
