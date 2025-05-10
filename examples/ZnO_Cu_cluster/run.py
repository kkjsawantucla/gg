from ase.io import read
from mace.calculators import mace_mp
from gg.gcbh import Gcbh
from gg.predefined_sites import FlexibleSites
from gg.modifiers import ClusterRotate, ClusterTranslate, Add, Remove, ModifierAdder, Replace

calc = mace_mp(model="medium-mpa-0",default_dtype="float64",device='cuda')
atoms = read('POSCAR')
atoms.calc = calc

for a in atoms:
    if a.symbol in ['Zn','O','H']:
        a.tag=-1

FS = FlexibleSites(tag=-1)
FS2 = FlexibleSites(constraints=True,com=0.75)

#Add,Remove,Swap O
addO = Add(FS2, "O", surf_coord=[2,3], ads_id = ["O"], surf_sym = ["Cu","Zn"], unique_method = "O")
remO = Remove(FS, "O", max_bond_ratio = 1.2, unique_method = "O")
swapO = ModifierAdder([remO,addO],unique_method = "O")

#Add,Remove,Swap H
addH = Add(FS2, "H", surf_coord=[1,2,3], ads_id = ["H"], surf_sym = ["Cu","O"], unique_method = "H")
remH = Remove(FS, "H", max_bond_ratio = 1.2, unique_method = "H")
swapH = ModifierAdder([remH,addH],unique_method = "H")

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

G = Gcbh(atoms,config_file='input.yaml',restart=True)
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

G.run(steps = 1000)
