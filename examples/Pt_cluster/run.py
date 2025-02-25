from ase.io import read
from ase.calculators.emt import EMT

atoms = read('POSCAR_Pt_TiO2')
for a in atoms:
    if a.symbol=='Pt':
        a.tag=-1
atoms.calc = EMT() #Add a calculator to the atoms object for geometric optimization

from gg.gcbh import Gcbh
G = Gcbh(atoms,config_file='input.yaml')

from gg.predefined_sites import FlexibleSites
FS = FlexibleSites(tag=True)

from gg.modifiers import ClusterRotate, ClusterTranslate, Add

clstr_rot = ClusterRotate(FS,contact_error=0.25,rotate_vector=(0,0,1))
clstr_trans = ClusterTranslate(FS,contact_error=0.2)

G.add_modifier(clstr_trans,'Cluster_Translate')
G.add_modifier(clstr_rot,'Cluster_Rotate')

G.run(steps = 25)
