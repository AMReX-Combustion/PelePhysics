import sys
import os

pele_physics_home = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),os.pardir,os.pardir, os.pardir, os.pardir, os.pardir))
dv_dir=os.path.join( pele_physics_home,"Support/Fuego/Pythia")
export_root = os.path.join(dv_dir,"pythia-0.4")
sys.path.append(os.path.join(export_root, "packages/fuego"))
sys.path.append(os.path.join(export_root, "packages/pyre"))
sys.path.append(os.path.join(export_root, "packages/journal"))
sys.path.append(os.path.join(export_root, "packages/weaver"))
    
