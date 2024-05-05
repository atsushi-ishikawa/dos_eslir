from ase.calculators.vasp import Vasp
from ase.io import read
from ase.build import surface, add_adsorbate
import yaml

bulk = read("test.cif")
with open("vasp_setting.yaml", "r") as f:
    vasp_setting = yaml.safe_load(f)

bulk.calc = Vasp(xc=vasp_setting["xc"],
                 ibrion=vasp_setting["ibrion"],
                 potim=vasp_setting["potim"],
                 pp="pbe",
                 lorbit=10,
                 nsw=vasp_setting["nsw"],
                 nelm=vasp_setting["nelm"],
                 kpts=vasp_setting["kpts"])

energy = bulk.get_potential_energy()

print(energy)
