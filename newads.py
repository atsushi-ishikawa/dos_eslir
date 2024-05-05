from ase.calculators.vasp import Vasp
from ase.io import read
from ase.build import surface, add_adsorbate
import yaml

bulk = read("test.cif")
with open("vasp_setting.yaml", "r") as f:
    vasp_setting = yaml.safe_load(f)

print(vasp_setting)
quit()

bulk.calc = Vasp(xc=vasp_setting["xc"],
                 kpts=vasp_setting["kpts"])

energy = bulk.get_potential_energy()

print(energy)