from ase.calculators.vasp import Vasp
from ase import Atoms
from ase.io import read
from ase.visualize import view
from ase.build import surface, add_adsorbate
import yaml
import numpy as np

vacuum = 6.0
bulk = read("test.cif")
surf = surface(lattice=bulk, indices=[1, 0, 0], layers=2,
               vacuum=vacuum, periodic=True)
surf = surf*[2, 2, 1]
surf.translate([0, 0, -vacuum+0.1])

adsorbate = Atoms("N")
adsorbate.pbc = True
adsorbate.cell = [10.0, 10.0, 10.0]
adsorbate.center()

surf_ads = surf.copy()
add_adsorbate(slab=surf_ads, adsorbate=adsorbate, height=2.0,
              position=[0, 0], offset=[0.25, 0.25])

with open("vasp_setting.yaml", "r") as f:
    vasp_setting = yaml.safe_load(f)

num_elec_tot = np.zeros(4)  # s, p, d, f, valence only
num_elec = {"O":  [2, 4, 0, 0],
            "Ba": [2, 0, 0, 0],
            "Zr": [2, 0, 2, 0]}

for iatom in surf:
    num_elec_tot += np.array(num_elec[iatom.symbol])

vasp_calc = Vasp(xc=vasp_setting["xc"],
                 ibrion=vasp_setting["ibrion"],
                 potim=vasp_setting["potim"],
                 pp="pbe",
                 lorbit=10,
                 nsw=vasp_setting["nsw"],
                 nelm=vasp_setting["nelm"],
                 kpts=vasp_setting["kpts"])

surf.calc = vasp_calc
adsorbate.calc = vasp_calc
surf_ads.calc = vasp_calc

ads_energy  = adsorbate.get_potential_energy()
surf_energy = surf.get_potential_energy()
surf_ads_energy = surf_ads.get_potential_energy()
adsorption_energy = surf_ads_energy - (ads_energy + surf_energy)

print(adsorption_energy)
