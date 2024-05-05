from ase.calculators.vasp import VaspDos
from ase.dft import get_distribution_moment
import numpy as np
import yaml
import matplotlib.pyplot as plt

doscar = "DOSCAR"
dos = VaspDos(doscar=doscar)
energy = dos.energy

with open(doscar, "r") as f:
    line = f.readline()
    num_atoms = int(line.split()[0])

pdos = np.zeros(len(energy))
orbital = 2
for iatom in range(num_atoms):
    pdos += dos.site_dos(atom=iatom, orbital=orbital)

d_center = get_distribution_moment(x=energy, y=pdos, order=1)

with open("descriptors.yaml", "r") as f:
    yaml_file = yaml.safe_load(f)

descriptors = yaml_file["system001"]
descriptors.update({"d_center": d_center})
print(descriptors)