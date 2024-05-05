from dos_analyzer.vaspdosplus import VaspDosPlus, get_efermi_from_doscar
import os
import glob
import argparse
from tinydb import TinyDB

def get_adsorption_energy(json=None, surface=None, adsorbate=None):
    from ase.db import connect
    db = connect(json)
    data = {}
    system = surface + "_" + adsorbate
    data.update({"system": system})
    id   = db.get(system=system).id
    row  = db.get(id=id)
    eads = row.data.E_ads
    data.update({"E_ads": eads})
    return data

parser = argparse.ArgumentParser()
parser.add_argument("--surf_json", default="surf_x.json", help="name of output json for surface descriptors")
parser.add_argument("--numpeaks",  default=1, type=int, help="number of peaks in surface dos")
args = parser.parse_args()

surf_json = args.surf_json
numpeaks = args.numpeaks
relative_to_fermi = True
with_centers = True

os.system("rm {}".format(surf_json))
db_surf = TinyDB(surf_json)

doscardir = "."
doscars = glob.glob(doscardir + "/" + "DOSCAR*")

surface_doscars = []
for idos in doscars:
    if idos.count("_") == 2:
        surface_doscars.append(idos)

#adsorbates = ["CO", "CH3", "NO", "N2", "H2"]
#adsorbates = ["CO", "CH3", "NO", "NH3"]
adsorbates = ["CO"]

# get descriptor for adsorbates
print("getting descriptors", flush=True)

for adsorbate in adsorbates:
    data = {}
    ads_json = adsorbate + "_x.json"
    os.system("rm {}".format(ads_json))
    db_ads = TinyDB(ads_json)

    dos = VaspDosPlus(doscar="DOSCAR" + "_" + adsorbate)

    dos.system = adsorbate
    dos.numpeaks = 1  # only one peak in adsorbate
    dos.relative_to_fermi = relative_to_fermi

    descriptor = dos.get_descriptors(adsorbate=True, with_centers=with_centers)
    data.update(descriptor)
    
    db_ads.insert(data)

for idoscar in surface_doscars:
    surface = idoscar.split("_")[1] + "_" + idoscar.split("_")[2]
    data = {}

    dos = VaspDosPlus(doscar=idoscar)
    dos.set_database("surf_data.json")
    dos.system = surface
    dos.numpeaks = numpeaks
    dos.relative_to_fermi = relative_to_fermi
    descriptor = dos.get_descriptors(adsorbate=False, with_centers=with_centers)

    data.update(descriptor)
    db_surf.insert(data)

Eads_json = "E_ads.json"
if os.path.exists(Eads_json):
    os.system("rm {}".format(Eads_json))

db_Eads = TinyDB(Eads_json)
for idoscar in surface_doscars:
    for adsorbate in adsorbates:
        surface = idoscar.split("_")[1] + "_" + idoscar.split("_")[2]
        system = surface + "_" + adsorbate
        data = get_adsorption_energy(json="surf_data.json", surface=surface, adsorbate=adsorbate)
        db_Eads.insert(data)
