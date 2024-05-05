import os
import sys

# no need to do separate calculation for pure metals when making alloy-including database,

add_pure_metals, alloy = True, False
env = "grand"  # "ito" or "grand"

#adsorbates = ["CO", "CH3", "NO", "N2", "H2"]
adsorbates = ["CO"]

adsorbates = "-".join(adsorbates)

dict = {
  "0":  {"element": "Rh", "lattice": "fcc"},
  "1":  {"element": "Pd", "lattice": "fcc"},
  "2":  {"element": "Ir", "lattice": "fcc"},
  "3":  {"element": "Pt", "lattice": "fcc"},
  "4":  {"element": "Cu", "lattice": "fcc"},
  "5":  {"element": "Ag", "lattice": "fcc"},
  "6":  {"element": "Au", "lattice": "fcc"},
  "7":  {"element": "Ru", "lattice": "fcc"},
  "8":  {"element": "Os", "lattice": "fcc"},
  "9":  {"element": "Cd", "lattice": "fcc"},
  "10": {"element": "Hg", "lattice": "fcc"}
  }

if alloy:
    for i in range(0, len(dict)):
        elem1 = dict[str(i)]["element"]
        for j in range(i+1, len(dict)):
            elem2 = dict[str(j)]["element"]

            for ads in adsorbates:
                # --- change by 50% --> gives 21 systems w.o. pure metals
                for comp1 in range(50, 100, 50):  # start, end, diff

                # --- change by 25% --> gives 63 systems w.o. pure metals
                #for comp1 in range(25,100,25): # start, end, diff

                # --- change by 10%
                #for comp1 in range(10,50,10):  # start, end, diff # part1
                #for comp1 in range(50,100,10): # start, end, diff # part2

                # --- change by 5%
                #for comp1 in range(25,50,5):  # start, end, diff # part2
                #for comp1 in range(50,75,5):  # start, end, diff # part3
                #for comp1 in range(75,100,5): # start, end, diff # part4
                #for comp1 in range(5, 25, 5):   # start, end, diff # part1

                    if env == "ito":
                        cmd = "pjsub -x INP1={0} -x INP2={1} -x INP3={2} -x INP4={3} run_ito.sh".format(elem1, elem2, comp1, ads)
                    elif env == "grand":
                        cmd = "pjsub -x INP1={0} -x INP2={1} -x INP3={2} -x INP4={3} run_grand.sh".format(elem1, elem2, comp1, ads)

                    os.system(cmd)

# do pure metal if the loop is first_time

if add_pure_metals:
    for i in range(0, len(dict)):
        elem = dict[str(i)]["element"]
        if env == "ito":
            cmd = 'pjsub run_ito.sh   -x INP1={0} -x INP2="{1}"'.format(elem, adsorbates)
        elif env == "grand":
            cmd = 'pjsub run_grand.sh -x INP1={0} -x INP2="{1}"'.format(elem, adsorbates)
        os.system(cmd)

