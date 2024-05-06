def lattice_info_guess(bulk):
    from ase import Atoms
    #
    # dict for already-known cases
    #
    dict = {
        "Os": {"lattice": "fcc", "a0": 2.734},
        "Ru": {"lattice": "fcc", "a0": 2.706},
        "Cd": {"lattice": "fcc", "a0": 2.979},
        "Hg": {"lattice": "fcc", "a0": 3.005},
        "Rh": {"lattice": "fcc", "a0": 3.803},
        "Pd": {"lattice": "fcc", "a0": 3.891},
        "Ir": {"lattice": "fcc", "a0": 3.839},
        "Pt": {"lattice": "fcc", "a0": 3.924},
        "Cu": {"lattice": "fcc", "a0": 3.615},
        "Ag": {"lattice": "fcc", "a0": 4.085},
        "Au": {"lattice": "fcc", "a0": 4.078},
    }
    element = bulk.get_chemical_symbols()[0]
    if element in dict:
        lattice = dict[element]["lattice"]
        a0 = dict[element]["a0"]
        print("I know this. element = {0:s}, a0 = {1:4.2f}".format(element, a0))
    else:
        print("I don't know this.")
        lattice = "fcc"
        a0 = 4.0

    return lattice, a0


def make_bulk(element1, element2=None, comp1=100, lattice="fcc", a0=4.0, repeat=1):
    from ase import Atoms
    from ase.build import bulk
    import numpy as np
    from ase.visualize import view

    # to fractional
    if element2 is not None:
        comp2  = 100.0 - comp1
        comp2 /= 100.0

    # form bulk first
    bulk = bulk(element1, lattice, a=a0, cubic=True)
    bulk = bulk.repeat(repeat)  # for consistency with alloys

    if element2 is not None:

        # make alloy if needed
        natom_tot = len(bulk.get_atomic_numbers())
        natom2    = int(comp2 * natom_tot)

        # replace atoms in bulk randomly
        list = bulk.get_chemical_symbols()

        np.random.seed(1)# set random seed for reproducability
        #np.random.seed()

        replace_list = np.random.choice(len(list), size=natom2, replace=False)
        for i in replace_list:
            list[i] = element2

        bulk.set_chemical_symbols(list)

    return bulk


def optimize_lattice_constant(bulk, lattice="fcc", a0=4.0, xc="PBEsol", isif=7,
                              encut=400, ediff=1.0e-5, ediffg=1.0e-6, nsw=100, npar=1, nsim=1):
    """
    function to do bulk optimization
    """
    from ase import Atoms
    from ase.calculators.vasp import Vasp
    import copy

    bulk_copy = bulk.copy()

    # compuational condition for Vasp
    prec   = "normal"
    potim  = 0.1
    nelmin = 5
    kpts   = [3, 3, 3]
    gamma  = True

    xc = xc.lower()
    if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
        pp = "pbe"
    elif xc == "pw91":
        pp = "pw91"
    elif xc == "lda":
        pp = "lda"
    else:
        print("xc error")

    calc = Vasp(prec=prec, xc=xc, pp=pp, ispin=1,
                ismear=1, sigma=0.2, isif=isif, nelmin=nelmin, encut=encut,
                ibrion=1, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
                kpts=kpts, gamma=gamma, isym=-1, npar=npar, nsim=nsim, lreal=False)
    calc.directory = "lattice_opt"

    bulk_copy.set_calculator(calc)
    bulk_copy.get_potential_energy()

    # note: returning Atoms may result in irregular positions

    return bulk_copy.get_positions(), bulk_copy.cell


def sort_atoms_by_z(atoms):
    from ase import Atoms, Atom
    import numpy as np

    # keep information for original Atoms
    tags = atoms.get_tags()
    pbc  = atoms.get_pbc()
    cell = atoms.get_cell()

    dtype = [("idx", int), ("z", float)]
    zlist = np.array([], dtype=dtype)

    for idx, atom in enumerate(atoms):
        tmp = np.array([(idx, atom.z)], dtype=dtype)
        zlist = np.append(zlist, tmp)

    zlist = np.sort(zlist, order="z")

    newatoms = Atoms()

    for i in zlist:
        idx = i[0]
        newatoms.append(atoms[idx])

    # restore
    newatoms.set_tags(tags)
    newatoms.set_pbc(pbc)
    newatoms.set_cell(cell)

    return newatoms


def make_lobsterin():
    """
    make input file for lobster (lobsterin) for alloy systems
    """
    f = open("lobsterin", "w")
    str = ["COHPstartEnergy -20", "COHPendEnergy 20",
           "basisSet pbeVaspFit2015", "includeOrbitals spd", "cohpbetween atom 1 and atom 2"]
    str = " \n".join(str)
    f.write(str)
    f.close()


def read_cohp(cohpfile="COHPCAR.lobster"):
    f = open(cohpfile, "r")
    f.readline()  # skip first
    ndos = int(f.readline().split()[2])

    [f.readline() for _ in range(2)]

    ene  = []
    cohp = []

    for i in range(ndos):
        tmp = f.readline().split()
        ene.append(float(tmp[0]))
        cohp.append(float(tmp[1]))

    f.close()
    return ene, cohp
