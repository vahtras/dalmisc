import re

ELEMENTS = [
    "X", "H", "He",
    "Li", "Be", 
            "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", 
            "Al", "Si", "p", "S", "Cl", "Ar",
    "K", "Ca", 
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
            "Hf", "Ta", "W" "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"
    ]
numeric=re.compile('\d')


def xyz2mol(xyzstring):
    lines = xyzstring.split('\n')
    nat = int(lines[0])
    comment = lines[1]

    atomtypes = {}
    elements = []
    for atom in lines[2:2+nat]:
        element = extract_element(atom)
        charge = ELEMENTS.index(element)
        if element not in atomtypes:
            elements.append(element)
            atomtypes[element] = {"charge":charge, "atoms":[]}
        atomtypes[element]["atoms"].append(atom)

    molstring = """BASIS
cc-pVDZ
%s
%s
Atomtypes=%d
"""%(comment, len(comment)*"=", len(atomtypes))

    for element in elements:
        charge = atomtypes[element]["charge"]
        atoms = atomtypes[element]["atoms"]
        molstring += "Charge=%d Atoms=%d\n"%(charge, len(atoms))
        for atom in atoms:
            molstring += "%s\n"%atom
    return molstring

def extract_element(atom_line):
    match_digit = re.compile('\d')
    element = re.split(match_digit, atom_line)[0].strip()
    return element

if __name__ == "__main__":
    import sys
    import os
    try:
        xyzfile = sys.argv[1]
        xyzstring = open(xyzfile).read()
    except IndexError, e:
        print "Usage: %s xyzfile"%sys.argv[0]

    n, e = os.path.splitext(xyzfile)
    molfile = open(n+'.mol', 'w')
    molfile.write(xyz2mol(xyzstring))
    molfile.close()
