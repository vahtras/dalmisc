#!/usr/bin/env python
import re
import os

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

def main(filename, units=None):
    xyz_string = read_xyzfile(filename)
    mol_string = xyz2mol(xyz_string, units=units)
    with open(rename(filename), "w") as mol_file:
        mol_file.write(mol_string)

def read_xyzfile(filename):
    return open(filename).read()

def xyz2mol(xyzstring, units=None, basis='cc-pVDZ'):
    xyz_lines = xyzstring.split('\n')
    nat = int(xyz_lines[0])
    comment = xyz_lines[1]

    atomtypes = {}
    elements = []
    for atom_line in xyz_lines[2:2+nat]:
        element = extract_element(atom_line)
        if element not in atomtypes:
            elements.append(element)
            atomtypes[element] = {"charge":ELEMENTS.index(element), "atoms":[]}
        atomtypes[element]["atoms"].append(atom_line)

    molstring = """BASIS
%s
%s
%s
Atomtypes=%d"""%(basis, comment, len(comment)*"=", len(atomtypes))
    if units:
        molstring += " Units=%s" % units
    molstring += '\n'

    for element in elements:
        charge = atomtypes[element]["charge"]
        atoms = atomtypes[element]["atoms"]
        molstring += "Charge=%d Atoms=%d\n"%(charge, len(atoms))
        for atom_line in atoms:
            molstring += "%s\n"%atom_line
    return molstring

def extract_element(atom_line):
    element = re.search("(^[A-Z][a-z]?)", atom_line).groups()[0]
    if element not in ELEMENTS:
        raise Exception("No such element:" + element)
    return element


def rename(xyz_filename):
    head, ext = os.path.splitext(xyz_filename)
    if ext != ".xyz": raise Exception("Input file is not xyz type:%s" % ext)
    return head + ".mol"



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('xyz')
    parser.add_argument('--units')
    args = parser.parse_args()

    main(args.xyz, units=args.units)


    

