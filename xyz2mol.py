import re

ELEMENTS = [
    "X", "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne"
    ]
numeric=re.compile('\d')


def xyz2mol(xyzstring):
    lines = xyzstring.split('\n')
    nat = int(lines[0])
    comment = lines[1]

    atomtypes = {}
    elements = []
    for atom in lines[2:2+nat]:
        label = atom.split()[0]
        # Separage chemical symbol from possible numeric label part
        element = re.split(numeric, label)[0]
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
