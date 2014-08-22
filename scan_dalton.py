#!/usr/bin/env python

import numpy as np

DIPLEN = ["XDIPLEN", "YDIPLEN", "ZDIPLEN"]

def get_total_dipole_moment(*args):
    RN = get_nuclear_dipole_moment(*args)
    Re = get_electronic_dipole_moment(*args)
    Rtot = [n + e for n, e in zip(RN, Re)]
    return tuple(Rtot)

def get_nuclear_dipole_moment(*args):
    R = get_coordinates(*args)
    Z = get_nuclear_charges(*args)
    ndip = [np.dot(z, r) for z, r in zip(Z, R)]
    return tuple(ndip)

def get_coordinates(*args):
    coordinates = []
    for output in args:
        print "output", output
        this = []
        lines = open(output)
        for line in lines:
            if "Total number of coordinates:" in line:
                number_of_atoms = int(line.split()[-1])//3
                for a in range(number_of_atoms):
                    words = lines.next().split()
                    this.append(map(float, [words[4], words[7], words[10]]))
                break
        print "this",  this
        coordinates.append(this)
        print tuple(coordinates)
    return tuple(coordinates)

def get_nuclear_charges(*args):
    charge_sets = []
    for output in args:
        charges = []
        lines = open(output)
        for line in lines:
            if "Atomic type no." in line:
                lines.next()
                charge = last_float(lines.next())
                number = last_int(lines.next())
                for n in range(number):
                    charges.append(charge)
        charge_sets.append(np.array(charges))
    return tuple(charge_sets)
            
def last_float(line):
    return float(last(line).replace('D', 'e'))

def last_int(line):
    return int(last(line))

def last(line):
    return line.split()[-1]

def get_electronic_dipole_moment(*args):
    diplens = []
    for output in args:
        diplen = np.zeros(3)
        for line in open(output):
            if "DIPLEN  total" in line:
                index, value = read_dipole_component(line)
                diplen[index] = -value
        diplens.append(diplen)
    return tuple(diplens)

def read_dipole_component(line):
    first = line.split()[0]
    last_double = line.split(':')[-1].replace('D', 'E')
    return DIPLEN.index(first), float(last_double)

def get_polarizability(*args):
    alphas = []
    for output in args:
        alpha = np.zeros((3,3))
        for line in open(output):
            if "@ QRLRVE:  <<" in line:
                i, j, a = read_pol_component(line)
                alpha[i, j] = a
                alpha[j, i] = a
        alphas.append(alpha)
    return tuple(alphas)

def read_pol_component(line):
    words = line.split()
    i = DIPLEN.index(words[3])
    j = DIPLEN.index(words[5])
    a = last_float(line)
    return i, j, a

def get_hyperpolarizability(*args):
    betas = []
    for output in args:
        beta = np.zeros((3, 3, 3))
        for line in open(output):
            if "beta(" in line and line.strip()[-1] != ")":
                i, j, k, b = read_hyp_component(line)
                beta[i, j, k] = b
                beta[j, k, i] = b
                beta[k, i, j] = b
                beta[j, i, k] = b
                beta[k, j, i] = b
                beta[i, k, j] = b
        betas.append(beta)
    return tuple(betas)

def read_hyp_component(line):
    words = line.split()
    beta_component = words[-3]
    i = "XYZ".index(beta_component[5])
    j = "XYZ".index(beta_component[7])
    k = "XYZ".index(beta_component[9])
    try:
        b = last_float(line)
    except ValueError:
        print "Cant do", line
        pass
    return i, j, k, b

def generate_pot_output(*args, **kwargs):
    fmt=kwargs.get('fmt', '%10.5f')
    coordinates = get_center_of_mass(*args)
    #charges = get_charges(*args)
    charges = [0.]*len(args)
    dipoles = get_total_dipole_moment(*args)
    alphas = get_polarizability(*args)
    betas = get_hyperpolarizability(*args)
    for n, (r, q, p, a, b) in enumerate(
        zip(coordinates, charges, dipoles, alphas, betas)
        ):
        retstr = "%d" % (n+1,)
        retstr += (3*fmt) % tuple(r)
        retstr += fmt % q
        retstr += (3*fmt) % tuple(p)
        retstr += (6*fmt) % tuple(upper_triangular(a))
        retstr += (10*fmt) % tuple(upper_triangular(b))
    return retstr

def upper_triangular(mat):
    if len(mat.shape) == 2:
        return np.array([mat[i, j] for i in range(mat.shape[0]) for j in range(i, mat.shape[1])])
    if len(mat.shape) == 3:
        return np.array([mat[i, j, k] for i in range(mat.shape[0]) for j in range(i, mat.shape[1]) for k in range(j, mat.shape[2])])

def get_center_of_mass(*args):
    cms = []
    for output in args:
        for line in open(output):
            if "Center-of-mass" in line:
                cm = [float(x) for x in line.split()[-3:]]
                cms.append(np.array(cm))
                break
    return tuple(cms)
        
def xyz_to_tuple(string):
    ints = tuple(['xyz'.index(i) for i  in string])
    return ints

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser("Scan Dalton output file")
    parser.add_argument('--dipole', action='store_true')
    parser.add_argument('--pol', action='store_true')
    parser.add_argument('--diagonal', action='store_true')
    parser.add_argument('--iso', action='store_true')
    parser.add_argument('--upper-triangular', action='store_true')
    parser.add_argument('--hyp', action='store_true')
    parser.add_argument('--select', nargs='+')
    parser.add_argument('--generate-potential', action='store_true')
    parser.add_argument('--cm', action='store_true')
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()
    if args.generate_potential:
        print "AU\n%d 1 22 1" % len(args.files)
    for n, filename in enumerate(args.files):
        #print filename,
        #print len(filename)*"="
        if args.dipole:
            #print "Nuclear   ", get_nuclear_dipole_moment(open(filename))
            #print "Electronic", get_electronic_dipole_moment(open(filename))
            print "Total     ", get_total_dipole_moment(open(filename))
        if args.pol:
            #print 'alpha',
            alpha = get_polarizability(open(filename))
            if args.diagonal:
                print alpha.diagonal()
            elif args.iso:
                print alpha.trace()/3
            elif args.upper_triangular:
                print upper_triangular(alpha)
            else:
                print alpha
        if args.hyp:
            #print 'beta',
            beta = get_hyperpolarizability(open(filename))
            if args.select:
                tuples = [xyz_to_tuple(xyz) for xyz in args.select]
                print [beta[i, j, k] for i, j, k in tuples]
            if args.upper_triangular:
                print upper_triangular(beta)
            else:
                print beta

        if args.generate_potential:
            print generate_pot_output(n+1, open(filename))

        if args.cm:
            print get_center_of_mass(open(filename))
