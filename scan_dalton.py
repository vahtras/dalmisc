#!/usr/bin/env python

import re
import os
import numpy as np
import scipy.constants
import matplotlib.pyplot as mpl

DIPLEN = ["XDIPLEN", "YDIPLEN", "ZDIPLEN"]
SECMOM = ["XXSECMOM", "XYSECMOM", "XZSECMOM", "YYSECMOM", "YZSECMOM", "ZZSECMOM"]

G_E = scipy.constants.physical_constants["electron g factor"][0]
ALPHA = scipy.constants.alpha

class NotFoundError(Exception):

    def __init__(self, pattern, filename):
        self.pattern = pattern
        self.filename = filename

def get_final_energy(*args):
    energies = []
    for output in args:
        for line in open(output):
            if "Final" in line and "energy:" in line:
                energies.append((last_float(line),))
                break
    return tuple(energies)

# with generators

def get_last_float(pattern, scale=1):
    def get_last_float(*args):
        floats = []
        for output in args:
            found_pattern = False
            for line in open(output):
                if re.search(pattern, line):
                    found_pattern = True
                    floats.append((scale*last_float(line),))
            if not found_pattern:
                raise NotFoundError, (pattern, output)
        return floats
    return get_last_float

get_final_energy = get_last_float("@.*Final.*energy:")
    
def get_total_dipole_moment(*args):
    RN = get_nuclear_dipole_moment(*args)
    Re = get_electronic_dipole_moment(*args)
    Rtot = [n + e for n, e in zip(RN, Re)]
    return tuple(Rtot)

def get_total_quadrupole_moment(*args):
    QN = get_nuclear_quadrupole_moment(*args)
    Qe = get_electronic_quadrupole_moment(*args)
    Qtot = [n + e for n, e in zip(QN, Qe)]
    return tuple(Qtot)

def get_nuclear_dipole_moment(*args):
    R = get_coordinates(*args)
    Z = get_nuclear_charges(*args)
    ndip = [np.dot(z, r) for z, r in zip(Z, R)]
    return tuple(ndip)

def get_coordinates(*args):
    coordinates = []
    for output in args:
        this = []
        lines = open(output)
        for line in lines:
            if "Total number of coordinates:" in line:
                number_of_atoms = int(line.split()[-1])//3
                for a in range(number_of_atoms):
                    words = lines.next().split()
                    this.append(map(float, [words[4], words[7], words[10]]))
                break
        coordinates.append(this)
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
    return line.split()[-1].split(':')[-1]

def get_electronic_dipole_moment(*args):
    diplens = []
    pattern = "DIPLEN  total"
    for output in args:
        diplen = np.zeros(3)
        found_diplen = False
        for line in open(output):
            if pattern in line:
                found_diplen = True
                index, value = read_dipole_component(line)
                diplen[index] = -value
        if not found_diplen:
            raise NotFoundError, (pattern, output)
        diplens.append(diplen)
    return tuple(diplens)

def read_dipole_component(line):
    first = line.split()[0]
    return DIPLEN.index(first), last_float(line)

def read_quadrupole_component(line):
    first = line.split()[0]
    return SECMOM.index(first), last_float(line)

def get_nuclear_quadrupole_moment(*args):
    quadrupoles = []
    for molecular_coordinates in get_coordinates(*args):
        rr = np.array([np.outer(r, r) for r in np.array(molecular_coordinates)])
        quadrupoles.append(
            upper_triangular(np.sum(rr, axis=0))
            )
    return tuple(quadrupoles)

def get_electronic_quadrupole_moment(*args):
    quadrupoles = []
    pattern = "SECMOM total"
    for output in (args):
        quadrupole = np.zeros(6)
        found_quadrupole = False
        for line in open(output):
            if pattern in line:
                found_quadrupole = True
                index, value = read_quadrupole_component(line)
                quadrupole[index] = -value
        if not found_quadrupole:
            raise NotFoundError, (pattern, output)
        quadrupoles.append(quadrupole)
    return tuple(quadrupoles)


def get_polarizability(*args):
    alphas = []
    pattern = "@ QRLRVE:  <<"
    for output in args:
        alpha = np.zeros((3,3))
        found_alpha = False
        for line in open(output):
            if pattern  in line:
                found_alpha = True
                i, j, a = read_pol_component(line)
                alpha[i, j] = a
                alpha[j, i] = a
        if not found_alpha: raise NotFoundError, (pattern, output)
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
    retstr = ""
    for n, (r, q, p, a, b) in enumerate(
        zip(coordinates, charges, dipoles, alphas, betas)
        ):
        retstr += "%d" % (n+1,)
        retstr += (3*fmt) % tuple(r)
        retstr += fmt % q
        retstr += (3*fmt) % tuple(p)
        retstr += (6*fmt) % tuple(upper_triangular(a))
        retstr += (10*fmt) % tuple(upper_triangular(b))
        retstr += "\n"
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

def get_g_rmc(*args, **kwargs):
    rmcs = []
    pattern = "KINENERG active"
    m_s = kwargs.get("M_S", 0.5)
    clebsh_gordan_factor = 1/m_s
    G_FAC = ALPHA**2*G_E/2*clebsh_gordan_factor
    return get_last_float(pattern, G_FAC)(*args, **kwargs)


def get_g_gc1(*args, **kwargs):
    gc1s = []
    M_S = kwargs.get("M_S", 0.5)
    SO_FAC = ALPHA**2/2
    B_x_R_FAC = 0.5
    TRIPLET_FAC = 0.5
    CG_FAC = 1.0/M_S
    MU_B = 0.5
    INTEGRAL_FAC = 2 # Dalton give 1/2 <a|(r.r_c_ - r_cr)/r^3|b>
    #G_FAC = ALPHA**2/(2*M_S)
    G_FAC = INTEGRAL_FAC*SO_FAC*B_x_R_FAC*CG_FAC*TRIPLET_FAC/MU_B
    for output in args:
        gc1 = np.zeros((3, 3))
        for line in open(output):
            gcmatch = re.search("D1-SO (\w)(\w) active", line)
            if gcmatch:
                i, j = ["XYZ".index(k) for k in gcmatch.groups()]
                gc1[i, j] =  last_float(line)
        gc1s.append(gc1*G_FAC)
    return tuple(gc1s)

def get_g_gc2(*args, **kwargs):
    gc2s = []
    M_S = kwargs.get("M_S", 0.5)
    SO_FAC = 1
    B_x_R_FAC = 1
    TRIPLET_FAC = 1
    CG_FAC = 1#1.0/M_S
    MU_B = 1 #0.5
    INTEGRAL_FAC = 1 # Dalton give 1/2 <a|(r.r_c_ - r_cr)/r^3|b>
    #G_FAC = ALPHA**2/(2*M_S)
    G_FAC = INTEGRAL_FAC*SO_FAC*B_x_R_FAC*CG_FAC*TRIPLET_FAC/MU_B
    for output in args:
        gc2 = np.zeros((3, 3))
        for line in open(output):
            gcmatch = re.search("^@G GC2", line)
            if gcmatch:
                values = map(float, line.split()[-7:])
                gc2[0, 0] =  values[0]
                gc2[1, 1] =  values[1]
                gc2[2, 2] =  values[2]
                gc2[0, 1] =  values[3]
                gc2[1, 0] =  values[4]
                gc2[0, 2] =  values[5]
                gc2[2, 0] =  values[6]
        gc2s.append(gc2*G_FAC)
    return tuple(gc2s)

def get_g_oz1(*args, **kwargs):
    oz1s = []
    M_S = kwargs.get("M_S", 0.5)
    SO_FAC = 1
    CG_FAC = 1.0/M_S
    G_FAC = SO_FAC*CG_FAC
    for output in args:
        oz1 = np.zeros((3, 3))
        for line in open(output):
            ozmatch = re.search("<< (\w)ANGMOM  ; (\w)1SPNORB >>", line)
            if ozmatch:
                i, j = ["XYZ".index(k) for k in ozmatch.groups()]
                oz1[i, j] =  last_float(line)
        oz1s.append(oz1*G_FAC)
    return tuple(oz1s)

def get_g_oz2(*args, **kwargs):
    oz2s = []
    M_S = kwargs.get("M_S", 0.5)
    SO_FAC = 1
    CG_FAC = 1.0/M_S
    G_FAC = SO_FAC*CG_FAC
    for output in args:
        oz2 = np.zeros((3, 3))
        for line in open(output):
            ozmatch = re.search("<< (\w)ANGMOM  ; (\w)2SPNORB >>", line)
            if ozmatch:
                i, j = ["XYZ".index(k) for k in ozmatch.groups()]
                oz2[i, j] =  last_float(line)
        oz2s.append(oz2*G_FAC)
    return tuple(oz2s)

def get_orbital_energies(*args, **kwargs):
    orbital_energies = []
    for output in args:
        file_ = open(output)
        for line in file_:
            if "Hartree-Fock orbital energies" in line:
                floats = []
                file_.next() #blank
                line = file_.next().split()[2:]
                while line != []:
                    floats += map(float, line)
                    line = file_.next().split()
                break
        orbital_energies.append(np.array(floats))
    return tuple(orbital_energies)

def get_excitation_energies(*args, **kwargs):
    pattern = r"@ Excitation energy : (.*) au"
    excitation_energies = []
    for output in args:
        file_ = open(output)
        exe = []
        for line in file_:
            wmatch = re.search(pattern, line)
            if wmatch:
                exe.append(float(wmatch.groups(1)[0]))
        if not exe:
            raise NotFoundError(pattern, output)
        excitation_energies.append(np.array(exe))
    return excitation_energies

def get_transition_moments(label, *args, **kwargs):
    pattern = transition_operator_pattern(label)
    moments = []
    for output in args:
        file_ = open(output)
        mom = []
        for line in file_:
            if re.search(pattern, line):
                next_line = file_.next()
                while re.match(XMOM_PATTERN, next_line):
                    mom.append(extract_transition_moment(next_line))
                    next_line = file_.next()
        if not mom: 
            raise NotFoundError, (pattern, output)
        moments.append(np.array(mom))
    return moments

def transition_operator_pattern(label):
    return r'@ Transition operator type: .*%s' % label

XMOM_PATTERN = r'@ STATE NO:.*NT: +(\S+)'

def extract_transition_moment(line):
    match = re.search(XMOM_PATTERN, line)
    return abs(float(match.groups(1)[0]))
    
        
def xyz_to_tuple(string):
    ints = tuple(['xyz'.index(i) for i  in string])
    return ints

def outline(floats, fmt="%10.6f"):
    retstr = ""
    for f in floats:
        retstr += fmt*len(f) % tuple(f) + "\n"
    return retstr

def single_plot(files, numbers):
    mpl.plot(numbers)
    filenames = [os.path.basename(f) for f in files]
    mpl.xticks(range(len(files)), filenames, rotation='vertical')
    mpl.margins(0.2)
    mpl.subplots_adjust(bottom=0.50)
    mpl.show()
    pass


if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser("Scan Dalton output file")
    parser.add_argument('--coordinates', action='store_true')
    parser.add_argument('--energy', action='store_true')
    parser.add_argument('--dipole', action='store_true')
    parser.add_argument('--nuclear-dipole', action='store_true')
    parser.add_argument('--quadrupole', action='store_true')
    parser.add_argument('--nuclear-quadrupole', action='store_true')
    parser.add_argument('--pol', action='store_true')
    parser.add_argument('--diagonal', action='store_true')
    parser.add_argument('--iso', action='store_true')
    parser.add_argument('--upper-triangular', action='store_true')
    parser.add_argument('--hyp', action='store_true')
    parser.add_argument('--select', nargs='+')
    parser.add_argument('--generate-potential', action='store_true')
    parser.add_argument('--cm', action='store_true')
    parser.add_argument('--g-rmc', action='store_true')
    parser.add_argument('--g-gc1', action='store_true')
    parser.add_argument('--g-gc2', action='store_true')
    parser.add_argument('--g-oz1', action='store_true')
    parser.add_argument('--g-oz2', action='store_true')
    parser.add_argument('--excitation-energy', action='store_true')
    parser.add_argument('--transition-moment')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--fmt', default='%10.6f')
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()
    if args.generate_potential:
        print "AU\n%d 1 22 1" % len(args.files)

    if args.coordinates:
        print get_coordinates(*args.files)

    if args.dipole:
        print get_total_dipole_moment(*args.files)

    if args.nuclear_dipole:
        print get_nuclear_dipole_moment(*args.files)

    if args.quadrupole:
        print get_total_quadrupole_moment(*args.files)

    if args.nuclear_quadrupole:
        print get_nuclear_quadrupole_moment(*args.files)

    if args.pol:
        alphas = get_polarizability(*args.files)
        if args.diagonal:
            print outline([alpha.diagonal() for alpha in alphas])
        elif args.iso:
            print outline([(alpha.trace()/3,) for alpha in alphas])
        elif args.upper_triangular:
            print outline([upper_triangular(alpha) for alpha in alphas])
        else:
            print outline([np.ravel(alpha) for alpha in alphas])

    if args.hyp:
        betas = get_hyperpolarizability(*args.files)
        if args.select:
            tuples = [xyz_to_tuple(xyz) for xyz in args.select]
            print outline([[beta[i, j, k] for i, j, k in tuples] for beta in betas])
        elif args.upper_triangular:
            print outline([upper_triangular(beta) for beta in betas])
        else:
            print outline([np.ravel(beta) for beta in betas])

    if args.generate_potential:
        print generate_pot_output(*args.files)

    if args.cm:
        print outline([cm for cm in get_center_of_mass(*args.files)])

    def  blob(files, func, fmt="%10.5f"):
        floats = func(*files)
        retstr = ""
        for f, v in zip(files, floats):
            retstr += "%-40s" % f + fmt*len(v) % tuple(x for x in v) + "\n"
        return retstr

    if args.energy:
        print blob(args.files, get_final_energy, fmt=args.fmt)
        if args.plot:
            floats = np.array([e[0] for e in get_final_energy(*args.files)])
            single_plot(args.files, floats)


    if args.g_rmc:
        rmcs = get_g_rmc(*args.files)
        print outline([np.ravel(rmc) for rmc in rmcs], fmt=args.fmt)

    if args.g_gc1:
        gc1s = get_g_gc1(*args.files)
        print outline([np.ravel(gc1) for gc1 in gc1s], fmt=args.fmt)

    if args.g_gc2:
        gc2s = get_g_gc2(*args.files)
        print outline([np.ravel(gc2) for gc2 in gc2s], fmt=args.fmt)

    if args.g_oz1:
        oz1s = get_g_oz1(*args.files)
        print outline([np.ravel(oz1) for oz1 in oz1s], fmt=args.fmt)

    if args.g_oz2:
        oz2s = get_g_oz2(*args.files)
        print outline([np.ravel(oz2) for oz2 in oz2s], fmt=args.fmt)

    if args.excitation_energy:
        print outline([ws for ws in get_excitation_energies(*args.files)])

    if args.transition_moment:
        try:
            print outline([ws for ws in get_transition_moments(args.transition_moment, *args.files)])
        except NotFoundError as e:
            print e.pattern, e.filename
