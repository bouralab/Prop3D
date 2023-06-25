"""Useful tables of protein properties
"""

hydrophobicity_scales = {
    "kd": {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
           'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
           'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
           'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 },
    "biological": {
           "A": -0.11,"C":  0.13,"D": -3.49,"E": -2.68,"F": 0.32,
           "G": -0.74,"H": -2.06,"I":  0.60,"K": -2.71,"L": 0.55,
           "M":  0.10,"N": -2.05,"P": -2.23,"Q": -2.36,"R": -2.58,
           "S": -0.84,"T": -0.52,"V":  0.31,"W": -0.30,"Y": -0.68},
    "octanal":{
           "A": -0.50, "C":  0.02, "D": -3.64, "E": -3.63,
           "F":  1.71, "G": -1.15, "H": -0.11, "I":  1.12,
           "K": -2.80, "L":  1.25, "M":  0.67, "N": -0.85,
           "P": -0.14, "Q": -0.77, "R": -1.81, "S": -0.46,
           "T": -0.25, "V":  0.46, "W":  2.09, "Y":  0.71,}
    }

#vdw_radii = {'Ru': 1.2, 'Re': 4.3, 'Ra': 2.57, 'Rb': 2.65, 'Rn': 2.3, 'Rh': 1.22, 'Be': 0.63, 'Ba': 2.41, 'Bi': 1.73, 'Bk': 1.64, 'Br': 1.85, 'D': 1.2, 'H': 1.2, 'P': 1.9, 'Os': 1.58, 'Es': 1.62, 'Hg': 1.55, 'Ge': 1.48, 'Gd': 1.69, 'Ga': 1.87, 'Pr': 1.62, 'Pt': 1.72, 'Pu': 1.67, 'C': 1.775, 'Pb': 2.02, 'Pa': 1.6, 'Pd': 1.63, 'Cd': 1.58, 'Po': 1.21, 'Pm': 1.76, 'Ho': 1.61, 'Hf': 1.4, 'K': 2.75, 'He': 1.4, 'Md': 1.6, 'Mg': 1.73, 'Mo': 1.75, 'Mn': 1.19, 'O': 1.45, 'S': 1.8, 'W': 1.26, 'Zn': 1.39, 'Eu': 1.96, 'Zr': 1.42, 'Er': 1.59, 'Ni': 1.63, 'No': 1.59, 'Na': 2.27, 'Nb': 1.33, 'Nd': 1.79, 'Ne': 1.54, 'Np': 1.71, 'Fr': 3.24, 'Fe': 1.26, 'Fm': 1.61, 'B': 1.75, 'F': 1.47, 'Sr': 2.02, 'N': 1.5, 'Kr': 2.02, 'Si': 2.1, 'Sn': 2.17, 'Sm': 1.74, 'V': 1.06, 'Sc': 1.32, 'Sb': 1.12, 'Se': 1.9, 'Co': 1.13, 'Cm': 1.65, 'Cl': 1.75, 'Ca': 1.95, 'Cf': 1.63, 'Ce': 1.86, 'Xe': 2.16, 'Lu': 1.53, 'Cs': 3.01, 'Cr': 1.13, 'Cu': 1.4, 'La': 1.83, 'Li': 1.82, 'Tl': 1.96, 'Tm': 1.57, 'Lr': 1.58, 'Th': 1.84, 'Ti': 1.95, 'Te': 1.26, 'Tb': 1.66, 'Tc': 2.0, 'Ta': 1.22, 'Yb': 1.54, 'Dy': 1.63, 'I': 1.98, 'U': 1.75, 'Y': 1.61, 'Ac': 2.12, 'Ag': 1.72, 'Ir': 1.22, 'Am': 1.66, 'Al': 1.5, 'As': 0.83, 'Ar': 1.88, 'Au': 1.66, 'At': 1.12, 'In': 1.93}
vdw_radii = {
    "H" : 1.2,
    "Li" : 1.82,
    "Na" : 2.27,
    "K" : 2.75,
    "C" : 1.7,
    "N" : 1.55,
    "O" : 1.52,
    "F" : 1.47,
    "P" : 1.80,
    "S" : 1.80,
    "Cl" : 1.75,
    "Br" : 1.85,
    "Se" : 1.90,
    "Zn" : 1.39,
    "Cu" : 1.4,
    "Ni" : 1.63,
}

vdw_volume = {
    'H': 7.235,
    'Li': 25.24,
    'Na': 48.972,
    'K': 87.07,
    'C': 20.569,
    'N': 15.591,
    'O': 14.703,
    'F': 13.299,
    'P': 24.417,
    'S': 24.417,
    'Cl': 22.438,
    'Br': 26.508,
    'Se': 28.716,
    'Zn': 11.244,
    'Cu': 11.488,
    'Ni': 18.131
}

vdw_r_volume = {
    1.2: 7.235,
    1.82: 25.24,
    2.27: 48.972,
    2.75: 87.07,
    1.7: 20.569,
    1.55: 15.591,
    1.52: 14.703,
    1.47: 13.299,
    1.8: 24.417,
    1.75: 22.438,
    1.85: 26.508,
    1.9: 28.716,
    1.39: 11.244,
    1.4: 11.488,
    1.63: 18.131
}

vdw_aa_radii = {
    'ALA': 2.52,
    'ARG': 3.60,
    'ASN': 3.0,
    'ASP': 2.94,
    'CYS': 2.92,
    'GLN': 3.25,
    'GLU': 3.21,
    'GLY': 2.25,
    'HIS': 3.308,
    'ILE': 3.22,
    'LEU': 3.20,
    'LYS': 3.42,
    'MET': 3.37,
    'PHE': 3.47,
    'PRO': 2.93,
    'SER': 2.67,
    'THR': 2.90,
    'TRP': 3.73,
    'TYR': 3.55,
    'VAL': 3.01
 }

from numpy import pi
total_surface_areas_with_solvent = {atom:4.*pi*((radius+1.4)**2) for atom, radius \
    in list(vdw_radii.items())}

maxASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0,
          "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "K": 201.0, "L": 236.0,
          "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0,
          "Y": 263.0, "V": 174.0}

import io
import pandas as pd

aa = pd.read_csv(io.StringIO("""aa3\taa1\tmass\tformula\tname
ALA\tALA\tA\t 71.078\t   C3 H5 N O1\t                         Alanine
ARG\tARG\tR\t157.194\t C6 H13 N4 O1\t                        Arginine
ASN\tASN\tN\t114.103\t  C4 H6 N2 O2\t                      Asparagine
ASP\tASP\tD\t114.079\t   C4 H4 N O3\t                   Aspartic Acid
CYS\tCYS\tC\t103.143\t C3 H5 N O1 S\t                         Cystein
GLN\tGLN\tQ\t117.126\t  C4 H9 N2 O2\t                       Glutamine
GLU\tGLU\tE\t128.106\t   C5 H6 N O3\t                   Glutamic Acid
GLY\tGLY\tG\t 57.051\t   C2 H3 N O1\t                         Glycine
HIS\tHIS\tH\t137.139\t  C6 H7 N3 O1\t                       Histidine
ILE\tILE\tI\t113.158\t  C6 H11 N O1\t                      Isoleucine
LEU\tLEU\tL\t113.158\t  C6 H11 N O1\t                         Leucine
LYS\tLYS\tK\t129.180\t C6 H13 N2 O1\t                          Lysine
MET\tMET\tM\t131.196\t C5 H9 N O1 S\t                      Methionine
PHE\tPHE\tF\t147.174\t   C9 H9 N O1\t                   Phenylalanine
PRO\tPRO\tP\t 97.115\t   C5 H7 N O1\t                         Proline
SER\tSER\tS\t 87.077\t   C3 H5 N O2\t                          Serine
THR\tTHR\tT\t101.104\t   C4 H7 N O2\t                       Threonine
TRP\tTRP\tW\t186.210\tC11 H10 N2 O1\t                      Tryptophan
TYR\tTYR\tY\t163.173\t   C9 H9 N O2\t                        Tyrosine
VAL\tVAL\tV\t 99.131\t   C5 H9 N O1\t                          Valine
ABA\tABA\tX\t 85.104\t  C4 H7 N1 O1\t         alpha-aminobutyric acid
ASH\tASH\tD\t115.087\t   C4 H5 N O3\t           Aspartic acid Neutral
CIR\tCIR\tR\t157.170\t C6 H11 N3 O2\t                      citrulline
CME\tCME\tC\t179.260\tC5 H9 N O2 S2\ts,s-(2-hydroxyethyl)thiocysteine
CMT\tCMT\tC\t115.154\t C4 H5 N O1 S\t                o-methylcysteine
CSD\tCSD\tC\t134.134\t C3 H4 N O3 S\t         s-cysteinesulfinic acid
CSO\tCSO\tC\t119.142\t C3 H5 N O2 S\t               s-hydroxycysteine
CSW\tCSW\tC\t135.142\t C3 H5 N O3 S\t              cysteine-s-dioxide
CSX\tCSX\tC\t119.142\t C3 H5 N O2 S\t                  s-oxy cysteine
CYM\tCYM\tC\t102.135\t C3 H4 N O1 S\t                Cystein Negative
CYX\tCYX\tC\t102.135\t C3 H4 N O1 S\t                  Cystein SSbond
DDE\tDDE\tH\t280.346\tC13 H22 N5 O2\t                     diphthamide
GLH\tGLH\tG\t129.114\t   C5 H7 N O3\t          Glutatmic acid Neutral
HID\tHID\tH\t137.139\t  C6 H7 N3 O1\t                       Histidine
HIE\tHIE\tH\t137.139\t  C6 H7 N3 O1\t                       Histidine
HIP\tHIP\tH\t138.147\t  C6 H8 N3 O1\t              Histidine Positive
HSD\tHSD\tH\t137.139\t  C6 H7 N3 O1\t                       Histidine
HSE\tHSE\tH\t137.139\t  C6 H7 N3 O1\t                       Histidine
HSP\tHSP\tH\t138.147\t  C6 H8 N3 O1\t              Histidine Positive
IAS\tIAS\tD\t115.087\t   C4 H5 N O3\t                   beta-aspartyl
KCX\tKCX\tK\t172.182\t C7 H12 N2 O3\t       lysine nz-carboxylic acid
LYN\tLYN\tK\t129.180\t C6 H13 N2 O1\t                  Lysine Neutral
MHO\tMHO\tM\t147.195\t C5 H9 N O2 S\t                 s-oxymethionine
MLY\tMLY\tK\t156.225\t C8 H16 N2 O1\t               n-dimethyl-lysine
MSE\tMSE\tM\t178.091\tC5 H9 N O1 SE\t                selenomethionine
OCS\tOCS\tC\t151.141\t C3 H5 N O4 S\t           cysteinesulfonic acid
PFF\tPFF\tF\t165.164\t C9 H8 F N O1\t        4-fluoro-l-phenylalanine
PTR\tPTR\tY\t243.153\tC9 H10 N O5 P\t               o-phosphotyrosine
SEP\tSEP\tS\t167.057\t C3 H6 N O5 P\t                   phosphoserine
TPO\tTPO\tT\t181.084\t C4 H8 N O5 P\t                phosphothreonine
"""), sep="\t")

def three_to_one(aa_name):
    try:
        return aa.loc[aa_name]["aa1"]
    except KeyError:
        return "X"

def to_int(row):
    out = {}
    for name, el in row.items():
        if el is None:
            out[name] = 0
        else:
            if len(el)>1:
                el = "".join([e for e in el if e.isdigit()])
                try:
                    out[name] = int(el)
                except ValueError:
                    out[name] =  1
            else:
                out[name] =  1
    return pd.Series(out)

formula_count = aa["formula"].str.strip().str.split(" ", expand=True).rename(columns=
    {i:e for i, e in enumerate("CHNOS")}).apply(to_int, axis=1)

def atoms_to_aa(atoms, raise_unknown=True):
    try:
        return _atoms_to_aa[atoms]
    except KeyError:
        pass
    elems = {e:0 for e in "CHNOS"}
    for a in atoms:
        elem = a[0]
        try:
            elems[elem] += 1
        except KeyError:
            pass

    result = formula_count[(formula_count.C==elems["C"])&(formula_count.H==elems["H"])&\
                           (formula_count.N==elems["N"])&(formula_count.O==elems["O"])&\
                           (formula_count.S==elems["S"])]

    if len(result)>0:
        if len(result)==1:
            return result.iloc[0].name
        else:
            return max(_aa_to_atoms.items(), key=lambda x: len(x[1].intersection(atoms)))[0]

    if raise_unknown:
        raise KeyError("Cannot find amino acid for specified atoms")

    return "UNK"



_atoms_to_aa = {frozenset({'C',
            'CA',
            'CB',
            'HA',
            'HB1',
            'HB2',
            'HB3',
            'HN',
            'N',
            'O'}): 'ALA',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CG',
            'CZ',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HD2',
            'HE',
            'HG1',
            'HG2',
            'HH11',
            'HH12',
            'HH21',
            'HH22',
            'HN',
            'N',
            'NE',
            'NH1',
            'NH2',
            'O'}): 'ARG',
 frozenset({'C',
            'CA',
            'CB',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HD21',
            'HD22',
            'HN',
            'N',
            'ND2',
            'O',
            'OD1'}): 'ASN',
 frozenset({'C',
            'CA',
            'CB',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HN',
            'N',
            'O',
            'OD1',
            'OD2'}): 'ASP',
 frozenset({'C',
            'CA',
            'CB',
            'HA',
            'HB1',
            'HB2',
            'HG1',
            'HN',
            'N',
            'O',
            'SG'}): 'CYS',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HE21',
            'HE22',
            'HG1',
            'HG2',
            'HN',
            'N',
            'NE2',
            'O',
            'OE1'}): 'GLN',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HG1',
            'HG2',
            'HN',
            'N',
            'O',
            'OE1',
            'OE2'}): 'GLU',
 frozenset({'C', 'CA', 'HA1', 'HA2', 'HN', 'N', 'O'}): 'GLY',
 frozenset({'C',
            'CA',
            'CB',
            'CD2',
            'CE1',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HD2',
            'HE1',
            'HE2',
            'HN',
            'N',
            'ND1',
            'NE2',
            'O'}): 'HSE',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CG1',
            'CG2',
            'HA',
            'HB',
            'HD1',
            'HD2',
            'HD3',
            'HG11',
            'HG12',
            'HG21',
            'HG22',
            'HG23',
            'HN',
            'N',
            'O'}): 'ILE',
 frozenset({'C',
            'CA',
            'CB',
            'CD1',
            'CD2',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HD11',
            'HD12',
            'HD13',
            'HD21',
            'HD22',
            'HD23',
            'HG',
            'HN',
            'N',
            'O'}): 'LEU',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CE',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HD2',
            'HE1',
            'HE2',
            'HG1',
            'HG2',
            'HN',
            'HZ1',
            'HZ2',
            'HZ3',
            'N',
            'NZ',
            'O'}): 'LYS',
 frozenset({'C',
            'CA',
            'CB',
            'CE',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HE1',
            'HE2',
            'HE3',
            'HG1',
            'HG2',
            'HN',
            'N',
            'O',
            'SD'}): 'MET',
 frozenset({'C',
            'CA',
            'CB',
            'CD1',
            'CD2',
            'CE1',
            'CE2',
            'CG',
            'CZ',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HD2',
            'HE1',
            'HE2',
            'HN',
            'HZ',
            'N',
            'O'}): 'PHE',
 frozenset({'C',
            'CA',
            'CB',
            'CD',
            'CG',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HD2',
            'HG1',
            'HG2',
            'N',
            'O'}): 'PRO',
 frozenset({'C',
            'CA',
            'CB',
            'HA',
            'HB1',
            'HB2',
            'HG1',
            'HN',
            'N',
            'O',
            'OG'}): 'SER',
 frozenset({'C',
            'CA',
            'CB',
            'CG2',
            'HA',
            'HB',
            'HG1',
            'HG21',
            'HG22',
            'HG23',
            'HN',
            'N',
            'O',
            'OG1'}): 'THR',
 frozenset({'CE1',
            'CE2',
            'CZ',
            'H2',
            'HE1',
            'HE2',
            'O2',
            'O3',
            'O4',
            'OH',
            'P1'}): 'TP1',
 frozenset({'C',
            'CA',
            'CB',
            'CD1',
            'CD2',
            'CE2',
            'CE3',
            'CG',
            'CH2',
            'CZ2',
            'CZ3',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HE1',
            'HE3',
            'HH2',
            'HN',
            'HZ2',
            'HZ3',
            'N',
            'NE1',
            'O'}): 'TRP',
 frozenset({'C',
            'CA',
            'CB',
            'CD1',
            'CD2',
            'CE1',
            'CE2',
            'CG',
            'CZ',
            'HA',
            'HB1',
            'HB2',
            'HD1',
            'HD2',
            'HE1',
            'HE2',
            'HH',
            'HN',
            'N',
            'O',
            'OH'}): 'TYR',
 frozenset({'C',
            'CA',
            'CB',
            'CG1',
            'CG2',
            'HA',
            'HB',
            'HG11',
            'HG12',
            'HG13',
            'HG21',
            'HG22',
            'HG23',
            'HN',
            'N',
            'O'}): 'VAL'}

_aa_to_atoms = {v:k for k,v in _atoms_to_aa.items()}
