import os
from itertools import groupby
from iotbx import file_reader
import iotbx.pdb
from collections import defaultdict
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter as one_three
from mmtbx.geometry import asa
from mmtbx.secondary_structure import build as ssb
from scitbx_array_family_flex_ext import size_t
from scitbx.array_family import flex
import mmtbx.secondary_structure as ss
from cctbx.eltbx.van_der_waals_radii import vdw
import numpy as np

from sklearn.decomposition import PCA

aas = one_three.keys()

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



class Structure(object):
  def __init__(self, pdb_hierarchy, xray_structure, pdb_in, expand_to_p1=False, gridding=None, parent=None, parent_iseqs=None):
    self.pdb_hierarchy = pdb_hierarchy
    self.xray_structure = xray_structure
    self.pdb_in = pdb_in
    self.unit_cell = self.xray_structure.crystal_symmetry().unit_cell()
    self.gridding = gridding or [1, 1, 1]
    if expand_to_p1:
      self.expand_to_p1()
    self.asa = None
    self.ss = None
    self.parent = parent
    self.parent_iseqs = parent_iseqs
    self.children = []
    self.mean_coord = self.pdb_hierarchy.atoms().extract_xyz().mean()
    
  @classmethod
  def from_pdb(cls, pdb, expand_to_p1=False, gridding=None):
    pdb_in, pdb_hierarchy = open_structure(pdb)
    xray_structure = pdb_in.xray_structure_simple()
    return cls(pdb_hierarchy, xray_structure, pdb_in, expand_to_p1, gridding)

  def get_parent(self):
    return self.parent if self.parent else self

  def make_grid(self, gridding):
    self.gridding = gridding

  def extract_peptides(self):
    return self.extract_selection_from_string("pepnames")

  def extract_chain(self, chain):
    return self.extract_selection_from_string("chain {}".format(chain))

  def extract_selection_from_string(self, sel_string, update_i_seq=False):
    sel_cache = self.pdb_hierarchy.atom_selection_cache()
    selection = sel_cache.selection(sel_string)
    return self.extract_selection(selection, update_i_seq=update_i_seq)

  def extract_selection(self, selection, update_i_seq=False):
    sel_hierarchy = self.pdb_hierarchy.select(selection) #Use parent so atom iseqs always match
    sel_xray = self.xray_structure.select(selection)
    sel_parent_iseqs = self.pdb_hierarchy.atoms().extract_i_seq().as_numpy_array().tolist()
    sel_hierarchy.reset_atom_i_seqs() 
    new = Structure(sel_hierarchy, sel_xray, self.pdb_in, gridding=self.gridding, parent=self.get_parent(), parent_iseqs=sel_parent_iseqs)
    self.children.append(new)
    return new

  def expand_to_p1(self):
    """
    """
    self.xray_structure = self.xray_structure.expand_to_p1()
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

  def orient_to_pai(self, flip_axis=(0.2, 0.2, 0.2)):
    #pai = self.xray_structure.principal_axes_of_inertia()
    coords = PCA(n_components = 3).fit_transform(self.pdb_hierarchy.atoms().extract_xyz())
    coords = flip_around_axis(coords, axis=flip_axis)

    self.xray_structure = self.pdb_hierarchy.extract_xray_structure()
    self.xray_structure.replace_sites_cart(flex.vec3_double(coords))
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

  def get_interface(self, other_structure, max_distance=5, residue_level=True):
    interface = defaultdict(list)
    res1 = set()
    res2 = set()
    iseqs1 = []
    iseqs2 = []
    if residue_level:
        get_unit = lambda atom: atom.parent().parent()
    else:
        get_unit = lambda atom: atom

    for atom1 in self.pdb_hierarchy.atoms():
        if atom1.parent().resname not in one_three.keys(): continue
        for atom2 in other_structure.pdb_hierarchy.atoms():
            if atom2.parent().resname not in one_three.keys(): continue
            distance = atom1.distance(atom2)
            if distance <= max_distance:
                residue1 = get_unit(atom1)
                residue2 = get_unit(atom2)
                interface[(residue1, residue2)].append(distance)
                res1.add(residue1.resseq_as_int())
                res2.add(residue2.resseq_as_int())

                new_atoms1 = residue1.atoms().extract_i_seq()
                if len(iseqs1) == 0:
                    iseqs1 = new_atoms1.as_numpy_array().tolist()
                if new_atoms1[0] <= iseqs1[-1] and new_atoms1[-1] > iseqs1[0]:
                    #Already included
                    pass
                if iseqs1[-1] < new_atoms1[0]:
                    iseqs1 += new_atoms1.as_numpy_array().tolist()
                elif iseqs1[0] > new_atoms1[-1]:
                    iseqs1 = new_atoms1.as_numpy_array().tolist()+iseqs1

                new_atoms2 = residue2.atoms().extract_i_seq()
                if len(iseqs2) == 0:
                    iseqs2 = new_atoms2.as_numpy_array().tolist()
                if new_atoms2[0] <= iseqs2[-1] and new_atoms2[-1] > iseqs2[0]:
                    #Already included
                    pass
                if iseqs2[-1] < new_atoms2[0]:
                    iseqs2 += new_atoms2.as_numpy_array().tolist()
                elif iseqs2[0] > new_atoms2[-1]:
                    iseqs2 = new_atoms2.as_numpy_array().tolist()+iseqs2

    #import pdb; pdb.set_trace()

    distances = {residues:min(distances) for residues, distances in interface.iteritems()}
    binding_site1 = self.extract_selection(size_t(iseqs1))
    binding_site2 = other_structure.extract_selection(size_t(iseqs2))
    return binding_site1, binding_site2, distances



  def get_fractional(self, grid):
    """Convert grid point to fractional coordinates

    Parameters:
    ___________
    grid : 3-tuple.
      grid coordinate

    Return:
    _______
    xyz : fractional coordinates
    """
    return [float(grid[j])/self.gridding[j] for j in xrange(3)]

  def get_cartesian(self, grid):
    """Convert grid point to cartesian coordinates

    Parameters:
    ___________
    grid : 3-tuple.
      grid coordinate

    Return:
    _______
    xyz : 3-tuple.
      cartesian coordinates
    """
    frac = self.get_fractional(grid)
    xyz = self.unit_cell.orthogonalize(frac)
    return xyz

  def get_grid(self, cart):
    frac = self.unit_cell.fractionalize(cart)
    result = [f*size for f, size in zip(frac, self.gridding)]
    return result

  def get_accessible_surface_area(self, atom=None, index=None):
    assert [atom, index].count(None) == 1
    if self.asa is None:
        self.asa = asa.calculate(self.pdb_hierarchy.atoms())

    if atom is not None:
        index = atom.i_seq

    try:
        return self.asa.areas[index]
    except IndexError:
        return 0

  def get_hydrophobicity(self, atom, scale="kd"):
    assert scale in hydrophobicity_scales.keys()
    try:
        return hydrophobicity_scales[scale][one_three[atom.parent().resname]]
    except KeyError:
        return 0

  def get_residue(self, atom, one_hot=False):
    """ADD:
    RESIDUE_NAME_IS_HOH
    RESIDUE_NAME_IS_OTHER
    """
    if one_hot:
        residue = [0 for _ in xrange(len(one_three)+1)] #one hot residue
        try:
            residue[aas.index(atom.parent().resname)] = 1
        except ValueError:
            residue[-1] = 1
        return residue
    else:
        return atom.parent().resname

  def get_features_for_atom(self, atom, i_seq=False):
    """Calculate FEATUREs"""
    if i_seq:
        atom = self.pdb_hierarchy.atoms()[atom]

    atom_type             = self.get_atom_type(atom)
    #partial_charge        = None
    element_type          = self.get_element_type(atom)
    #hydroxyl              = None
    #amide                 = None
    #amine                 = None
    #carbonyl              = None
    #ring_system           = None
    #peptide               = None
    vdw_volume            = self.get_vdw(atom)
    charge                = atom.charge_as_int()
    neg_charge            = int(atom.charge_as_int() < 0)
    pos_charge            = int(atom.charge_as_int() >= 0)
    #charge_with_his       = None
    hydrophobicity        = self.get_hydrophobicity(atom)
    #mobility              = None
    solvent_accessibility = self.get_accessible_surface_area(atom)

    residue               = self.get_residue(atom, one_hot=True)
    # residue_class1        = self.get_residue_class2(atom, one_hot=True)
    # residue_class2        = self.get_residue_class2(atom, one_hot=True)
    ss_class             = self.get_ss(atom)

    return atom_type+element_type+[vdw_volume, charge, neg_charge, pos_charge, hydrophobicity, solvent_accessibility]+residue+ss_class

  def get_atom_type(self, atom):
    """
    ATOM_TYPE_IS_C
    ATOM_TYPE_IS_CT
    ATOM_TYPE_IS_CA
    ATOM_TYPE_IS_N
    ATOM_TYPE_IS_N2
    ATOM_TYPE_IS_N3
    ATOM_TYPE_IS_NA
    ATOM_TYPE_IS_O
    ATOM_TYPE_IS_O2
    ATOM_TYPE_IS_OH
    ATOM_TYPE_IS_S
    ATOM_TYPE_IS_SH
    ATOM_TYPE_IS_OTHER"""
    atom_types = ["C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH"]
    return [int(atom.name.strip() == a) for a in atom_types]+[int(atom.name.strip() not in atom_types)]

  def get_element_type(self, atom):
    """ELEMENT_IS_ANY
    ELEMENT_IS_C
    ELEMENT_IS_N
    ELEMENT_IS_O
    ELEMENT_IS_S
    ELEMENT_IS_OTHER"""
    elems = "CNOS"
    return [int(atom.element.strip() == e) for e in elems]+[int(atom.element.strip() not in elems)]

  def get_residue_class1(self, atom):
    """RESIDUE_CLASS1_IS_HYDROPHOBIC
    RESIDUE_CLASS1_IS_CHARGED
    RESIDUE_CLASS1_IS_POLAR
    RESIDUE_CLASS1_IS_UNKNOWN"""
    pass

  def get_residue_class2(self, atom):
    """RESIDUE_CLASS2_IS_NONPOLAR
    RESIDUE_CLASS2_IS_POLAR
    RESIDUE_CLASS2_IS_BASIC
    RESIDUE_CLASS2_IS_ACIDIC
    RESIDUE_CLASS2_IS_UNKNOWN"""
    pass

  def get_ss(self, atom):
    """Returns a 3-vector (is_h, is_b, is_unk)

    TODO: Update for FEATURE:
      SECONDARY_STRUCTURE1_IS_3HELIX
      SECONDARY_STRUCTURE1_IS_4HELIX
      SECONDARY_STRUCTURE1_IS_5HELIX
      SECONDARY_STRUCTURE1_IS_BRIDGE
      SECONDARY_STRUCTURE1_IS_STRAND
      SECONDARY_STRUCTURE1_IS_TURN
      SECONDARY_STRUCTURE1_IS_BEND
      SECONDARY_STRUCTURE1_IS_COIL
      SECONDARY_STRUCTURE1_IS_HET
      SECONDARY_STRUCTURE1_IS_UNKNOWN
      SECONDARY_STRUCTURE2_IS_HELIX
      SECONDARY_STRUCTURE2_IS_BETA
      SECONDARY_STRUCTURE2_IS_COIL
      SECONDARY_STRUCTURE2_IS_HET
      SECONDARY_STRUCTURE2_IS_UNKNOWN
    """
    if self.ss is None:
      #sec_str_from_pdb_file = self.pdb_in.input.extract_secondary_structure()
      sec_str_from_pdb_file = None
      ssm = ss.manager(pdb_hierarchy=self.pdb_hierarchy,
        sec_str_from_pdb_file=sec_str_from_pdb_file)
      self.ss = np.zeros((self.pdb_hierarchy.atoms_size(), 3)).astype(int)
      self.ss[ssm.helix_selection().iselection().as_numpy_array(), 0] = 1
      self.ss[ssm.beta_selection().iselection().as_numpy_array(), 1] = 1

      unknown = ~(ssm.beta_selection() & ssm.helix_selection())
      self.ss[unknown.iselection().as_numpy_array(), 2] = 1

    try:
        return list(self.ss[atom.i_seq])
    except IndexError:
        return [0, 0, 0]

  def get_vdw(self, atom):
    return vdw.table.get(atom.element.strip().title(), 0)

  def get_coords_extended(self, p = 5):
    'Adds the coordinates from interpolation betweens atoms of the backbone'
    # Initialization
    
    C_coords = self.pdb_hierarchy.atoms().extract_xyz()
    new_coords = np.zeros((p*(C_coords.shape[0]-1), C_coords.shape[1]))

    # Computations
    for i in range(1, C_coords.shape[0]):
        for k in range(p, 0, -1):
            new_coords[p*i-k,:] = ((p-k+1)*C_coords[i-1,:] + k*C_coords[i,:])/(p+1)

    # Store results
    return new_coords

  def get_grid_coord(self, atom, vsize = 144, max_radius = 128):
    center = np.array(atom.xyz) - self.mean_coord
    adjusted = center*(vsize/2.-1)/float(max_radius)
    translated = adjusted + (vsize-1)/2. # Translate center
    rounded = translated.astype(int) # Round components
    return rounded

  def get_atoms(self, grids, vsize = 144, max_radius = 128):
    for atom in self.pdb_hierarchy.atoms():
        if self.get_coords(atom) in grids:
            yield atom

def download_pdb(id):
    import urllib
    if not os.path.isfile("{}.pdb".format(id)):
        urllib.urlretrieve("http://files.rcsb.org/download/{}.pdb".format(id), "{}.pdb".format(id))
    return "{}.pdb".format(id)

def open_structure(fname):
    #pdb_in = file_reader.any_file(fname, force_type="pdb").file_object
    pdb_in = iotbx.pdb.input(file_name="/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/pdb/pdb/{}/pdb{}.ent.gz".format(fname[1:3].lower(), fname.lower()))
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    return pdb_in, pdb_hierarchy

def get_interface_by_distance(subunit1_atoms, subunit2_atoms, max_distance=5, residue_level=True):
    interface = defaultdict(list)
    if residue_level:
        get_unit = lambda atom: atom.parent().parent()
    else:
        get_unit = lambda atom: atom

    for atom1 in subunit1_atoms:
        if atom1.parent().resname not in one_three.keys(): continue
        for atom2 in subunit2_atoms:
            if atom2.parent().resname not in one_three.keys(): continue
            distance = atom1.distance(atom2)
            if distance <= max_distance:
                residue1 = get_unit(atom1)
                residue2 = get_unit(atom2)
                interface[(residue1, residue2)].append(distance)

    return {residues:min(distances) for residues, distances in interface.iteritems()}

def get_residue_groups(hierarchy):
    """Must shift down to 0. 4 set for 1kzq"""
    return [rg.resseq_as_int()-3 for rg in hierarchy.residue_groups()]

def get_residue_names(hierarchy, format_string=None):
    if format_string is None:
        format_string = "{index} {name}"
    return [format_string.format(index=rg.resseq_as_int(), name=rg.unique_resnames()[0]) \
        for rg in hierarchy.residue_groups()]

def get_domains_from_scop(allowed_scope_ids=None, allowed_pdbs=None, dir_cla_scopefile=None):
    if dir_cla_scopefile is None:
        dir_cla_scopefile = "/web/public/data/projects/6CysDB/dir.cla.scope.2.06-stable.txt"
    
    toxo = False
    print allowed_pdbs
    if allowed_pdbs == ["1kzq"]:
        import StringIO
        toxo = True
        toxo_scop = StringIO.StringIO("""d1kzqa1   1kzq    A:3-131 b.6.2.1 73374   cl=48724,cf=49502,sf=74877,fa=74878,dm=74879,sp=74880,px=73374
d1kzqa2 1kzq    A:132-255   b.6.2.1 73375   cl=48724,cf=49502,sf=74877,fa=74878,dm=74879,sp=74880,px=73375
d1kzqb1 1kzq    B:3-131 b.6.2.1 73376   cl=48724,cf=49502,sf=74877,fa=74878,dm=74879,sp=74880,px=73376
d1kzqb2 1kzq    B:132-255   b.6.2.1 73377   cl=48724,cf=49502,sf=74877,fa=74878,dm=74879,sp=74880,px=73377
""")

    with open(dir_cla_scopefile) as cla:
        for _ in xrange(4): cla.next()
        for pdb, scop_entries in groupby(cla if not toxo else toxo_scop, key=lambda l: l.split()[1]):
            domain_ranges = {}
            scop_ids = {}
            for entry in scop_entries:
                fields = entry.rstrip().split()
                pdb, scop_id, chain_domain_range = fields[1], fields[3], fields[2]

                if allowed_pdbs and pdb.lower() not in allowed_pdbs: continue
                
                if chain_domain_range.endswith(":"): continue

                try:
                    chain, domain_range = chain_domain_range.split(":")
                    domain_range.replace("-", ":")
                except ValueError:
                    continue
                
                if allowed_scope_ids is not None and not any([scop_id.startswith(i) for i in allowed_scope_ids]):
                    continue

                if not chain in domain_ranges:
                    domain_ranges[chain] = []
                    scop_ids[chain] = []

                domain_ranges[chain].append(domain_range)
                scop_ids[chain].append(scop_id)
            yield pdb, scop_ids, domain_ranges

def split_hierarchy_by_ranges(hierarchy, domain_ranges):
    sel_cache = hierarchy.atom_selection_cache()
    for chain, domain_ranges in domain_ranges.iteritems():
        domains = []
        for domain_range in domain_ranges:
            domain_selection = sel_cache.selection("chain {} and resseq {}".format(chain, domain_range.replace("-", ":")))
            domain_hierarchy = hierarchy.select(domain_selection)
            domains.append(domain_hierarchy)

        yield chain, domains

def split_hierarchy_by_chain(hierarchy, chain):
    sel_cache = hierarchy.atom_selection_cache()
    selection = sel_cache.selection("chain {}".format(chain))
    chain_hierarchy = hierarchy.select(selection)
    return chain_hierarchy

def get_domains_from_hierarchy(pdb, hierarchy, allowed_scope_ids=None):
    for _, scop_ids, domain_ranges in get_domains_from_scop(allowed_pdbs=[pdb],allowed_scope_ids=allowed_scope_ids):
        for chain, domains in split_hierarchy_by_ranges(hierarchy, domain_ranges):
            yield chain, domains

def get_domains_resseqs_from_hierarchy(pdb, hierarchy, allowed_scope_ids=None):
    for _, scop_ids, domain_ranges in get_domains_from_scop(allowed_pdbs=[pdb],allowed_scope_ids=allowed_scope_ids):
        for chain, domains in split_hierarchy_by_ranges(hierarchy, domain_ranges):
            resseqs = [set([r.resseq_as_int() for r in d.residue_groups()]) \
                for d in domains]
            yield chain, resseqs

def get_domain_hierarchies_from_scop(allowed_scope_ids=None, allowed_pdbs=None, dir_cla_scopefile=None):
    for pdb, scop_ids, domain_ranges in get_domains_from_scop(allowed_scope_ids=allowed_scope_ids, allowed_pdbs=allowed_pdbs, dir_cla_scopefile=dir_cla_scopefile):
        pdb_file = download_pdb(pdb)
        if len(domain_ranges) > 0 and all([len(r) > 1 for c, r in domain_ranges.iteritems()]):
            yield pdb, pdb_file, scop_ids, domain_ranges

def get_direction_from_hierarchy(hierarchy, start_axis=None):
    min_axis = hierarchy.extract_xray_structure().principal_axes_of_inertia().min()
    if start_axis is not None:
        return min_axis==start_axis, min_axis
    else:
        return None, min_axis

extended_asas = {}
def get_extended_asa_for_residue(resname):
    # extended_asa = get_extended_asa_for_residue(rg.parent().resname)

     #  #One hot encoding for Levy interface types
      # del_rASA = rASAm-rASAc
      # res_type = [0,0,0,0,0,0]
      # if rASAc < 0.25 and del_rASA == 0:
     #    #interior
     #    res_type[0] = 1
     #  elif rASAc > 0.25 and del_rASA == 0:
     #    #surface
     #    res_type[1] = 1
     #  elif del_rASA > 0 and rASAm < 0.25:
     #    #support
     #    res_type[2] = 1
     #  elif del_rASA > 0 and rASAc > 0.25:
     #    #rim
     #    res_type[3] = 1
     #  elif del_rASA > 0 and rASAm > 0.25 and rASAc < 0.25:
     #    #core
     #    res_type[4] = 1
     #  else:
     #    #unknown
     #    res_type[5] = 1
    if len(resname) == 3:
        resname = one_three[resname]

    if resname in extended_asas:
        return extended_asas[resname]

    #Create initial structure (linear chain of amino acids)
    pdb_hierarchy = ssb.secondary_structure_from_sequence(
        ssb.beta_pdb_str, 
        "G{}G".format(resname))
    pdb_hierarchy.reset_atom_i_seqs()

    for atom in pdb_hierarchy.atoms():
      atom.b = 15.
      atom.occ = 1.0

    extended_asa = asa.calculate(pdb_hierarchy.atoms())

    area = sum([area for atom, area in zip(pdb_hierarchy.atoms(), extended_asa.areas) \
        if int(atom.parent().parent().resseq.lstrip()) == 2])
    extended_asas[resname] = area
    return area

def flip_around_axis(coords, axis = (0.2, 0.2, 0.2)):
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in xrange(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_around_axis(coords, factor = 0, axis = [1,0,0]):
    'Rotation of coords around axis'
    # Adapted from stackoverflow.com/a/6802723
    theta = factor * 2 * np.pi
    axis = np.asarray(axis)
    axis = axis/float(np.sqrt(np.dot(axis, axis)))
    a = np.cos(theta/2.0)
    b, c, d = axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    M = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                  [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                  [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    return np.dot(coords, M)


