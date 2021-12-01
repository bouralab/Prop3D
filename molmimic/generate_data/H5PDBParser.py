# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for PDB files."""


import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use the PDB parser."
    ) from None

from Bio.File import as_handle

from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.parse_pdb_header import _parse_pdb_header_list


# If PDB spec says "COLUMNS 18-20" this means line[17:20]


class H5PDBParser:
    """Parse a PDB file and return a Structure object."""

    def __init__(
        self,
        PERMISSIVE=True,
        get_header=False,
        structure_builder=None,
        QUIET=False,
        is_pqr=False,
        distributed=False
    ):
        """Create a PDBParser object.

        The PDB parser call a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instanciated by the
        PDBParser object itself, but if the user provides his/her own
        StructureBuilder object, the latter is used instead.

        Arguments:
         - PERMISSIVE - Evaluated as a Boolean. If false, exceptions in
           constructing the SMCRA data structure are fatal. If true (DEFAULT),
           the exceptions are caught, but some residues or atoms will be missing.
           THESE EXCEPTIONS ARE DUE TO PROBLEMS IN THE PDB FILE!.
         - get_header - unused argument kept for historical compatibilty.
         - structure_builder - an optional user implemented StructureBuilder class.
         - QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
           the SMCRA data will be suppressed. If false (DEFAULT), they will be shown.
           These warnings might be indicative of problems in the PDB file!
         - is_pqr - Evaluated as a Boolean. Specifies the type of file to be parsed.
           If false (DEFAULT) a .pdb file format is assumed. Set it to true if you
           want to parse a .pqr file instead.
         - distributed - Use disributed h5 reading (from h5serv or HSDS using h5pyd)
           to load PDB info. ENDPOINTS must be set as env vars according to h5pyd.

        """
        # get_header is not used but is left in for API compatibility
        if structure_builder is not None:
            self.structure_builder = structure_builder
        else:
            self.structure_builder = StructureBuilder()
        self.header = None
        self.trailer = None
        self.line_counter = 0
        self.PERMISSIVE = bool(PERMISSIVE)
        self.QUIET = bool(QUIET)
        self.is_pqr = bool(is_pqr)
        self.disributed = disributed

        if self.distributed:
            try:
                import h5pyd
            except ImportError:
                raise RuntimeError("h5pyd not instlled. Please run 'pip install h5pyd' and set correct environment variables.")
                self.h5py = h5pyd
        else:
            self.h5py = h5py

    # Public methods

    def get_structure(self, id, file):
        """Return the structure.

        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the PDB file OR an open filehandle

        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            self.header = None
            self.trailer = None
            # Make a StructureBuilder instance (pass id of structure as parameter)
            self.structure_builder.init_structure(id)

            with self.h5py.File(file, mode='r', use_cache=False) as f:
                try:
                    data = f[id]
                except KeyError:
                    raise ValueError(f"key '{id}' does not exist in {file}")

                try:
                    self.header = data.attrs["header"]
                except KeyError:
                    self.header = {}

                if len(data) == 0:
                    raise ValueError("Empty file.")

                self._parse(lines)

            self.structure_builder.set_header(self.header)
            # Return the Structure instance
            structure = self.structure_builder.get_structure()

        return structure

    def get_header(self):
        """Return the header."""
        return self.header

    def get_trailer(self):
        """Return the trailer."""
        return self.trailer

    # Private methods

    def _parse(self, header_coords_trailer):
        """Parse the PDB file (PRIVATE)."""
        # Extract the header; return the rest of the file
        self.header, coords_trailer = self._get_header(header_coords_trailer)
        # Parse the atomic data; return the PDB file trailer
        self.trailer = self._parse_coordinates(coords_trailer)

    def _parse_coordinates(self, dataset):
        """Parse the atomic data in the PDB file (PRIVATE)."""
        allowed_records = {
            "ATOM  ",
            "HETATM",
            "MODEL ",
            "ENDMDL",
            "TER   ",
            "ANISOU",
            # These are older 2.3 format specs:
            "SIGATM",
            "SIGUIJ",
            # bookkeeping records after coordinates:
            "MASTER",
        }

        local_line_counter = 0
        structure_builder = self.structure_builder


        current_chain_id = None
        current_segid = None
        current_residue_id = None
        current_resname = None

        if isinstance(data, self.h5py.Dataset):
            #one model
            models = [data]
        else:
            models = data.values()

        for model_id, model in enumerate(models):
            # Initialize the Model - there was no explicit MODEL record
            structure_builder.init_model(model_id)
            current_chain_id = None
            current_segid = None
            current_residue_id = None
            current_resname = None


            for record_num, record in enumerate(dataset):
                global_line_counter = f'(Model:{model_id}, Record:{record_num})'
                structure_builder.set_line_counter(global_line_counter)

                fullname = record[1]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name = fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name = split_list[0]
                altloc = record[2]
                resname = record[3]
                chainid = record[4]
                try:
                    serial_number = int(record[0])
                except Exception:
                    serial_number = 0
                resseq = int(record[5].split()[0])  # sequence identifier
                icode = record[6]  # insertion code
                if record_type == "HETATM":  # hetero atom flag
                    if resname == "HOH" or resname == "WAT":
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H"
                else:
                    hetero_flag = " "
                residue_id = (hetero_flag, resseq, icode)
                # atomic coordinates
                coord = record[7:10]

                # occupancy & B factor
                if not self.is_pqr:
                    occupancy = record[10]
                    if occupancy is not None and occupancy < 0:
                        # TODO - Should this be an error in strict mode?
                        # self._handle_PDB_exception("Negative occupancy",
                        #                            global_line_counter)
                        # This uses fixed text so the warning occurs once only:
                        warnings.warn(
                            "Negative occupancy in one or more atoms",
                            PDBConstructionWarning,
                        )
                    bfactor = record[11]

                elif self.is_pqr:
                    # Attempt to parse charge and radius fields
                    pqr_charge = record[10]
                    radius = record[11]

                    if radius is not None and radius < 0:
                        # In permissive mode raise fatal exception.
                        message = "Negative atom radius"
                        self._handle_PDB_exception(message, global_line_counter)
                        radius = None

                segid = record[12]
                element = record[13]
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

                if not self.is_pqr:
                    # init atom with pdb fields
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            bfactor,
                            occupancy,
                            altloc,
                            fullname,
                            serial_number,
                            element,
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif self.is_pqr:
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            pqr_charge,
                            radius,
                            altloc,
                            fullname,
                            serial_number,
                            element,
                            pqr_charge,
                            radius,
                            self.is_pqr,
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

        return []

    def _handle_PDB_exception(self, message, line_counter):
        """Handle exception (PRIVATE).

        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE), or raises it again, this time adding the
        PDB line number to the error message.
        """
        message = "%s at line %i." % (message, line_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn(
                "PDBConstructionException: %s\n"
                "Exception ignored.\n"
                "Some atoms or residues may be missing in the data structure."
                % message,
                PDBConstructionWarning,
            )
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message) from None
