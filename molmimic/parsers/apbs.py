import os

from molmimic.parsers.psize import Psize
from molmimic.parsers.container import Container
from molmimic.parsers.pdb2pqr import Pdb2pqr
from molmimic.util.pdb import get_first_chain, replace_chains

class APBS(Container):
    IMAGE = 'docker://edraizen/apbs:latest'
    LOCAL = ["apbs"]
    PARAMETERS = [("in_file", "path:in", ["{}"])]#, (":out_file", "path:out")]
    RETURN_FILES = True

    def _pdb_to_pqr(self, pdb_file, force_field="amber", **pdb2pqr_kwds):
        remove_pdb = False
        if get_first_chain(pdb_file).isdigit():
            #APBS requires chains A-Z
            pdb_file = replace_chains(pdb_file, new_chain="A")
            remove_pdb = True

        pdb2pqr = Pdb2pqr(work_dir=self.work_dir, job=self.job)
        pdb2pqr_kwds["whitespace"] = pdb2pqr_kwds.get("whitespace", True)
        pdb2pqr_kwds["chain"] = pdb2pqr_kwds.get("chain", True)
        pdb2pqr_kwds["reres"] = pdb2pqr_kwds.get("reres", True)
        pqr_file = pdb2pqr.create_pqr(pdb_file, force_field=force_field,
            **pdb2pqr_kwds)

        print(os.listdir(os.path.dirname(pqr_file)))
        apbs_in_file = f"{os.path.splitext(pdb_file)[0]}.in"
        if not pdb2pqr_kwds.get("apbs_input", False) or not os.path.isfile(apbs_in_file):
            apbs_in_file =  None

        print("apbs_in_file", apbs_in_file, pdb_file, f"{os.path.splitext(pdb_file)[0]}.in")

        if remove_pdb:
            self.files_to_remove.append(pdb_file)

        return pqr_file, apbs_in_file

    def atom_potentials_from_pdb(self, pdb_file, force_field="amber", with_charge=True, **pdb2pqr_kwds):
        pqr_file, apbs_in_file = self._pdb_to_pqr(pdb_file, force_field=force_field, **pdb2pqr_kwds)
        atom_pot_file = self.atom_potentials_from_pqr(pqr_file, apbs_in_file=apbs_in_file)
        return atom_pot_file

    def atom_potentials_from_pqr(self, pqr_file, apbs_in_file=None):
        self.full_pqr_path = pqr_file

        pqr_path = self.format_in_path(None, pqr_file)
        input_file, out_file = self.write_apbs_electrostatics_input(pqr_file, apbs_input=apbs_in_file)
        # else:
        #     #Pdb2pqr was run with '--apbs-input'?
        #     filename_base = os.path.splitext(apbs_in_file)[0]
        #     input_file = f"{filename_base}.in.real"
        #
        #     #Make sure output name is same as pqr file (from MaSIF, computeAPBS)
        #     with open(apbs_in_file) as old, open(f"{filename_base}.in.real", "w") as new:
        #         for line in old:
        #             if "write pot dx pot" in line:
        #                 line = line.replace("write pot dx pot", f"write pot dx {filename_base}")
        #             print(line.rstrip(), file=new)

        #self.files_to_remove.append(input_file)

        try:
            self.running_atom_potentials = True
            atom_pot_file = self(in_file=input_file)
            self.running_atom_potentials = False
        except RuntimeError as e:
            if "Problem opening virtual socket" in str(e):
                try:
                    with open(pqr_file) as f:
                        raise RuntimeError("{} \n PQR file: {}".format(str(e), f.read()))
                except:
                    pass
            raise e
        print(atom_pot_file)
        print(input_file)

        print()
        return atom_pot_file

    def get_atom_potentials_from_pqr(self, pqr_file, apbs_in_file=None, with_charge=True):
        atom_pot_file = self.atom_potentials_from_pqr(pqr_file, apbs_in_file=apbs_in_file)

        if with_charge:
            output = {atom:(charge, potential) for (atom, charge), potential in \
                zip(Pdb2pqr.parse_pqr(pqr_file),
                    self.parse_atom_pot_file(atom_pot_file))}
        else:
            output = {atom:potential for (atom, charge), potential in \
                zip(Pdb2pqr.parse_pqr(pqr_file),
                    self.parse_atom_pot_file(atom_pot_file))}

        self.files_to_remove += [atom_pot_file]

        self.clean()

        return output

    def get_atom_potentials_from_pdb(self, pdb_file, force_field="amber", with_charge=True, **pdb2pqr_kwds):
        pqr_file, apbs_in_file = self._pdb_to_pqr(pdb_file, force_field=force_field, **pdb2pqr_kwds)
        return self.get_atom_potentials_from_pqr(pqr_file, apbs_in_file=apbs_in_file, with_charge=with_charge)

    def local(self, *args, **kwds):
        self.is_local = True
        if hasattr(self, "running_atom_potentials") and self.running_atom_potentials:
            pqr_path = self.format_in_path(None, self.full_pqr_path)
            input_file, out_file = self.write_apbs_electrostatics_input(self.full_pqr_path)
            kwds["in_file"] = input_file
        super(self, APBS).local(*args, **kwds)

    @staticmethod
    def parse_atom_pot_file(atom_pot_file):
        with open(atom_pot_file) as atom_pots:
            #Skip first 4 rows of atompot file
            for _ in range(4):
                next(atom_pots)

            for atom_pot in atom_pots:
                yield float(atom_pot)

    def write_apbs_electrostatics_input(self, pqr_file, ouput_input_prefix=None,
      output_prefix=None, apbs_input=None):
        base = os.path.basename(pqr_file)
        ouput_input_prefix = os.path.join(self.work_dir, ouput_input_prefix or base+".apbs_input")
        output_prefix = os.path.join(self.work_dir, output_prefix or base+".apbs_output")

        if isinstance(apbs_input, str) and os.path.isfile(apbs_input):
            #dx pot has the .dx extension
            pot_path = self.format_out_path(None, output_prefix+".dx")
            with open(apbs_input) as f1, open(ouput_input_prefix, "w") as f2:
                for line in f1:
                    if "mol pqr" in line:
                        line = f"    mol pqr {self.format_in_path(None, pqr_file)}\n"
                    if "write pot" in line:
                        line = f"    write pot {line.strip().split()[2]} {pot_path[:-4]}\n"
                    print(line.rstrip(), file=f2)
                    print(line.rstrip())
        else:
            #Atom pot has the .txt extension
            pot_path = self.format_out_path(None, output_prefix+".txt")
            with open(ouput_input_prefix, "w") as f:
                contents = self.make_apbs_electrostatics_input(pqr_file, pot_path[:-4])
                f.write(contents)
        return ouput_input_prefix, pot_path

    def make_apbs_electrostatics_input(self, pqr_file, output_prefix):
        ps = Psize()
        ps.runPsize(pqr_file)
        pqr_path = self.format_in_path(None, pqr_file)

        cglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getCoarseGridDims())
        fglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getFineGridDims())
        dime = "{:d} {:d} {:d}".format(*ps.getFineGridPoints())

        contents = """read
        mol pqr {pqr}
    end
    elec
        mg-auto # Specify the mode for APBS to run
        dime {dime} # The grid dimensions
        cglen {cglen}
        fglen {fglen}
        cgcent mol 1
        fgcent mol 1
        mol 1 # Perform the calculation on molecule 1
        lpbe # Solve the linearized Poisson-Boltzmann equation
        bcfl sdh # Use all single moments when calculating the potential
        pdie 2.00 # Solute dielectric
        sdie 2.00 #78.54 # Solvent dielectric
        chgm spl2 # Spline-based discretization of the delta functions
        srfm smol # Molecular surface definition
        srad 1.4 # Solvent probe radius (for molecular surface)
        swin 0.3 # Solvent surface spline window (not used here)
        sdens 10.0 # Sphere density for accessibility object
        temp 298.15 # Temperature
        calcenergy total # Calculate energies
        calcforce no # Do not calculate forces
        write atompot flat {out}
    end
    quit
    """.format(
        pqr=pqr_path,
        dime=dime,
        cglen=cglen,
        fglen=fglen,
        out=output_prefix
        )

        return contents
