#!/usr/bin/env python

"""
Extracts one or more chains from a PDB file.

usage: python pdb_selchain.py -<chain> <pdb file>
example: python pdb_selchain.py -A 1CTF.pdb

Author: {0} ({1})

This program is part of the PDB tools distributed with HADDOCK
or with the HADDOCK tutorial. The utilities in this package
can be used to quickly manipulate PDB files, with the benefit
of 'piping' several different commands. This is a rewrite of old
FORTRAN77 code that was taking too much effort to compile. RIP.
"""

import os
import re
import sys

__author__ = "Joao Rodrigues"
__email__ = "j.p.g.l.m.rodrigues@gmail.com"

USAGE = __doc__.format(__author__, __email__)


def check_input(args):
    """
    Checks whether to read from stdin/file and validates user input/options.
    """

    if not len(args):
        # No model, from pipe
        if not sys.stdin.isatty():
            pdbfh = sys.stdin
            model = ' '
        else:
            sys.stderr.write(USAGE)
            sys.exit(1)
    elif len(args) == 1:
        # Model & Pipe _or_ file & no chain
        if re.match('\-[0-9]+', args[0]):
            model = args[0][1:]
            if not sys.stdin.isatty():
                pdbfh = sys.stdin
            else:
                sys.stderr.write(USAGE)
                sys.exit(1)
        else:
            if not os.path.isfile(args[0]):
                sys.stderr.write('File not found: ' + args[0] + '\n')
                sys.stderr.write(USAGE)
                sys.exit(1)
            pdbfh = open(args[0], 'r')
            model = ' '
    elif len(args) == 2:
        # Chain & File
        if not re.match('\-[0-9]+', args[0]):
            sys.stderr.write('Invalid model ID: ' + args[0] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        if not os.path.isfile(args[1]):
            sys.stderr.write('File not found: ' + args[1] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        model = args[0][1:]
        pdbfh = open(args[1], 'r')
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)

    return (model, pdbfh)


def _select_model(fhandle, model_id):
    """Enclosing logic in a function to speed up a bit"""

    coord_re = re.compile('^(ATOM|HETATM|TER)')
    fhandle = fhandle

    parse_model = False
    single_model = False

    for line in fhandle:
        if line.startswith('MODEL') and line.split()[-1].strip() in model_id:
            parse_model = True
            continue
        elif coord_re.match(line) and not parse_model and model_id in ["1", " "]:
            single_model = True
            yield line
            continue

        if coord_re.match(line) and (parse_model or single_model):
            yield line
        elif parse_model and line.startswith('ENDMDL'):
            yield "END"
            break

    if single_model:
        yield "END"

if __name__ == '__main__':
    # Check Input
    model, pdbfh = check_input(sys.argv[1:])

    # Do the job
    new_pdb = _select_model(pdbfh, model)

    try:
        sys.stdout.write(''.join(new_pdb))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)
