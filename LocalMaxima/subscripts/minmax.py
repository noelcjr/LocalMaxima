import argparse
import sys
import LocalMaxima.structure.structure3D as strct

def register_parser(subparsers):
    parser = subparsers.add_parser('minmax', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", help="""Input structure 
                        file.""", required=True)
    parser.add_argument("--expression", metavar="SE", type=str, 
                        default="m[*]c[*]r[*]a[*]", help="""A Structural 
                        expression (SE) that must be a subset of --structure 
                        is an option. Results will correspond to the subset 
                        structure, but the structure will remain intact. If 
                        the structural expression includes elements not found 
                        in the --structure, the program will not run. The 
                        structural expression must be in \" quotation marks."""
                        , required=False)
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.structure)
    if ((structure.file_sufix.lower() == 'pdb') or 
       (structure.file_sufix.lower() == 'cif')):
        print(options.expression)
        structure.min_max_se(options.expression)
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return '''This command gives the maximum and minim XYZ coordinate values of
     all atoms in the structure. This information is used to add a water box 
    around the protein that extends a defined number of agnstroms from the 
    protein. This values are used to determine the right box dimensions for
    minimum-image conversion.'''

def usage():
    return '\npdb_cif.py minmax --input 1BRS.pdb\n or type strctr minmax --help'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
