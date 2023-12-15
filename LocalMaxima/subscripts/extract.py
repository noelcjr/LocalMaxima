import argparse
import sys
import LocalMaxima.structure.structure3D as strct

def register_parser(subparsers):
    parser = subparsers.add_parser('extract', usage=usage(), 
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", help="""Initial structire
                        input file in PDF or CIF format.""", 
                        required=True)
    parser.add_argument("--expression", type=str,help="""A Structural expresion that 
                        must be a subset of --structure. If the structural 
                        expression includes elements not found in the --structure, 
                        the program will not run. The structural expression 
                        must be in \'\"\' quotation marks.""", required=True)
    parser.add_argument("--output", metavar="FILE", type=str, default="OUT_2_TERMINAL",
                        help="""Path filename for output structure file in PDB format. If
                        no --output is provided, structure will be output to the
                        standard output in the shell prompt.""",required=False)
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.structure)
    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        structure.get_se_output(options.expression, options.output)
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return """This command extracts parts of the structure into separate PDB 
           files defined by a Structural Expression."""

def usage():
    return """\nmodify extract --structure 1BRS.pdb --expression \
"[0]c[A]r[3:5]a[*]"\nmodify extract --structure 1BRS.pdb --expression \
"m[0]c[A]r[3:5]a[*]" --out 1BRS_A.pdb\n"""

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
