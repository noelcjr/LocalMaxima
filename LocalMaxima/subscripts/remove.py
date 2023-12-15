import argparse
import sys
import LocalMaxima.structure.add_remove_join as ARJ

def register_parser(subparsers):
    parser = subparsers.add_parser('remove', usage=usage(), 
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--xyz", metavar="FILE", help="""Initial structure 
                        input file in PDF or CIF format.""", 
                        required=True)
    parser.add_argument("--se", type=str,help="""A structural expresion (se) 
                        that must be a subset of --xyz. The structural 
                        expression must be in \'\"\' quotation marks.""", 
                        required=True)
    parser.add_argument("--out", metavar="FILE", type=str, help="""Path and 
                        filename to output modified structure.""",
                        required=False)
    parser.add_argument("--ow", action="store_true", default=False, 
                        help="""If atom removal is to be done on the input --xyz
                        file, this flag must be included to prevent accidental
                        overwrite (ow).""")
    parser.set_defaults(func=run)

def run(options):
    structure = ARJ.mod_structure(options.xyz)
    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        if options.out is None:
            if options.ow:
                structure.remove_from_se(options.se, options.xyz)
            else:
                print("ERROR: Inlcue a path and filename for the structure from which")
                print("       the atoms are removed after --out. If you want to remove")
                print("       atoms from the input structure, include the overwrite")
                print("       flag --ow instead of --out to permanently modify the")
                print("       input file. Program will exit without results")
                sys.exit(1)
        else:
            structure.remove_from_se(options.se, options.out)
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return """This command removes atoms specified by a structural expression.
Results are output to a separate file.\n

WARNING: The format for structures with multiple models will be corrupted if
         the models are not identical in their atom, chain and amino acid
         conformations. This could happen if atoms are removed from only one
         model and not the remaining models. This could cause problems when the 
         modified structure is loaded into other programs.
"""

def usage():
    return """\nmodifr remove --xyz 1BRS.pdb --se "m[0]c[A]r[3:5]a[*]"\n\
modifr remove --xyz 1BRS.pdb --se "m[0]c[A]r[3:5]a[*]" --out 1BRS_1_AD.pdb\n"""

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
