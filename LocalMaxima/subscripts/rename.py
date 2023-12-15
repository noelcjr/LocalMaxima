import argparse
import sys
import LocalMaxima.structure.structure3D as strct
import LocalMaxima.structure.structural_expression as SE
import LocalMaxima.structure.rename_struct as RNM

def register_parser(subparsers):
    parser = subparsers.add_parser('rename', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-s","--structure", metavar="FILE", help="""Path to structure
                              in which entities are renamed.""", required=True)
    requiredArgs.add_argument("--se1", type=str, help="""Structural expression for 
                        selecting model, chain, residue or atom to be renamed.
                        """, required=True)
    requiredArgs.add_argument("--se2", type=str, help="""Structural expression with
                        specify how entities are added from one structure to another.
                        To avoid clashes 
                        """, required=True)
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-o","--out", metavar="FILE",
                               help="""Output file with renamed entities. If no output
                               options and file are given, renaming will be done on 
                               --structure files.""",
                               required=False)
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.structure)
    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        renm_o = RNM.rename(structure.dir_path, structure.file_name)
        se1 = SE.structural_expression(options.se1)
        se1_atoms = se1.atom_list(renm_o.strct)
        # Don't get atoms list for the second structural expression.
        se2 = SE.structural_expression(options.se2)
        renm_o.check_renaming_expressions(se1, se2)
        if options.out:
            renm_o.output_new_s(options.out)
        else:
            renm_o.output_new_s(options.structure)
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return """This comands renames entities (i.e models, chains, residues or atoms) 
           within one strucuture. The two structural expressions must have an equivalent
           number of entities.
           """
def usage():
    return '''\nmodifr rename --structure 1GIA.cif --se1 m[0]c[A]r[3:6]a[CA] --se2
m[0]c[A]r[3:6]a[CA] --out 1GDD_aligned_completed_with_1GIA.pdb'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
