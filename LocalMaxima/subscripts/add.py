import argparse
import sys
import LocalMaxima.structure.structure3D as strct
import LocalMaxima.structure.MolecularRigidManipulations as MRP

def register_parser(subparsers):
    parser = subparsers.add_parser('add', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--from", metavar="FILE", help="""Path to structure
                        to be added to another one.""", required=True)
    parser.add_argument("--se1", type=str, help="""Structural expression for 
                        selecting model, chain or residue location to be added. 
                        Individual atoms can't be added for now (tricky), so a[] in strucutral
                        expresion can only have an asterisc inside. m[], c[], r[] can
                        have anything inside, but some restrictions apply to meet
                        PDB/CIF formating constrains.volunteer_channel
                                Ex: m[0]c[*]r[*]a[*] to add a model;
                                    m[*]c[B]r[*]a[*] to add a chain to all models
                                    m[*]c[B]r[6]a[*] to add residue in position 6 of chain b.
                        """, required=True)
    parser.add_argument("--to", metavar="FILE", type=str, help="""Path to structure 
                        fitted structure.""", required=True)
    parser.add_argument("--se2", type=str, help="Structural expression for moved structure.", required=True)
    parser.add_argument("--output", metavar="FILE", type=str, help="Path for moved structure.",required=False)
    parser.add_argument("--ver", action="store_true", default=False, help="Verbose output for debugging.")
    parser.set_defaults(func=run)

def run(options):
    fix = strct.structure(options.fix)
    if (fix.file_sufix.lower() != 'pdb') and (fix.file_sufix.lower() != 'cif'):
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)
    mov = strct.structure(options.mov)
    if (mov.file_sufix.lower() != 'pdb') and (mov.file_sufix.lower() != 'cif'):
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)
    alignment = MRP.structure_manipulations()
    alignment.add_structure(fix)
    alignment.add_structure(mov)
    # fix = 0, and mov = 1 in a list inside aligment class. output = 1.
    alignment.align_structures(0,options.fixatoms,1,options.movatoms)
    if options.output is None:
        alignment.output_structure(1)
    else:
        alignment.output_structure(1,options.output)

def description():
    return '''This command adds atoms defined by a regular expression from one structure to another.
              It is recommended that structurs are aligned before adding to another structure.
           '''
def usage():
    return '''\nstrctr add --fix 1GIA.cif --fixatoms m[0];c[A];r[3:6];a[CA,CG] --mov 1GDD.pdb --movatoms 
           m[0];c[A];r[3:6];a[CA,CG] --out 1GDD_aligned_completed_with_1GIA.pdb'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
