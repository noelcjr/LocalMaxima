import argparse
import sys
import LocalMaxima.structure.structure3D as strct
import LocalMaxima.structure.MolecularRigidManipulations as MRP

def register_parser(subparsers):
    parser = subparsers.add_parser('copy', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--frm", metavar="FILE", help="""Path to structure
                        FRoM which entities are copied to another one.""", required=True)
    parser.add_argument("--se1", type=str, help="""Structural expression for 
                        selecting model, chain or residue location to be added
                        to another structure. Individual atoms can't be added 
                        for now (tricky), so a[] in strucutral expresion can 
                        only have an asterisc inside. m[], c[], r[] can have 
                        anything inside, but some restrictions apply to meet
                        PDB/CIF formating constrains.
                                Ex: m[0]c[*]r[*]a[*] to add a model;
                                    m[*]c[B]r[*]a[*] to add a chain to all models
                                    m[*]c[B]r[6]a[*] to add residue in position 6 of chain b.
                        """, required=True)
    parser.add_argument("--to", metavar="FILE", type=str, help="""Path to structure 
                        to which models, chains or residues will be added.""", required=True)
    parser.add_argument("--se2", type=str, help="""Structural expression to 
                        specify how entities are added from one structure to another.
                        To avoid clashes between model numbers, chain identifiers
                        or residue numbers. --se1 and --se2 must be identical with the
                        exception of entities that will be renamed when added.
                            Ex: --se1 m[0]c[*]r[*]a[*] and --se2 m[1]c[*]r[*]a[*]
                                This will rename model 0 to model 1 to avoid a clash.
                                --se1 m[*]c[B]r[*]a[*] and --se2 m[*]c[C]r[*]a[*]
                                This will rename chain B to chain C to avoid a clash.
                                --se1 m[*]c[B]r[6]a[*] and --se2 m[*]c[B]r[40]a[*]
                                This will rename model 0 to model 1 to avoid a clash.
                        """, required=True)
    parser.add_argument("--output", metavar="FILE", type=str, help="""Optional. If present,
                        A third structure file is generated with the original and copied
                        entities.""",required=False)
    parser.add_argument("--ver", action="store_true", default=False, help="Verbose output for debugging.")
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.frm)
    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        print("modifr copy has nor been implemented yet.")
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return """This comands copies entities (i.e models, chains, or residues) from
           from one strucuture to another. It allows for renaming the added entities
           to avoid name conflicts in the receiving structures.
           """
def usage():
    return '''\nstructure.py align --fix 1GIA.cif --fixatoms m[0];c[A];r[3:6];a[CA,CG] --mov 1GDD.pdb --movatoms 
           m[0];c[A];r[3:6];a[CA,CG] --out 1GDD_aligned_completed_with_1GIA.pdb'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
