import argparse
import sys
import LocalMaxima.structure.structure3D as strct

def register_parser(subparsers):
    parser = subparsers.add_parser('gaps', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", \
                        help="PDB or CIF Structure File.", required=True)
    parser.add_argument("--sequence", default='\0', metavar="FILE", \
                        help="Sequence FASTA File. This file is considered \
                            only when a PDB structure file is used.", \
                                required=False)
    parser.add_argument('--num', type=int, action="store", default=60,\
                        help="This option takes an integer value for the width\
                            of the alignment displays. Default is 60.")
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.structure)
    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        if options.sequence != '\0':
            print("Number of columns:"+str(options.num))
            structure.sequence_path = options.sequence
            structure.add_sequence('fasta')
        structure.gap_report(options.num)
    else:
        print("ERROR: Unrecognized structure file format.")
        print("       Program will exit without results.")
        sys.exit(1)

def description():
    return """This command detects gaps in the crystal structure of a protein. 
    The search requires at least a structure PDB or CIF file. When only a CIF 
    or PDB file is provided, gaps are checked by tracking monotonically increa-
    sing residue numbers from 1 to the last residue in the structure; this 
    search has the potential to not report missing residues in the C-terminal 
    of the structure as gaps because the last amino acid in the structure is 
    assumed to be the last one in the sequence.
                                                             
    All gaps in crystalographic structures can be detected when a FASTA file
    is provided in the arguments. An alignment is done between a sequence of 
    amino acids present in the structure with the ones in the FASTA file. This
    will detect missing amino acids in the C-terminal
    """

def usage():
    return '\nstrctr gaps --structure 1BRS.pdb --sequence 1BRS.fasta.txt\n or type strctr gaps --help'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
