import argparse
import LocalMaxima.structure.structure3D as strct

def register_parser(subparsers):
    parser = subparsers.add_parser('summary', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", help="Structure Input File.", required=True)
    parser.add_argument("--sequence", default='\0', metavar="FILE", help="Optional Sequence FASTA File.", required=False)
    parser.set_defaults(func=run)

def run(options):
    structure = strct.structure(options.structure)

    if (structure.file_sufix.lower() == 'pdb') or (structure.file_sufix.lower() == 'cif'):
        if options.sequence != '\0':
            structure.sequence_path = options.sequence
            structure.add_sequence('fasta')
        structure.summary()

def description():
    return '''This command gives a brief or summarized output of the structure as a guide for other commands.
        The output gives information that is relevant to biopython, and that options that will help decide
        command line inputs for other options for this program such as --align and --extract. The output gives
        information on number of models and chains.(TODO: This option is not coded yet)'''

def usage():
    return '\npdb_cif.py summary --structure 1HIU.pdb --sequence 1HIU.fasta.txt\n or type strctr summary --help'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
