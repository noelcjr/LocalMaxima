import argparse

import LocalMaxima.ThirdPartySoftware.charmm.GenInputs as IO

def register_parser(subparsers):
    parser = subparsers.add_parser('pdb', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", help="Input PDB file.", 
                        required=True)
    parser.add_argument('--terminals', type=str, 
                        help="""Format: A,none,CTER:B,ACE,none (e.i. for chain 
                        A, not ACE and a CTER. For chain B, an Ace, and no 
                        CTER.)""",required=True)
    parser.set_defaults(func=run)

def run(options):
    gen_charmm = IO.GenInputs()
    gen_charmm.prepare_pdb_for_charmm(options.structure,options.terminals)

def description():
    return '''This command takes a PDB file, and it prepares it so hydrogen 
        and heavy atoms are added if missing. It identifies HIS residues as 
        any of the three histidine types in CHARMM, HSD, HSE or HSP, depending 
        on their probable protonation state. It adds missing atoms to the 
        protein as well as terminals to the protein ends. Only ACE and NTERM
        terminals available for now.'''

def usage():
    return '\nprepr pdb --structure 1BRS.pdb\n'\
           'prepr pdb --structure 1BRS.pdb --terminals A,none,CTER:B,ACE,none'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
