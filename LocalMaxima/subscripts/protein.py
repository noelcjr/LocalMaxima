import argparse

import LocalMaxima.ThirdPartySoftware.charmm.GenInputs as IO
import sys

def register_parser(subparsers):
    parser = subparsers.add_parser('protein', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--name", type=str, help="""Required for both --file and
                        --seq options. It gives a distinguishing name to the new
                        structure. Useful for tracking structures for multiple 
                        structures in a folder. Enter just a prefix. Ex: 2hiu_,
                        12345_ or AATTFF_. Names will be converted to lowercase.
                        """, required=True)
    parser.add_argument("--seq", type=str, help="""Input a sequence in text 
                        format. Ex: 'ALA,GLY,VAL,ILE' or 'AGVI'. The short 
                        notation for protein sequences can only be greater
                        than three amino acids long to avoid confusion with a
                        single three letter aa identifier.""", 
                        required=False)
    parser.add_argument("--chn_id", type=str, help="""Required for --seq option only. 
                        Select a chain_id when generating and unfolded structure
                        from sequence.
                        """, required=False)
    parser.add_argument("--file", type=str, help="""Input a structure PDB or . 
                        Sequence fasta file. Fasta header format: >name:chain
                        Files must have '.pdb' or '.fasta' suffixs.
                        """, 
                        required=False)
    parser.add_argument('--terminals', type=str,
                         help="""Format: A,ACE,CTER (e.i. for chain 
                         A, not ACE and a CTER.) If terminals are added, it must be
                         either terminals or none for both N an C proteins.
                         Only one terminal in either end is not allowed due
                         to CHARMMS options avilable. This argument is optional.
                         """,required=False)

    parser.set_defaults(func=run)

def run(options):
    if options.seq and options.file:
        print("ERROR: You can only generate a structure from a file or a")
        print("       sequence not both. Exit now.")
        sys.error(1)
    gen_charmm = IO.GenInputs()
    if options.seq:
        sequence = gen_charmm.gen_unfolded_xtruct_seq_wrapper(options)
    elif options.file:
        sequence = gen_charmm.gen_unfolded_xtruct_file_wrapper(options)

def description():
    return """This command takes a PDB, fasta file or string to create and
    and unfolded protein structure. The sequence string can be typed directly
    on the terminal following --seq or entered in a file after --file.
    Files will multiple chains will generate separate structure files.
    Non amino acid molecules will be ignore (e.g. water, ions, nucleotides.)
    and all HIS will be protonated as HSE.
    The --seq and --file options are both optional but mutally exclusive, only
    one is allowed.
    PREREQUISITE: Prepare sturcture with 'prepr pdb' and use output pdb as input
    for 'prepr protein' to add possible missing atoms or script will crash.
    """
    #lns = mult1+'\n'+mult2+mult3+mult4+mult5+mult6+mult7+mult8+mult9
    #lns += mult10+mult11+'\n'
    #return print(lns)

def usage():
    return '\nprepr protein --seq 1BRS.pdb\n'\
           'prepr protein --seq ALA,GLY,VAL,ILE\n'\
           'prepr protein --seq AGVI\n'\
           'prepr protein --seq A.fasta\n'\
           'prepr protein --seq 1BRS.pdb --terminals A,none,CTER:B,ACE,none'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
