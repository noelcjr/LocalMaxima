#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 21:26:58 2018

@author: noel
"""
import argparse
import LocalMaxima.energy.CA as ca
import LocalMaxima.ThirdPartySoftware.charmm.GenInputs as IO
import os

def register_parser(subparsers):
    parser = subparsers.add_parser('folding', usage=usage(), \
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):  
    parser.add_argument("--crd", metavar="FILE",help="Enter a CRD file generate\
                        d by CHARMM. Filename will be used as output suffix."
                        ,required=True)
    parser.add_argument("--psf", metavar="FILE",help="Enter PSF file generate\
                        d by CHARMM NOT in XPLOR format.", required=True)
    parser.add_argument("--CA", action="store_true", default=False, help="Ou\
                        tput component analysis atom-number indexes.")
    parser.add_argument('--terminals', type=str,
                         help="""Format: A,none,CTER (e.i. for chain 
                         A, not ACE and a CTER.) This is needed even if it's 
                         already present in the
                         structure to create an equivalent unfolded protein.
                         """,required=True)
    parser.set_defaults(func=run)

def run(options):
    gen_charmm = IO.GenInputs()
    if options.CA:
        gen_charmm.get_gbsa_folding(options)
        analysis = ca.mmgbsa_ca_analysis(os.getcwd())
        analysis.mmgbsa_CA_bindingMatrix(options)
    else:
        # gen_charmm.get_binding_mmgbsa(options)
        gen_charmm.get_folding_mmgbsa(options)

def description():
    return '''Generates de delta_MMGBSA between the folded and unfolded states.
              Generates a file with delta (folded - unfolded) MMGBSA CA of binding
              It works on structures with only one chain because bound chains
              affect born radii of atoms near protein interphases, and results would
              be confused with binding***. Use 'strctr extract ...' for each chain
              of a protein complex before calculating bonding energy.

              *** The shape of a protein chain not in complex is probably different from
                  the same protein in complex. Extracting the chain from the complex will
                  not provide a folding energy strictly speaking, but it is better
                  than getting folding of a chain while still bound to another chain.
           '''
def usage():
    return '''\nmmgbsa bind_ca --xyz, --psf'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)

