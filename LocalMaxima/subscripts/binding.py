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
    parser = subparsers.add_parser('binding', usage=usage(), \
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):  
    parser.add_argument("--xyz", metavar="FILE",help="Enter a CRD file generate\
                        d by CHARMM. Filename will be used as output suffix."
                        ,required=True)
    parser.add_argument("--psf", metavar="FILE",help="Enter PSF file generate\
                        d by CHARMM NOT in XPLOR format.", required=True)
    parser.add_argument("--mov", type=str,help="""Enter chain IDs separated by\
                        comas to be moved in the z-axis. For example: 'A' or \
                        "A,B'.""", required=True)
    parser.add_argument("--trans", type=str,action="store", default=5000,help=""\
                        "Distance chain is move for binding energy calcualtion\
                        s in Angstroms. Default 5000 A in Z-axis direction.""",\
                        required=False)
    parser.add_argument("--CA", action="store_true", default=False, help="Ou\
                        tput component analysis atom-number indexes.")
    parser.set_defaults(func=run)

def run(options):
    gen_charmm = IO.GenInputs()
    if options.CA:
        gen_charmm.get_gbsa_binding(options)
        analysis = ca.mmgbsa_ca_analysis(os.getcwd())
        analysis.mmgbsa_CA_bindingMatrix(options)
    else:
        gen_charmm.get_binding_mmgbsa(options)

def description():
    return '''Generates de delta_MMGBSA between the bound and unbound state.
              Generates a file with delta of binding contributions from
              electrostatics, Van der Waals, SASA and self and Pair GB MMGBSA 
              when the --CA is used.
           '''
def usage():
    return '''\nmmgbsa bind_ca --xyz, --psf, --mov, --trans, --CA'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)

