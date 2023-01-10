#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 21:46:18 2021

@author: noel
"""
import sys
import LocalMaxima.structure.XYZ_Formats as xyzf

def register_parser(subparsers):
    parser = subparsers.add_parser('gen_psf', usage=usage(), description=description())
    add_arguments(parser)
    
def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-o","--out", type=str,
                               help="Output file name. Must have psf sufix.",
                               required=True)   
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-n","--mols", type=str,
                               help="A comma separated sets of ChainID,MoleculeName,MoleculeNumber,\
                               to identify molecules to be constructed in the PSF file.\
                                  ChainID: A single letter chain ID for each molecule entry.\
                                  MoleculeName: As they appear in CHARMM top_27 parameter file.and number to build.\
                                  MoleculeNumber: The number of molecules to be generated.",
                               required=False, nargs='+')
    optionalArgs.add_argument("-l","--list", action="store_true",
                              help="List molecules that can be built.",
                              required=False)
    parser.set_defaults(func=run)
def description():
    return '''Generates PSF files for simulations by listing molecules in triplets:
           ChainID,MoleculeName,MoleculeNumber. Ex, A,tip5p,216 B,CLA,108. 
           Long polymers are not dealt by this program yet. This program will
           not add bonds or angles to of amino or nucleic acid polimers.
           This capability will be added later.
           '''
def usage():
    return '''\naguan gen_psf --mols '''
def run(options):
    if options.mols:
        if options.list:
            print("ERROR: options --mols and --list selected. Only one allowed.")
            sys.exit()
        psf_file = xyzf.psf(options.out)
        psf_file.gen_psf_file(options)
        # TODO collect all needed parameters and send it to XYZ_Formats.psf to out a psf file.
        # consider loading parameters inside PSF object.
    else:
        CP.read_charmm_FF()
if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)