#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:40:04 2021

@author: noel
"""
import argparse
import sys
import LocalMaxima.structure.structure3D as strct
import LocalMaxima.structure.MolecularRigidManipulations as MRP

def register_parser(subparsers):
    parser = subparsers.add_parser('align1D', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--fix", metavar="FILE", help="Path to fixed Reference Structure.", required=True)
    parser.add_argument("--fixatoms", type=str, help="Structural expression for fixed structure.", required=True)
    parser.add_argument("--mov", metavar="FILE", type=str, help="Path to moved or fitted structure.", required=True)
    parser.add_argument("--movatoms", type=str, help="Structural expression for moved structure.", required=True)
    parser.add_argument("--output", metavar="FILE", type=str, help="Path for moved structure.",required=False)
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
    return '''This command aligns two structures's sets of atoms specified by 
           structural expressions. The --fix structure is the reference structure. 
           The --mov structure option is translated over the --fix structure. 
           Both --fixatoms and --moveatoms are structural expression that define
           the set of atoms to be aligned for both structures. The number of atoms 
           defined by these two structural expressions must be equal. \n\n
           The --out specifies a filename for a new aligned structure, 
           and it is optional. If no --out is present, the structure defined 
           by --mov will be overwritten.
           '''
def usage():
    return '''\nstructure.py align --fix 1GIA.cif --fixatoms m[0]c[A]r[3:6]a[CA,CG] --mov 1GDD.pdb --movatoms 
           m[0]c[A]r[3:6]a[CA,CG] --out 1GDD_aligned_completed_with_1GIA.pdb\n or type strctr align --help'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
