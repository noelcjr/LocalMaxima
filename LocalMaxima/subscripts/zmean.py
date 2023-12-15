#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 21:26:58 2018

@author: noel
"""
import argparse
import sys
import os
import LocalMaxima.structure.XYZ_Formats as xyzf

def register_parser(subparsers):
    parser = subparsers.add_parser('zmean', usage=usage(), \
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-p","--psf", metavar="FILE",
                              help="A PSF parameter file.", 
                              required=True)
    requiredArgs.add_argument("-d","--dcd", type=str,
                               help="DCD files to analyze. Multiple groups of DCDs\
                               are allowed as separate arguments. Individual groups\
                               within a directory can be liste separeted by comas and\
                               following the path to the directory. Ex: \
                               2.0X_2.0Y_2.0Z/A_1.dcd,A_2.dcd 2.0X_2.0Y_2.2Z/B_1.dcd\
                               will calculate zdist on merged A_1.dcd and A_2.dcd, \
                               and will be placed in separate columns in a dataframe\
                               from B_1.dcd.",
                               required=True, nargs='+')
    requiredArgs.add_argument("-o","--outcsv", metavar="FILE",
                              help="Enter output file name including the suffix csv.",
                              required=True)
    requiredArgs.add_argument("-a","--atoms", type=str,
                               help="A coma separated string of atom names for \
                               zmean calculations.",
                               required=True)
    requiredArgs.add_argument("-i","--init", type=int,
                               help="The INIT variable in dcds is set wrong. \
                               I won't fix it because VMD doesn't need this is\
                               a patch for now to get the plots.\
                               zmean calculations.",
                               required=True)
    parser.set_defaults(func=run)

def run(options):
    if not os.path.exists(options.psf):
        print("ERROR: PSF input file does not exist.")
        sys.exit(1)
    dcd_zmeans = xyzf.dcd(options)
    dcd_zmeans.read_OH_LJ_Z(options.init)

def description():
    return '''Gets the mean z-axis location and a list of z-axis distributions for
              selected atoms of one or more DCD files. The output of this function
              can be used to build moving z-distributions histograms with, for example:
              ntral zhist --csv A_0.csv --outmp4 A_0.mp4
           '''
def usage():
    return '''\nntraj zmean --psf A_0.psf --dcd A_0.dcd --outcsv A_0.csv --atoms OH,LJ0 --init 0'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
