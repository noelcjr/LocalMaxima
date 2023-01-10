#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 21:26:58 2018

@author: noel
"""
import argparse
import sys
import os
import LocalMaxima.utilities.plototypes as pt

def register_parser(subparsers):
    parser = subparsers.add_parser('zmplot', usage=usage(), \
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-i","--input", metavar="FILE",
                              help="Enter a character dellimited file, default comma.", 
                              required=True)
    requiredArgs.add_argument("-o","--out", metavar="FILE",
                              help="Enter file out put name.",
                              required=True)
    # TODO: for more control of X and Y axis ranges. develop this. Otherwise,
    #       default is min and max of index for X, and of the whole DF for Y
    #requiredArgs.add_argument("-x","--xy", type=str,
    #                           help="X and Y axis ranges and intervals.\
    #                           Ex: 0,10,10,ns:5,15,1,nm Separate X and Y with ':'.\
    #                           If only one entered, it applies to the whole DF.\
    #                           Otherwise enter one for each column. Error, oherwise.\
    #                           Values for X will override index in data frame.",
    #                           required=True, nargs='+')
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-d","--delimiter", type=str,
                               help="A delimiter other than coma. Ex: \t",
                               required=False)
    optionalArgs.add_argument("-f","--flag", type=str, 
                               help="A flag indicates type of plot wanted",
                               required=False)
    parser.set_defaults(func=run)

def run(options):
    if not os.path.exists(options.input):
        print("ERROR: input file does not exist.", options.input)
        sys.exit(1)
    pltyps = pt.plototypes(options,"zmplot")
    pltyps.OH_LJ_Average()  # Check min max parameters are selected right.
    #pltyps.OH_LJ_Zhistog() # Time consuming and not right output.

def description():
    return '''Makes plots for special analysis, not a general purpose ploter, but
              some of its functions my be general enough that can be applied to
              other data/analysis.
           '''
def usage():
    return '''\nnplot --input'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
