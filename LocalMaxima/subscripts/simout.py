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
    parser = subparsers.add_parser('simout', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-i","--input", type=str,
                               help="OUT files to analyze. Multiple groups of OUTs\
                               are allowed as separate arguments. Individual groups\
                               within a directory can be liste separeted by comas and\
                               following the path to the directory. Ex: \
                               2.0X_2.0Y_2.0Z/A_1.out,A_2.out 2.0X_2.0Y_2.2Z/B_1.out\
                               will calculate zdist on merged A_1.out and A_2.out, \
                               and will be placed in separate columns in a dataframe\
                               from B_1.out.",
                               required=True, nargs='+')
    requiredArgs.add_argument("-o","--out", metavar="FILE",
                              help="Enter file out put name.",
                              required=True)
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-c","--columns", type=str,
                              help="Simulations outputs have no column headers.\
                              The default is a comma separated string:\
                              dt,E,PE,PEvdw,PEee,PErf1,PErf2,PEef,KE,KErot,KEtr,T,Tt,Tr,P",
                              required=False, default="t(ns),E,PE,PEvdw,PEee,PErf1,PErf2,PEef,KE,KErot,KEtr,T,Tt,Tr,P")
    parser.set_defaults(func=run)

def run(options):
    pltyps = pt.plototypes(options,"simout")
    pltyps.plot_sim_outputs(options)

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
