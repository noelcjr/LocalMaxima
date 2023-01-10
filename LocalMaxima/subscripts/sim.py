#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 21:51:31 2018

@author: noel
"""
import argparse
import optparse
import sys
from LocalMaxima.simulations import monoatomic

def register_parser(subparsers):
    parser = subparsers.add_parser('sim', usage=usage(), description=description())
    add_arguments(parser)
    
def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-m","--mono", type=str,
                               help="OUT files to analyze. Multiple groups of OUTs",
                               required=False, nargs='+')
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-c","--columns", type=str,
                              help="Simulations outputs have no column headers.",
                              required=False, default="somedefault")
    parser.set_defaults(func=run)
def description():
    return '''runs monoatimic simulations.'''
def usage():
    return '''\naguan sim --mono'''
def run(options):
    print("Inside Monoatomic")
if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
