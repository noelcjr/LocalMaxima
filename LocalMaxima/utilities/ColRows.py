#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:11:05 2020

@author: noel
"""
import os
import sys

class ColumnsRows(object):
    def __init__(self, str_path):
        self.file_name = os.path.basename(str_path).split('.')[0]
        self.file_sufix = os.path.basename(str_path).split('.')[1]
        self.dir_path = os.path.dirname(str_path)
        if self.file_sufix == 'pdb':
            if self.dir_path == "":
                self.PDB_file = self.file_name+'.'+self.file_sufix
            else:
                self.PDB_file = self.dir_path+'/'+self.file_name+'.'+self.file_sufix
            if not os.path.exists(self.PDB_file):
                print("Error: No valid --input option entered. Type -h or --help.")
                print("       PDB input file path:"+self.PDB_file)
                sys.exit(1)
            self.inFile = open(self.PDB_file, 'r')
        else:
            print("Error: No valid --input option entered. Type -h or --help.")
            print("       pdb input file required: "+self.file_name+"."+
                  self.file_sufix)
            sys.exit(1)
    
    def fix_pdb_from_CHARMM(self,a=22, b=73):
        ''' This fixes the problem of CHARMM not generating a chain id. 
        But needs to be tested in more systems '''
        contents = self.inFile.read()
        contents_list = contents.split('\n')
        for i in range(len(contents_list)):
            if contents_list[i][0:4] == 'ATOM':
                temp = list(contents_list[i])
                # The following two indexes are off by one becasue list index 
                # start at 0, and lines in files at 1.
                temp[a-1] = temp[b-1]
                contents_list[i] = ''.join(temp)
        outFile = open(self.PDB_file, 'w')
        for i in contents_list:
            outFile.write(i + '\n')
        outFile.close()

    def delete_column(self):
        pass
    
    def move_row(self):
        pass

    def delete_row(self):
        pass
