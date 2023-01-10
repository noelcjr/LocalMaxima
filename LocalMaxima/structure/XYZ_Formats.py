#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 20:05:07 2021

@author: noel
"""
import sys
import numpy as np
import Bio.PDB as PDB
import struct
import LocalMaxima.utilities.CHARMM_Parser as CP
import os
import pandas as pd
import datetime as dt

class dcd(object):
    def __init__(self,options):
        self.NSET = 0
        self.ISTART = 0 
        self.NSAVC = 0
        self.NAMNF = 0
        self.DELTA = 0
        self.title_size = 0
        self.title_size2 = 0
        self.N = 0
        self.cord_idx = 0
        self.XYZ = []
        self.INIT = options.init
        self.freeindexes= []
        self.psf = options.psf
        self.dcds = options.dcd
        self.atms = options.atoms
        self.outcsv = options.outcsv
    def read_OH_LJ_Z(self, init):
        """ Reads multpiple DCDs and produces average Z distance for atoms,
        in this ase OH and LJ, and it also outputs the Z coordoordinate of these
        atoms into a list that is stored in a DataFrame, and output as csv file.
        """
        psf_file = psf(self.psf)
        psf_file.read_psf_file()
        LJ = []
        OH = []
        count = 0
        a = self.atms.split(',')
        for i in psf_file.atmtyp1:
            if i ==  a[0]:
                OH.append(count)
            elif i == a[1]:
                LJ.append(count)
            count += 4
        columns_dict = {}
        count2 = 0
        for i in self.dcds:
            count2 += 1
            dcd_path = os.path.dirname(i)
            dcd_box = dcd_path.split('/')[-1]
            individual_files = os.path.basename(i).split(",")
            dcd_files = [dcd_path+"/"+jj for jj in individual_files]
            dcdID = self.check_file_names(individual_files)
            #print(count2,i, dcdID)
            OH_mZ = []
            LJ_mZ = []
            OHs_Zl = []
            LJs_Zl = []
            for fdcd in dcd_files:
                with open(fdcd, mode='rb') as file: # b is important -> binary
                    self.fileContent = file.read()
                self._read_header()
                count = self.N*4
                #print("   ",fdcd,self.INIT,self.NSET)
                for ii in range(self.INIT,self.NSET):
                    OH_Zloc = []
                    LJ_Zloc = []
                    for k in OH:
                        a = self.cord_idx+k
                        a += count + 8
                        a += count + 8
                        Z = struct.unpack(">f", self.fileContent[a:a+4])[0]
                        OH_Zloc.append(Z)
                    for j in LJ:
                        a = self.cord_idx+j
                        a += count + 8
                        a += count + 8
                        Z = struct.unpack(">f", self.fileContent[a:a+4])[0]
                        LJ_Zloc.append(Z)
                    OH_mZ.append(sum(OH_Zloc)/len(OH_Zloc))
                    LJ_mZ.append(sum(LJ_Zloc)/len(LJ_Zloc))
                    OHs_Zl.append(OH_Zloc)
                    LJs_Zl.append(LJ_Zloc)
                    self.cord_idx = (a+12)
                self.INIT = self.NSET
            self.INIT = init
            columns_dict[dcd_box+"_OH_"+dcdID] = pd.Series(OH_mZ)
            columns_dict[dcd_box+"_LJ_"+dcdID] = pd.Series(LJ_mZ)
            columns_dict[dcd_box+"_OH_"+dcdID+"_L"] = pd.Series(OHs_Zl)
            columns_dict[dcd_box+"_LJ_"+dcdID+"_L"] = pd.Series(LJs_Zl)
        DF = pd.DataFrame(columns_dict)
        print("Getting Z distribution OH and monoatomic.")
        print("Indexed scaled by multiplying by 0.002 and get nanosecons")
        DF.index = (DF.index+1)*0.02
        DF.index.name = 'ns'
        print(DF.shape)
        print(DF.columns," output ",self.outcsv,"\n")
        DF.to_csv(self.outcsv)
    
    def check_file_names(self,individual_files):
        # This same class also exist in utilities/plototypes.
        first_file = True
        for indF in individual_files:
            current_prefix_label = indF.split('.')[0].split('_')[0]
            if first_file:
                prefix_label = current_prefix_label
                first_file = False
            else:
                if current_prefix_label != prefix_label:
                    print("""WARNING:        The first underscored separated 
                    string from DCDs must be identical. This will help identify
                    multiple DCDs as part of the same simulation. For example:
                    A naming conventions such as A_1.dcd, A_2.dcd ... will 
                    represent two DCDs belonging to the same contiguous 
                    simulation. If this lables are different, the first one
                    will be used to identified the concatenated trajectories.""")
        return prefix_label
        
    def _read_header(self):
        if struct.unpack(">i", self.fileContent[:4])[0] != 84:
            print("ERROR")
        d = ''
        for i in [t.decode("utf-8") for t in struct.unpack("cccc", self.fileContent[4:8])]:
            d += i
        if d != 'CORD':
            print("ERROR: Only DCDs with CORD can be read for now.")
        self.NSET, self.ISTART, self.NSAVC = struct.unpack(">iii", self.fileContent[8:20])
        self.NAMNF = struct.unpack(">i", self.fileContent[40:44])[0]
        self.DELTA = struct.unpack(">f", self.fileContent[44:48])[0]
        if struct.unpack(">i", self.fileContent[88:92])[0] != 84:
            print("ERROR: 88.")
            sys.exit(1)
        self.title_size = struct.unpack(">i", self.fileContent[92:96])[0] # = 164
        self.cord_idx = 96+self.title_size  # 260
        self.title_size2 = struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]
        if self.title_size != self.title_size2:
            print("ERROR: title wrappers are not equal")
            sys.exit(1)
        self.cord_idx += 4
        if struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0] != 4:
            print("ERROR: 88.")
            sys.exit(1)
        self.cord_idx += 4
        self.N = struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]
        self.cord_idx += 4  
        if struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0] != 4:
            print("ERROR: 88.")
            sys.exit(1)
        if self.NAMNF !=0:
            self.cord_idx += 4
            self.freeindexes = [0 for i in range(self.N-self.NAMNF)]
            freeindx = struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]
            self.cord_idx += 4
            if freeindx != ((self.N-self.NAMNF)*4):
                print("ERROR: bad DCD format.")
                sys.exit(1)
            for i in range(self.N-self.NAMNF):
                self.freeindexes[i] = struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]
                self.cord_idx += 4
            freeindx2 = struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]
            self.cord_idx += 4
            if freeindx2 != ((self.N-self.NAMNF)*4):
                print("ERROR: bad DCD format.")
                sys.exit(1)
        self.cord_idx +=  4
        if (self.N*4) != struct.unpack(">i", self.fileContent[self.cord_idx:self.cord_idx+4])[0]:
            print("ERROR: DCD reading out of sync.")
            sys.exit(1)
        self.cord_idx += 4
    def _read_body(self):
        count = self.N*4
        for i in range(self.INIT,self.NSET):
            for j in self.N:
                a = self.cord_idx+j
                X = struct.unpack(">f", self.fileContent[a:a+4])[0]
                a += count + 8
                Y = struct.unpack(">f", self.fileContent[a:a+4])[0]
                a += count + 8
                Z = struct.unpack(">f", self.fileContent[a:a+4])[0]
                self.XYZ.append(np.array(X,Y,Z))
            self.cord_idx = (a+12)
        self.INIT = self.NSET

class crd(object):
    '''
     This fixes when there are multiple chains and the residue number gets
     reset at the beginning of each chain. With this fix, residue numbers
     will be renumbered
      1 - 5 Integer Atom no. sequential
     6 - 10 Integer ires Residue position from file beginning
    11 - 11 Space
    12 - 15 Achar resname Residue name
    16 - 16 Space
    17 - 20 Achar type Atom type, IUPAC name of atom left justified
    21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
    31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
    41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
    51 - 51 Space
    52 - 55 Achar segid Segment identifier
    56 - 56 Space
    57 - 60 Achar resid Residue number within the segment
    61 - 70 Real(10.5) Weighting array value
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    '''
    # https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/io/#IOFORM
    def __init__(self, crd_file):
        if crd_file:
            self.CRD_file = crd_file
            self.crd_lines = self.read_crd_file()
            self.col_names = ['atm_num','res_num','res_name','sp1','atm_nam','x','y','z','sp2','chn_id','sp3','res_num2','holder']
            self.form_line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
            self.crd_lines = self.read_crd_file()
        else:
            pass

    def read_crd_file(self):
        """Read each line in a CRD file. A line example: '     1    1 ALA  CAY  -24.37800  14.19800  -2.71000 A    1        0.00000'
           will be processed according to CRD column groups designated for fields in 'col_length'.
           This function also reads the header/comments of the file, and the number of atoms right below.
        """
        self.col_length = [(0,5),(5,11),(11,15),(15,16),(16,20),(20,30),(30,40),(40,50),(50,51),(51,55),(55,56),(56,60),(60,72)]
        inFile = open(self.CRD_file, 'r')
        contents = inFile.read()
        self.comments = []
        self.aa = []  
        self.aaid = []  
        self.atmtyp1 = [] 
        self.chain_id = []  
        self.xyz = []  
        self.occupancy = []  
        B_iso = [] 
        #  the next loop could probably be simpler, but I made sure it can handle atoms with 
        atom_count = 0
        number_atoms = 0
        for i in contents.split('\n')[:-1]:
            if i[0] == '*':
                self.comments.append(i)
            elif len(i) == 5:
                number_atoms = int(i.strip())
            else:
                self.aa.append(int(i[self.col_length[1][0]:self.col_length[1][1]].strip()))
                self.aaid.append(i[self.col_length[2][0]:self.col_length[2][1]].strip())
                self.atmtyp1.append(i[self.col_length[4][0]:self.col_length[4][1]].strip())
                self.chain_id.append(i[self.col_length[9][0]:self.col_length[9][1]].strip())
                self.xyz.append(np.array((float(i[self.col_length[5][0]:self.col_length[5][1]].strip()),
                                          float(i[self.col_length[6][0]:self.col_length[6][1]].strip()),
                                          float(i[self.col_length[7][0]:self.col_length[7][1]].strip()))))
                self.occupancy.append(float(i[self.col_length[12][0]:self.col_length[12][1]].strip()))
                atom_count += 1
    # TODO check /home/noel/Datasets_nobackup/Constructs/Insulin_n_Barnase_barstar_2/proximity_optimization.py
    # for a corrected line to out put PDB directly from coordinates.
    # also check that lines do not come pre-formatted and that is why they are
    # different from proximity_optimization
    def write_pdb_from_crd(self,basedir,filename,lines):
        ''' TODO: write pdbdirectly with all chains the right way.
        http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  6        Record name   "ATOM  "
         7 - 11        Integer       serial       Atom  serial number.
        13 - 16        Atom          name         Atom name.
        17             Character     altLoc       Alternate location indicator.
        18 - 20        Residue name  resName      Residue name.
        22             Character     chainID      Chain identifier.
        23 - 26        Integer       resSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     occupancy    Occupancy.
        61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        77 - 78        LString(2)    element      Element symbol, right-justified.
        79 - 80        LString(2)    charge       Charge  on the atom.
        '''
        line = '{:6}{:>4} {:5}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
        pdb_converted_lines = []
        for i in lines:
            pdb_converted_lines.append(line.format('ATOM',i[0],i[4],'',i[2],i[9],int(i[1]),'',float(i[5]),float(i[6]),\
                                                   float(i[7]),float('0'),float('0'),i[9],''))
        if basedir == '':
            basedir = '.'
        if basedir[-1] != "/":
            basedir += "/"
        outFile = open(basedir+filename, 'w')
        for i in pdb_converted_lines:
            outFile.write(i + '\n')
        outFile.close()
        
    def write_crd_from_pdb(self,basedir,first,last=''):
        '''
         This fixes when there are multiple chains and the residue number gets
         reset at the beginning of each chain. With this fix, residue numbers  
         will be renumbered
          1 - 5 Integer Atom no. sequential
         6 - 10 Integer ires Residue position from file beginning
        11 - 11 Spacevi 
        12 - 15 Achar resname Residue name
        16 - 16 Space
        17 - 20 Achar type Atom type, IUPAC name of atom left justified
        21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
        31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
        41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
        51 - 51 Space
        52 - 55 Achar segid Segment identifier
        56 - 56 Space
        57 - 60 Achar resid Residue number within the segment
        61 - 70 Real(10.5) Weighting array value
        '''
        line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
        if last == '':
            file_name = first.split('.')
            if len(file_name) != 2:
                print("Error: A full PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                print("       Please enter a correct file name. Script will exit now.")
                sys.exit(1)
            elif len(file_name) == 2:
                if file_name[1] != 'pdb':
                    print("Error: A correct PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                    print("       Please enter a correct file name with lower case suffix. Script will exit now.")
                    sys.exit(1)
                else:
                    structure = PDB.PDBParser(QUIET=True).get_structure('pdb', basedir+str(first))
                    count_atoms = 1
                    count_resid = 1
                    lines = []
                    for j in structure.get_atoms():
                        if count_atoms == 1:
                            current_resid = j.get_parent().id[1]
                        if current_resid != j.get_parent().id[1]:
                            current_resid = j.get_parent().id[1]
                            count_resid += 1
                        lines.append(line.format(count_atoms, count_resid, j.get_parent().resname, j.get_full_id()[4][0], j.get_coord()[0], 
                              j.get_coord()[1], j.get_coord()[2], j.get_full_id()[2], j.get_parent().id[1], 0.00))
                        count_atoms += 1
                    count_atoms -= 1    
                    len(str(count_atoms))
                    lines.insert(0,(5-len(str(count_atoms)))*' '+str(count_atoms))
                    lines.insert(0,'*')
                    lines.insert(0,'* Date: XX/XX/XXXX by Noel Carrascal')
                    lines.insert(0,'* Converting PDB to CRD')
                    
                    outFile = open(basedir+str(file_name[0])+'.crd', 'w')
                    for j in lines:
                        outFile.write(j+'\n')
                    outFile.close()
            else:
                print("Error: A correct PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                print("       Please enter a correct file name. Script will exit now.")
                sys.exit(1)               
        else:
            if (type(int(first)) is int) and (type(int(last)) is int):
                for i in range(0,990):
                    structure = PDB.PDBParser(QUIET=True).get_structure('pdb', basedir+str(i)+'.pdb')
                    count_atoms = 1
                    count_resid = 1
                    lines = []
                    for j in structure.get_atoms():
                        if count_atoms == 1:
                            current_resid = j.get_parent().id[1]
                        if current_resid != j.get_parent().id[1]:
                            current_resid = j.get_parent().id[1]
                            count_resid += 1
                        lines.append(line.format(count_atoms, count_resid, j.get_parent().resname, j.get_full_id()[4][0], j.get_coord()[0], 
                              j.get_coord()[1], j.get_coord()[2], j.get_full_id()[2], j.get_parent().id[1], 0.00))
                        count_atoms += 1
                    count_atoms -= 1    
                    len(str(count_atoms))
                    lines.insert(0,(5-len(str(count_atoms)))*' '+str(count_atoms))
                    lines.insert(0,'*')
                    lines.insert(0,'* Date: XX/XX/XXXX by Noel Carrascal')
                    lines.insert(0,'* Converting PDB to CRD')
                    
                    outFile = open(basedir+str(i)+'.crd', 'w')
                    for j in lines:
                        outFile.write(j+'\n')
                    outFile.close()

    # TODO this method coulfd be giving a wrong output.
    def write_crd(self,structure, output_name, model=0):
        line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
        atom_count = 0
        resi_count = 0
        lines = []
        for a in structure.get_models():
            if a.id == model:
                for b in a.get_chains():
                    for c in b.get_residues():
                        resi_count += 1
                        for d in c.get_atoms():
                            atom_count += 1
                            # The following uncomment lines seem to be needed only
                            # if the first residue is not 1. do we need them?
                            #if count_atoms == 1:
                            #    current_resid = d.get_parent().id[1]
                            # TODO could this cause a bug?
                            #if current_resid != d.get_parent().id[1]:
                            #    current_resid = d.get_parent().id[1]
                            #    count_resid += 1
                            lines.append(line.format(atom_count, resi_count, \
                                                     c.get_resname(), \
                                                     d.get_name(), \
                                                     d.get_coord()[0], \
                                                     d.get_coord()[1], \
                                                     d.get_coord()[2], b.id, \
                                                     c.get_full_id()[3][1], 0))
        outFile = open(output_name+".crd", 'w')
        outFile.write('* CRD generated with fixed HS \n')
        outFile.write('* Name: NN\n')
        outFile.write('* Title: New Molecule Generated by pdb_preparation.py\n')
        outFile.write('*\n')
        outFile.write('{:>5}'.format(len(lines)) + '\n')
        for i in lines:
            outFile.write(i + '\n')
        outFile.close()

class psf(object):
    # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html
    def __init__(self, psf_file):
        self.PSF_file = psf_file
        # Must be called separately to accomodate other psf functions like generation
        #self.psf_lines = self.read_psf_file()

    def read_psf_file(self):
        """ This function reads the !NATOM section of a PSF file in xplor format.
            A line example: '         1 A        1        ALA      CAY      CT3   -0.270000       12.0110           0   0.00000     -0.301140E-02' 
            It gets from the atmtyp2, which is the second atom type after the residue name (e.g ALA   CAY  CT3 <-this last one->),
            and right after that the charge and molecular weight.
        """
        self.title = []
        self.atomnumber = []
        self.segment = []  
        self.resnumber = []
        self.resname = []
        self.atmtyp1 = []
        self.atmtyp2 = []
        self.charge = []
        self.weight = []
        with open(self.PSF_file, 'r') as f:
            contents = f.readlines()
        atom_count = 0 
        number_atoms = 0
        current_field = ''
        self.firstLine = contents[0]
        #print("DDD",self.firstLine.strip().split()[1])
        valid_fields = ['!NTITLE','!NATOM','!NBOND: bonds']
        for i in contents:
           if i[i.find("!"):].strip() in valid_fields:
               div1 = i.find("!")
               current_field = i[div1:].strip()
               if current_field == "!NATOM":
                   number_atoms = int(i[0:div1])
               elif current_field == "!NTITLE":
                   pass
           else:
               if current_field == '!NTITLE':
                   self.title.append(i)     
               elif current_field == '!NATOM' and len(i) != 1:
                   if atom_count < number_atoms:
                       if self.firstLine.strip().split()[1] == 'EXT':
                           self.atomnumber.append(int(i[:11].strip()))
                           self.segment.append(i[11:20].strip())
                           self.resnumber.append(int(i[20:29].strip()))
                           self.resname.append(i[29:38].strip())
                           self.atmtyp1.append(i[38:47].strip())                       
                           self.atmtyp2.append(i[47:55].strip())
                           self.charge.append(float(i[55:71].strip()))
                           self.weight.append(float(i[71:82].strip()))
                       else:
                           self.atomnumber.append(int(i[:9].strip()))
                           self.segment.append(i[9:14].strip())
                           self.resnumber.append(int(i[14:19].strip()))
                           self.resname.append(i[19:24].strip())
                           self.atmtyp1.append(i[24:29].strip())
                           self.atmtyp2.append(i[29:33].strip())
                           self.charge.append(float(i[33:49].strip()))
                           self.weight.append(float(i[49:59].strip()))
                       atom_count += 1
               elif current_field == '!NBOND: bonds':
                   pass
        if atom_count != number_atoms:
            print("ERROR: Number of atoms in PSF is different from the nunmber read.")
            print("       PSF file might be corrupted. Exiting now.")
            print("       number_atoms =",number_atoms," atom_count =",atom_count)
            sys.exit()

    def gen_psf_file(self,options):
        """ Generates a PSF file"""
        #if not os.path.exists(file_path):
        #    print(f"ERROR: path {file_path} does not exist")
        allowed_to_build = ["TIP3","CLA"]
        entries = []
        num_entries = 0
        for i in range(len(options.mols)):
            entries.append(options.mols[i].split(','))            
        PAR = CP.read_charmm_FF()
        # Count atoms in Entries
        for i in entries:
            atms = 0
            for j in PAR.AA[i[1]].atoms:
                atms += len(j)
            i.append(atms)
        file_name = os.path.basename(options.out).split('.')[0]
        file_suffix = os.path.basename(options.out).split('.')[1]
        dir_path = os.path.dirname(options.out)
        if file_suffix == 'psf':
            outFile = open(options.out, 'w')
            outFile.write("PSF CMAP\n")
            outFile.write("\n")
            outFile.write("{:>8d} !NTITLE".format(len(options.mols)+1)+"\n")
            for i in range(len(entries)):
                if isinstance(int(entries[i][-1]),int):
                    num_entries += int(entries[i][-2])*entries[i][-1]
                else:
                    print(f"ERROR: instance {entries[i][-1]} in i must be an integer.")
                    print("      No PSF file generated.")
                    sys.exit()
                outFile.write(f"REMARKS PSF file with {entries[i][2]} Molecules of {entries[i][1]}\n")
            outFile.write(f"REMARKS Date: {str(dt.date.today())} Created with 'aguan gen_psf --help'.\n\n")
            outFile.write("{:>8d} !NATOM".format(num_entries))
            outFile.write("\n")
            num_atoms = 1
            num_molecs = 1
            line = "{:>8d}{:>2}    {:<5d}{:>4} {:<3}  {:<3} {:>12.6f}{:>14.5f}{:>41}"
            bond_indx = []
            for i in entries:
                if i[1] in PAR.AA:
                    chain = i[0]
                    molec = i[1]
                    for j in range(int(i[-2])):
                        for k in PAR.AA[i[1]].atoms:
                            molec = {}
                            for l in k:
                                outFile.write(line.format(num_atoms,
                                                          chain,
                                                          num_molecs,
                                                          i[1],
                                                          l,
                                                          PAR.AA[i[1]].atom_type[l],
                                                          PAR.AA[i[1]].atom_chrg[l],
                                                          float(PAR.am.MASS[PAR.AA[i[1]].atom_type[l]][1]),
                                                          "0   0.00000     -0.301140E-026")+"\n")
                                molec[l] = num_atoms
                                num_atoms += 1
                            for l in PAR.AA[i[1]].bonds:
                                if len(l) > 0:
                                    for m in l:
                                        if (m[0] in molec) and (m[1] in molec):
                                            bond_indx.append(molec[m[0]])
                                            bond_indx.append(molec[m[1]]) 
                        num_molecs += 1
                else:
                    print("ERROR:",i[1],"NOT Found in parameter files.")
                    print("      No PSF file generated.")
                    sys.exit(1)
            outFile.write("\n")
            outFile.write("{:>8d} !NBOND\n".format(int(len(bond_indx)/2)))
            for i in range(0,len(bond_indx),8):
                for j in range(i,i+8,2):
                    outFile.write("{:>8d}{:>8d}".format(bond_indx[j],bond_indx[j+1]))
                outFile.write("\n")
            outFile.write("\n")
            outFile.write("{:>8d} !NTHETA\n".format(0))
            outFile.write("\n")
            outFile.write("{:>8d} !NPHI\n".format(0))
            outFile.write("\n")
            outFile.write("{:>8d} !NIMPHI\n".format(0))
            outFile.write("\n")
            outFile.write("{:>8d} !NDON\n".format(0))
            outFile.write("\n")
            outFile.write("{:>8d} !NACC\n".format(0))
            if os.path.exists(options.out):
                pass
            outFile.close()
        else:
            print("ERROR: Output file must have psf suffix.")
            sys.exit()