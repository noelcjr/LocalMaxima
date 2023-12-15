#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 21:20:28 2023

@author: noel
"""
import LocalMaxima.structure.XYZ_Formats as xyzf

pth = "/home/noel/Projects/Phase_Diag_Aguan_JAVA/Second_Poster_EF_effects/tip5p/"

in_json  = {"NSET":0,"ISTART":0,"NSAVC":0,"NAMNF":0,"DELTA":0,"title_size":0,
           "title_size2":0,"N":0,"cord_idx":0,"XYZ":[],
           "INIT":"0", #options.init
           "freeindexes":[],
           "psf":pth+"A_0.psf", #options.psf
           "dcds":"", #options.dcd
           "atms":"OH2", #options.atoms
           "outcsv":""} #options.outcsv

psf_file = xyzf.psf(in_json['psf'])
psf_file.read_psf_file()
OH = []
count = 0
a = in_json['atms'].split(',')
for i in psf_file.atmtyp1:
    if i ==  a[0]:
        OH.append(count)
    #elif i == a[1]:
    #    LJ.append(count)
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
    columns_dict[dcd_box+"_OH_"+dcdID] = OH_mZ #pd.Series(OH_mZ)
    columns_dict[dcd_box+"_LJ_"+dcdID] = LJ_mZ #pd.Series(LJ_mZ)
    columns_dict[dcd_box+"_OH_"+dcdID+"_L"] = OHs_Zl #pd.Series(OHs_Zl)
    columns_dict[dcd_box+"_LJ_"+dcdID+"_L"] = LJs_Zl #pd.Series(LJs_Zl)
DF = pd.DataFrame(columns_dict)
print("Getting Z distribution OH and monoatomic.")
print("Indexed scaled by multiplying by 0.002 and get nanosecons")
DF.index = (DF.index+1)*0.02
DF.index.name = 'ns'
print(DF.shape)
print(DF.columns," output ",self.outcsv,"\n")
DF.to_csv(self.outcsv)

