#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 18:53:21 2020

@author: noel
"""
import LocalMaxima.structure.XYZ_Formats as strct
import os
import sys
import LocalMaxima.utilities.CHARMM_Parser as CP
import LocalMaxima.energy.CA as CA
import LocalMaxima.utilities.AA_Info as AA

# analysis = ca.mmgbsa_ca_analysis(os.getcwd())
# __init__
selfdirpath = "/home/noel/Datasets_nobackup/Structures/1B27/brs_A41C_A83C/structures/"
components = {}
df_pair = ''
df_self = ''

# mmgbsa bind_ca
optionscrd = '0.crd'
optionspsf = '1b27_0_ad_a41c_a83c_box_ions_xplor.psf'
optionsgb = 'radi_z0.txt'
optionsgbz = 'radi_z500.txt'
optionssa = 'sasa_z0.txt'
optionssaz = 'sasa_z500.txt'
optionsmov = 'A,D'
optionsdis = 500

# analysis.mmgbsa_CA_bindingMatrix(options) 
crd_file = strct.crd(selfdirpath+'/'+optionscrd)
psf_file = strct.psf(selfdirpath+'/'+optionspsf)
psf_file.read_psf_file()
out_file = os.path.basename(optionscrd).split('.')[0]

params = CP.read_charmm_FF()

ent = 1

sa_z0 = [0]*len(crd_file.atmtyp1)
for i in open(selfdirpath+'/'+optionssa, 'r').read().split('\n'):
    ii = i.strip().split()
    if len(ii) > 0:  
        sa_z0[int(ii[0])-1] = float(ii[1])
sa_z500 = [0]*len(crd_file.atmtyp1)
for i in open(selfdirpath+'/'+optionssaz, 'r').read().split('\n'):
    ii = i.strip().split()
    if len(ii) > 0:
        sa_z500[int(ii[0])-1] = float(ii[1])    
#else:
gb_z0 = []
for i in open(selfdirpath+'/'+optionsgb, 'r').read().split('\n'):
    ii = i.strip().split()
    if len(ii) > 0:  
        gb_z0.append(ii[1])
gb_z500 = []
for i in open(selfdirpath+'/'+optionsgbz, 'r').read().split('\n'):
    ii = i.strip().split()
    if len(ii) > 0:  
        gb_z500.append(ii[1])
entity_id = []
if len(crd_file.atmtyp1) == len(psf_file.charge):
    ent = 0 
    for i in range(len(crd_file.chain_id)):
        if i == 0:
            ent = 1 
            chn = crd_file.chain_id[i]
        else:
            if chn != crd_file.chain_id[i]:
                ent += 1
                chn = crd_file.chain_id[i]
        entity_id.append(ent)

nuc = ['GUA','ADE','CYT','THY','URA']
pro = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HSE','HSD','HSP',
       'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
wat = ['TIP3']
ion = ['SOD','CLA']
component_count = 0
water_count = 0
ion_count = 0
carbonyl_Carbon_CTER_loc = []
component = []
component_id = []
for j in range(len(crd_file.atmtyp1)):
    resi = crd_file.aaid[j]
    atom = crd_file.atmtyp1[j]
    chan = crd_file.chain_id[j]
    if resi in nuc:
        if atom in params.AA[resi].atoms[0]:
            component.append('NUC1')
            component_id.append(component_count)
            if atom == "H5''":
                component_count += 1
        elif atom in params.AA[resi].atoms[1]:
            component.append('NUC2')
            component_id.append(component_count)
            if atom == "H1'":
                component_count += 1            
        elif atom in params.AA[resi].atoms[2]:
            component.append('NUC3')
            component_id.append(component_count)
            if atom in ['H8','H62','H42','H53','H5']:
                component_count += 1   
        elif atom in params.AA[resi].atoms[3]:
            component.append('NUC4')
            component_id.append(component_count)
            if atom == "H2'":
                component_count += 1   
        elif atom in params.AA[resi].atoms[4]:
            component.append('NUC5')
            component_id.append(component_count)
            if atom == "O3'":
                component_count += 1   
        else:
            # TODO check if when these atoms are found they are at the end of the component
            if atom == 'H5T':
                component.append('NUC1')
                component_id.append(component_count)
            elif atom == 'H3T':
                component.append('NUC5')
                component_id.append(component_count)
        component_count += 1
        print("IN NUC")
    elif resi in pro:
        if atom in params.AA[resi].atoms[0]:
            component.append('AMIN')
            component_id.append(component_count)
            if atom in ["HA","HA2"]:
                component_count += 1
        elif atom in params.AA[resi].atoms[-1]:
            component.append('CARB')
            component_id.append(component_count)
            if atom == "O":
                component_count += 1
        else:
            if atom in params.AA['ACE'].atoms[0]:
                component.append('ACE1T')
                component_id.append(component_count)
                if atom == "HY3":
                    component_count += 1
            elif atom in params.AA['ACE'].atoms[1]:
                component.append('ACE2T')
                component_id.append(component_count)
                if atom == "OY":
                    component_count += 1
            elif atom in params.AA['NTER'].atoms[0]:
                component.append('NTER')
                component_id.append(component_count)
                if atom == "HA":
                    component_count += 1
            elif atom in params.AA['CTER'].atoms[0]:
                component.append('CTER')
                component_id.append(component_count)
                if atom == "OT2":
                    carbonyl_Carbon_CTER_loc.append(j-2)
                    component_count += 1
            else:
                component.append('SIDE')
                component_id.append(component_count)
                if resi == "ALA":
                    if atom == "HB3":
                        component_count += 1
                elif resi == "ARG":
                    if atom == "HH22":
                        component_count += 1
                elif resi == "ASN":
                    if atom == "HD22":
                        component_count += 1
                elif resi == "ASP":
                    if atom == "OD2":
                        component_count += 1
                elif resi == "CYS":
                    if atom == "HG1":
                        component_count += 1
                elif resi == "GLU":
                    if atom == "OE2":
                        component_count += 1
                elif resi == "GLN":
                    if atom == "HE22":
                        component_count += 1
                elif resi == "GLY":
                    print("ERROR: GLY atom not recognized.")
                    sys.exit()
                elif resi == "HSE":
                    if atom == "HD2":
                        component_count += 1
                elif resi == "HSD":
                    if atom == "HD2":
                        component_count += 1
                elif resi == "HSP":
                    if atom == "HE1":
                        component_count += 1
                elif resi == "ILE":
                    if atom == "HD3":
                        component_count += 1
                elif resi == "LEU":
                    if atom == "HD23":
                        component_count += 1
                elif resi == "LYS":
                    if atom == "HZ3":
                        component_count += 1
                elif resi == "MET":
                    if atom == "HE3":
                        component_count += 1
                elif resi == "PHE":
                    if atom == "HE2":
                        component_count += 1
                elif resi == "PRO":
                    if atom == "HG2":
                        component_count += 1
                elif resi == "SER":
                    if atom == "HG1":
                        component_count += 1
                elif resi == "THR":
                    if atom == "HG23":
                        component_count += 1
                elif resi == "TRP":
                    if atom == "HH2":
                        component_count += 1
                elif resi == "TYR":
                    if atom == "HE1":
                        component_count += 1
                elif resi == "VAL":
                    if atom == "HG23":
                        component_count += 1
    elif resi in wat:
        if atom in params.AA['TIP3'].atoms[0]:
            component.append('WATR'+str(water_count))
            component_id.append(component_count)
            if atom == 'H2':
                component_count += 1
                water_count += 1
    elif resi in ion:
        if atom in params.AA['SOD'].atoms[0]:
            component.append('IONS'+str(ion_count))
            component_id.append(component_count)
            component_count += 1
            ion_count += 1
        if atom in params.AA['CLA'].atoms[0]:
            component.append('IONS'+str(ion_count))
            component_id.append(component_count)
            component_count += 1        
            ion_count += 1
    else:
        print('ERROR: Amino Acid '+resi+' not found in parameters. Exit now.')
        sys.exit(1)

for i in carbonyl_Carbon_CTER_loc:
    print(i,"Converting ",component[i]," to ",component[i+1])
    component[i] = "CTER"
    
    
print(len(crd_file.atmtyp1),len(component),len(component_id))
# 25001 25001 25001
# Component analysis check count for Barnase Barstar
A = 110         # AA 
B = 238 - A     # TIP3
D = 90          # AA
E = 154 - D     # TIP3
W = 7067
SOD = 22
CLA = 16
AAa = ((A+D)*3) - 15  # 15 GLYs
TC = AAa + B + E + W + SOD + CLA + 4  # 4 terminal components for 2 ACEs
print(TC)

out_file = 'out.txt'
dirpath = '/home/noel/Datasets_nobackup/Structures/1B27/brs_A41C_A83C/structures'
idxFile = open(dirpath+'/'+out_file, 'w')
for j in range(len(crd_file.atmtyp1)):
    idxFile.write(str(j)+' '+str(crd_file.aa[j])+' '+crd_file.aaid[j]+
                  ' '+crd_file.atmtyp1[j]+' '+crd_file.chain_id[j]+
                  ' '+str(entity_id[j])+' '+component[j]+' '+
                  str(component_id[j])+'\n')
idxFile.close()

mass = []
atmNum = []
atmtyp3 = []
epsilon = []
rmin_half = []
atminfo = []
for i in psf_file.atmtyp2:
    atmNum.append(params.am.MASS[i][0])
    mass.append(float(params.am.MASS[i][1]))
    atmtyp3.append(params.am.MASS[i][2])
    epsilon.append(params.NONBONDED[i][1])
    rmin_half.append(params.NONBONDED[i][2])
    atminfo.append(False)

# Component Indexing
comp_indx = []
comp = []
current_comp = ''
comp_end = []
comp_chain_id = []
comp_entity_id = []
comp_aaid = []
comp_aa = []
comp_cnt = 0
for i in range(len(crd_file.aa)):
    if i == 0:
        comp_indx.append(i)
        comp_chain_id.append(crd_file.chain_id[i])
        comp_entity_id.append(entity_id[i])
        comp_aaid.append(crd_file.aaid[i])
        comp_aa.append(crd_file.aa[i])
        comp.append(component[i])
        current_comp = component[i]
        comp_cnt = 0
    else:
        if component[i] != current_comp:
            comp_indx.append(i)
            comp_chain_id.append(crd_file.chain_id[i])
            comp_entity_id.append(entity_id[i])
            comp_aaid.append(crd_file.aaid[i])
            comp_aa.append(crd_file.aa[i])
            comp.append(component[i])
            current_comp = component[i]
            comp_end.append(comp_cnt)
            comp_cnt += 1
        else:
            comp_cnt += 1
comp_end.append(comp_cnt)
comp_indx_tuple = []
for i in range(len(comp_indx)):
    comp_indx_tuple.append((comp_chain_id[i],comp_entity_id[i],comp_aaid[i],
                            comp_aa[i],comp[i],comp_indx[i],comp_end[i]))

out_file = 'out2.txt'
dirpath = '/home/noel/Datasets_nobackup/Structures/1B27/brs_A41C_A83C/structures'
idxFile = open(dirpath+'/'+out_file, 'w')
for i in range(0,len(comp_indx_tuple)):
    idxFile.write(str(comp_indx_tuple[i])+'\n')
idxFile.close()

# Number of component-component calculations, or one per line, which is # of lines.
count = 0

for i in range(0,7885):
    for j in range(i+1,7886):
        count += 1
