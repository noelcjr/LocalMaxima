# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 23:01:32 2017

The classes here correspond to MMGBSA_CA_L.py in Entropy Maxima.

@author: noel
"""
import os
import sys
import LocalMaxima.utilities.CHARMM_Parser as CP
import LocalMaxima.utilities.AA_Info as AA
#import em.code.Energy_Functions as EF
#import em.code.Super_Structures as SS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from LocalMaxima.structure.XYZ_Formats import crd,psf

class hydration_shell_from_gb(object):
    def __init__(self,crd_path, gb0_path, gb500_path):
        self.coor = crd(crd_path)
        self.coor.read_crd_file()
        self.gb0_list = []
        for i in open(gb0_path, 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0:
                self.gb0_list.append(ii[1])
        self.gb500_list = []
        for i in open(gb500_path, 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0:
                self.gb500_list.append(ii[1])
        print("Check number of entries in each file:")
        print("     crd = "+str(self.coor.xyz.shape[0]))
        print("     gb0 = "+str(len(self.gb0_list)))
        print("   gb500 = "+str(len(self.gb500_list)))
        
    def map_waters_ions(self):
        crd_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B27/brs_A41C_A83C/structures/'
        params = CP.read_charmm_FF()
        for k in range(2):
            coor = crd(crd_path+'0.crd')
            coor.read_crd_file()
            gb0_list = []
            for i in open(crd_path+'radi_z0.txt', 'r').read().split('\n'):
                ii = i.strip().split()
                if len(ii) > 0:
                    gb0_list.append(float(ii[1]))
            gb500_list = []
            for i in open(crd_path+'radi_z500.txt', 'r').read().split('\n'):
                ii = i.strip().split()
                if len(ii) > 0:
                    gb500_list.append(float(ii[1]))
            if coor.xyz.shape[0] == len(gb0_list):
                if len(gb0_list) == len(gb500_list):
                    #print('yes')
                    # 1. Map Waters and Ions
                    waters = []
                    ions = []
                    for i in coor.crd_contents:
                        if i[2] == 'TIP3':
                            if i[4] == 'OH2':
                                waters.append(i[0])
                        if i[2] in ['SOD','CLA']:
                            ions.append(i[0])
                    # 2. Confir that water go before or after Ions.
                    EE_K = -327.90
                    EE_K2 = -163.95
                    hyd_gb = pd.DataFrame(columns=['atm_num','res_num','mol_typ','atm_type','chn','DGB'])
                    index = 0
                    if (waters[-1]+3) == ions[0]:
                        qH = params.AA['TIP3'].atom_chrg['H1']
                        qO = params.AA['TIP3'].atom_chrg['OH2']
                        EE_KqOqH = EE_K*qO*qH
                        EE_KqHqH = EE_K*qH*qH
                        for j in waters:
                            i = j - 1
                            GB = EE_K2*((qO*qO/gb0_list[i])+(qH*qH/gb0_list[i+1])+(qH*qH/gb0_list[i+2]))
                            GB_Z = EE_K2*((qO*qO/gb500_list[i])+(qH*qH/gb500_list[i+1])+(qH*qH/gb500_list[i+2]))
                            
                            d = np.power(np.linalg.norm(coor.xyz[i]-coor.xyz[i+1]),2)
                            temp3 = gb0_list[i]*gb0_list[i+1]
                            temp3_Z = gb500_list[i]*gb500_list[i+1]
                            GB += EE_KqOqH/np.sqrt(d + (temp3*np.exp((-1*d)/(4*temp3))))
                            GB_Z += EE_KqOqH/np.sqrt(d + (temp3_Z*np.exp((-1*d)/(4*temp3_Z))))
                            
                            d = np.power(np.linalg.norm(coor.xyz[i]-coor.xyz[i+2]),2)
                            temp3 = gb0_list[i]*gb0_list[i+2]
                            temp3_Z = gb500_list[i]*gb500_list[i+2]
                            GB += EE_KqOqH/np.sqrt(d + (temp3*np.exp((-1*d)/(4*temp3))))
                            GB_Z += EE_KqOqH/np.sqrt(d + (temp3_Z*np.exp((-1*d)/(4*temp3_Z))))
                            
                            d = np.power(np.linalg.norm(coor.xyz[i+1]-coor.xyz[i+2]),2)
                            temp3 = gb0_list[i+1]*gb0_list[i+2]
                            temp3_Z = gb500_list[i+1]*gb500_list[i+2]
                            GB += EE_KqHqH/np.sqrt(d + (temp3*np.exp((-1*d)/(4*temp3))))
                            GB_Z += EE_KqHqH/np.sqrt(d + (temp3_Z*np.exp((-1*d)/(4*temp3_Z))))
                            hyd_gb.loc[index] = [coor.crd_contents[i][0],coor.crd_contents[i][1],coor.crd_contents[i][2],
                                       coor.crd_contents[i][4],coor.crd_contents[i][9],(GB-GB_Z)]
                            index += 1
                        for j in ions:
                            i = j - 1
                            GB = EE_K2/gb0_list[i]
                            GB_Z = EE_K2/gb500_list[i]
                            hyd_gb.loc[index] = [coor.crd_contents[i][0],coor.crd_contents[i][1],coor.crd_contents[i][2],
                                                 coor.crd_contents[i][4],coor.crd_contents[i][9],(GB-GB_Z)]
                            index += 1
                        segments = [float(xx)/10 for xx in range(10)]
                        histo = plt.hist(hyd_gb.DGB,bins=segments)
                        solvent_indexes = {}
                        init = histo[1][0]
                        for i in histo[1][1:]:
                            solvent_indexes[i] = list(hyd_gb.atm_num[(hyd_gb.DGB >= init) & (hyd_gb.DGB < i)])
                            init = i
                elif (ions[-1]+1) == waters[0]:
                    #TODO Not common 
                    print("Error: Ions before waters in structure")
                    sys.exit(1)
                
            else:
                print("Error: gb0 and gb500 files have different number of elements.")
                sys.exit(1)
        else:
            print("Error: xyz and gb500 files have different number of elements.")
            sys.exit(1)
        for h in [float(xx)/10 for xx in range(1,10)]:
            new_crd = []
            for i in coor.crd_contents:
                if i[2] == 'TIP3':
                    if i[0] in solvent_indexes[h]:
                        new_crd.append(i[0])
                        new_crd.append(i[0]+1)
                        new_crd.append(i[0]+2)
                elif i[2] in ['SOD','CLA']:
                    new_crd.append(i[0])
                else:
                    new_crd.append(i[0])
            new_crd2 = []
            for i in new_crd:
                new_crd2.append(coor.crd_contents[i-1])
            pdb_converted_lines = []
            for i in new_crd2:
                if len(i[4]) == 4:
                    if len(i[2]) == 4:
                        line = '{:6}{:>5} {:4} {:4}{:1}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>10}'
                    else:
                        line = '{:6}{:>5} {:4} {:3} {:1}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>10}'
                else:
                    if len(i[2]) == 4:
                        line = '{:6}{:>5}  {:3} {:4}{:1}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>10}'
                    else:
                        line = '{:6}{:>5}  {:3} {:3} {:1}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>10}'
                pdb_converted_lines.append(line.format('ATOM',i[0],i[4],i[2],i[9],int(i[11]),float(i[5]),\
                                                       float(i[6]),float(i[7]),float('0'),float('0'),'IONS',''))
            outFile = open(crd_path+'0_'+str(h)+'.pdb', 'w')
            for i in pdb_converted_lines:
                outFile.write(i + '\n')
            outFile.close()
            
    def get_conserved_waters(self):
        print("Getting Conserved waters.")
# This class manages the search for a linker or loop in search of some outcome. 
# It could be moved later if it is too large to be here. It uses 
# mmgbsa_ca_analysis.pyc and bash scripts to automate the search for a 
# structure. Start by placinf in this class Component_Analysis_linkers.py
class mmgbsa_ca_searcher(object):
    def __init__(self):
        pass

class mmgbsa_ca_analysis(object):
    def __init__(self,oscwd):
        self.dirpath = oscwd
        self.component = []
        self.component_id = []
        self.comp_indx = []
        self.df_pair = ''
        self.df_self = ''
        
    def rank_mmgbsa(self):
        folder = '/home/noel/Projects/Protein_design/EntropyMaxima/examples/CA_analysis/analyze_CA/'
        file_name = '2hiu_1rr_min_mmgbsa.txt'
        mmgbsa = pd.read_csv(folder+file_name)
        mmgbsa_v =mmgbsa.values
        mmgbsa_dict = {}
        for i in range(mmgbsa_v.shape[0]):
            if mmgbsa_v[i][0] in mmgbsa_dict:
                mmgbsa_dict[mmgbsa_v[i][0]] = mmgbsa_dict[mmgbsa_v[i][0]] + mmgbsa_v[i][2:]
            else:
                mmgbsa_dict[mmgbsa_v[i][0]] = mmgbsa_v[i][2:]
            if mmgbsa_v[i][1] in mmgbsa_dict:
                mmgbsa_dict[mmgbsa_v[i][1]] = mmgbsa_dict[mmgbsa_v[i][1]] + mmgbsa_v[i][2:]
            else:
                mmgbsa_dict[mmgbsa_v[i][1]] = mmgbsa_v[i][2:]
                
    # TODO https://en.wikipedia.org/wiki/Connected-component_labeling
    def connected_component_labeling(self):
        pass
    
    def mmgbsa_raw_summary(self,options):
        # TODO: NOT USED. For now only working with delta free energies. 
        #    In the future I plan to use direct free
        #    energies of components in search of paterns or new insights. 
        #    Direct free energies need 1-4 interactions to
        #    be more meaningful, and the 1-4 interactions code has a bug!
        if not options.range:
            options.range = '(-10,10,200):(-5,5,300)'
        filename_prefix_for_output = options.pair2.split('.')[0][0:-5]
        fullpath = self.dirpath+'/'+filename_prefix_for_output
        self.df_pair = pd.read_csv(self.dirpath+'/'+options.pair2)
        self.df_self = pd.read_csv(self.dirpath+'/'+options.self2)
        self.fix_cols_components()
        r1 = 'pos_ene(>+'+str(options.cut)+'kcal/mol)'
        r2 = 'neg_ene(<-'+str(options.cut)+'kcal/mol)'
        rows = ['min_val','max_val','mean','std','total_components','nan',r1,\
                'pos_eng(%)',r2,'neg_eng(%)','included',\
                'included(%)','excluded','excluded(%)']
        self_energy_terms = 's'+self.df_self.columns[2:]
        pair_energy_terms = 'p'+self.df_pair.columns[2:]
        summary = pd.DataFrame(index=rows,columns=pair_energy_terms.union(self_energy_terms))
        for i in pair_energy_terms:
            summary[i] = summary[i].astype(float)
        for i in self_energy_terms:
            summary[i] = summary[i].astype(float)
        leng_self = len(self.df_self)
        for i in self_energy_terms:
            d = self.df_self.describe()[i[1:]]
            temp = []
            temp.append(d['min'])
            temp.append(d['max'])
            temp.append(d['mean'])
            temp.append(d['std'])
            temp.append(leng_self)
            np_self = self.df_self.loc[:,i[1:]].values
            temp.append(np.isnan(np_self).sum())
            neg_eng = np.where(np_self < -options.cut)
            pos_eng = np.where(np_self > options.cut)
            temp.append(len(pos_eng[0]))
            temp.append(len(pos_eng[0])/float(leng_self))
            temp.append(len(neg_eng[0]))
            temp.append(len(neg_eng[0])/float(leng_self))
            included = 100*(len(pos_eng[0])+len(neg_eng[0]))/float(leng_self)
            excluded = 100-included
            pos_neg_tot = len(pos_eng[0])+len(neg_eng[0])
            temp.append(pos_neg_tot)
            temp.append(included)
            temp.append(leng_self-pos_neg_tot)
            temp.append(excluded)
            summary[i] = np.array(temp)
        leng_pair = len(self.df_pair)
        for i in pair_energy_terms:
            d = self.df_pair.describe()[i[1:]]
            temp = []
            temp.append(d['min'])
            temp.append(d['max'])
            temp.append(d['mean'])
            temp.append(d['std'])
            temp.append(leng_pair)
            np_pairs = self.df_pair.loc[:,i[1:]].values
            temp.append(np.isnan(np_pairs).sum())
            neg_eng = np.where(np_pairs < -options.cut)
            pos_eng = np.where(np_pairs > options.cut)
            temp.append(len(pos_eng[0]))
            temp.append(len(pos_eng[0])/float(leng_pair))
            temp.append(len(neg_eng[0]))
            temp.append(len(neg_eng[0])/float(leng_pair))
            included = 100*(len(pos_eng[0])+len(neg_eng[0]))/float(leng_pair)
            excluded = 100-included
            pos_neg_tot = len(pos_eng[0])+len(neg_eng[0])
            temp.append(pos_neg_tot)
            temp.append(included)
            temp.append(leng_pair-pos_neg_tot)
            temp.append(excluded)
            summary[i] = np.array(temp)
        print(summary)
        summary.to_csv(fullpath+'_rawsum.csv')
        fig, axes = plt.subplots(nrows=2, ncols=2)
        ax0, ax1, ax2, ax3 = axes.flatten()
        pair_labels = ['pGBE', 'pGBE_Z', 'pEE','pEE_Z','pVWE','pVWE_Z']
        self_labels = ['sGBE', 'sGBE_Z', 'sSAE','sSAE_Z']
        x = self.df_pair.loc[:,self.df_pair.columns[2:]]
        y = self.df_self.loc[:,self.df_self.columns[2:]]
        tuples = []
        for i in options.range.split(':'):
            tuples.append(tuple(i[1:-1].split(',')))
        if len(tuples) != 2:
            print("ERROR: The string parameters for histograms must be two tuples separated by columns. For example:\
            \"(-10,10,200):(-5,5,300)\". First two values of the tuples can be real numbers corresponding to a minimum \
            and maimum range for a histogram. The las value is an integer corresponding to the number of bins.")
            sys.exit(1)
        # First histogram's ranges and numbers of bins are selcted automatically to contain the minimum and maximum
        # values of the distribution and 100 bins.
        n_bins = 100
        ax0.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax0.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax0.legend(prop={'size': 10})
        ax0.set_title('Full Energy Range. Bins=100')
        # Second histogram uses the options.cut as limits on the x-axis and 10000 bins
        n_bins = 10000
        ax1.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax1.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax1.legend(prop={'size': 10})
        ax1.set_title("Energy Range within options.cut. Bins=10000' ")
        ax1.set_xlim((-options.cut, options.cut))
        # Third and Fourth Histograms are defulted two (-10,10,200):(-5,5,300), but can be changed by user.
        n_bins = int(tuples[0][2])
        ax2.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax2.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax2.legend(prop={'size': 10})
        ax2.set_title('Energy Range '+str(tuples[0][0:2])+". Bins="+str(tuples[0][2]))
        ax2.set_xlim((float(tuples[0][0]), float(tuples[0][1])))
        #axes.set_ylim((-10, 10))
        n_bins = int(tuples[1][2])
        ax3.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax3.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax3.legend(prop={'size': 10})
        ax3.set_title('Energy Range '+str(tuples[1][0:2])+". Bins="+str(tuples[1][2]))
        ax3.set_xlim((float(tuples[1][0]), float(tuples[1][1])))
        #axes.set_ylim((-10, 10))
        fig.set_size_inches(12.0, 12.0)
        fig.savefig(self.dirpath+'/'+filename_prefix_for_output+'_rawhis.png')
        
    def mmgbsa_binding(self, options):
        filename_prefix_for_output = options.pair1.split('.')[0][0:-5]
        fullpath = self.dirpath+'/'+filename_prefix_for_output
        self.df_pair = pd.read_csv(self.dirpath+'/'+options.pair1)
        self.df_self = pd.read_csv(self.dirpath+'/'+options.self1)
        self.fix_cols_components()
        np_self = self.df_self.loc[:,self.df_self.columns].values
        np_pair = self.df_pair.loc[:,self.df_pair.columns].values
        # Binding between the protein and DNA and RNA
        MMGBSA = np.zeros((np_pair.shape[0],12))
        for i in range(np_pair.shape[0]):
            MMGBSA[i][0] = int(np_pair[i][0])                                           # Component 1
            MMGBSA[i][1] = int(np_pair[i][1])                                           # Component 2
            self1 = np_self[int(MMGBSA[i][0])][2] - np_self[int(MMGBSA[i][0])][3]       # Self 1 GB
            self2 = np_self[int(MMGBSA[i][1])][2] - np_self[int(MMGBSA[i][1])][3]       # Self 2 GB
            MMGBSA[i][2] = self1 + self2 + (np_pair[i][2] - np_pair[i][3])              # pair GB
            MMGBSA[i][3] = np_pair[i][4] - np_pair[i][5]                                # pair EE
            MMGBSA[i][4] = np_pair[i][6] - np_pair[i][7]                                # pair VDW folled by SA
            MMGBSA[i][5] = (np_self[int(MMGBSA[i][0])][4]-np_self[int(MMGBSA[i][0])][5])+(np_self[int(MMGBSA[i][1])][4]-np_self[int(MMGBSA[i][1])][5])
            MMGBSA[i][6] = MMGBSA[i][2] + MMGBSA[i][3]                                  # Polar
            MMGBSA[i][7] = MMGBSA[i][4] + MMGBSA[i][5]                                  # Non-Polar
            MMGBSA[i][8] = MMGBSA[i][6] + MMGBSA[i][7]                                  # MMGBSA
        # First rank unfiltered. Output results and decide later/
        mmgbsa = []
        for i in range(MMGBSA.shape[0]):
            mmgbsa.append([self.components[int(MMGBSA[i][0])],self.components[int(MMGBSA[i][1])],MMGBSA[i][2],\
                           MMGBSA[i][3],MMGBSA[i][4],MMGBSA[i][5],MMGBSA[i][6],MMGBSA[i][7],MMGBSA[i][8]])
        outFile = open(fullpath+'_mmgbsa.txt', 'w')
        outFile.write('comp1,comp2,dGB,dEE,dVDW,dSA,Polar,Nonpolar,MMGBSA\n')
        for i in mmgbsa:
            temp = ''
            temp += i[0]+','+i[1]+','
            for j in i[2:]:
                # TODO add formating to other outputs to a text file. 
                # Reduce significant figures.
                temp += "{:.3e}".format(j)+','
            outFile.write(temp[0:-1]+'\n')
        outFile.close()

    def load_data(self,options):
        """ 
        This function loads CRD, PSF and parameter files in order to do
        Component analysis. It creates a component list with the indexes to
        atoms for each component in the protein.
        """ 
        self.crd_file = crd(self.dirpath+'/'+options.xyz)
        # This PSF was written by GenInputs/get_gb_sa.
        self.psf_file = psf(self.dirpath+'/temp_xplo.psf')
        self.out_file = os.path.basename(options.xyz).split('.')[0]
        ent = 1
        self.sa_z0 = [0]*len(self.crd_file.atmtyp1)
        for i in open(self.dirpath+'/sasa_z0.txt', 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0: 
                self.sa_z0[int(ii[0])-1] = float(ii[1])
        self.sa_z500 = [0]*len(self.crd_file.atmtyp1)
        for i in open(self.dirpath+'/sasa_z500.txt', 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0:
                self.sa_z500[int(ii[0])-1] = float(ii[1])   
        #else:
        self.gb_z0 = []
        for i in open(self.dirpath+'/gb_z0.txt', 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0: 
                self.gb_z0.append(float(ii[1]))
        self.gb_z500 = []
        for i in open(self.dirpath+'/gb_z500.txt', 'r').read().split('\n'):
            ii = i.strip().split()
            if len(ii) > 0: 
                self.gb_z500.append(float(ii[1]))
        self.entity_id = []
        if len(self.crd_file.atmtyp1) == len(self.psf_file.charge):
            ent = 0
            for i in range(len(self.crd_file.chain_id)):
                if i == 0:
                    ent = 1
                    chn = self.crd_file.chain_id[i]
                else:
                    if chn != self.crd_file.chain_id[i]:
                        ent += 1
                        chn = self.crd_file.chain_id[i]
                self.entity_id.append(ent)
        self.mass = []
        self.atmNum = []
        self.atmtyp3 = []
        self.epsilon = []
        self.rmin_half = []
        self.atminfo = []
        for i in self.psf_file.atmtyp2:
            self.atmNum.append(self.params.am.MASS[i][0])
            self.mass.append(float(self.params.am.MASS[i][1]))
            self.atmtyp3.append(self.params.am.MASS[i][2])
            self.epsilon.append(float(self.params.NONBONDED[i][1]))
            self.rmin_half.append(float(self.params.NONBONDED[i][2]))
            self.atminfo.append(False)
            
    def component_indexing(self):
        ## After reading files, Generate component Indexes #########
        nuc = ['GUA','ADE','CYT','THY','URA']
        pro = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HSE','HSD','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        wat = ['TIP3']
        ion = ['SOD','CLA']
        component_count = 0
        water_count = 0
        ion_count = 0
        carbonyl_Carbon_CTER_loc = []
        self.component = []
        self.component_id = []
        for j in range(len(self.crd_file.atmtyp1)):
            resi = self.crd_file.aaid[j]
            atom = self.crd_file.atmtyp1[j]
            if resi in nuc:
                if atom in self.params.AA[resi].atoms[0]:
                    self.component.append('NUC1')
                    self.component_id.append(component_count)
                    if atom == "H5''":
                        component_count += 1
                elif atom in self.params.AA[resi].atoms[1]:
                    self.component.append('NUC2')
                    self.component_id.append(component_count)
                    if atom == "H1'":
                        component_count += 1            
                elif atom in self.params.AA[resi].atoms[2]:
                    self.component.append('NUC3')
                    self.component_id.append(component_count)
                    if atom in ['H8','H62','H42','H53','H5']:
                        component_count += 1   
                elif atom in self.params.AA[resi].atoms[3]:
                    self.component.append('NUC4')
                    self.component_id.append(component_count)
                    if atom == "H2'":
                        component_count += 1   
                elif atom in self.params.AA[resi].atoms[4]:
                    self.component.append('NUC5')
                    self.component_id.append(component_count)
                    if atom == "O3'":
                        component_count += 1   
                else:
                    # TODO check if when these atoms are found they are at the
                    # end of the component
                    if atom == 'H5T':
                        self.component.append('NUC1')
                        self.component_id.append(component_count)
                    elif atom == 'H3T':
                        self.component.append('NUC5')
                        self.component_id.append(component_count)
                component_count += 1
            elif resi in pro:
                if atom in self.params.AA[resi].atoms[0]:
                    self.component.append('AMIN')
                    self.component_id.append(component_count)
                    if atom in ["HA","HA2"]:
                        component_count += 1
                elif atom in self.params.AA[resi].atoms[-1]:
                    self.component.append('CARB')
                    self.component_id.append(component_count)
                    if atom == "O":
                        component_count += 1
                else:
                    if atom in self.params.AA['ACE'].atoms[0]:
                        self.component.append('ACE1T')
                        self.component_id.append(component_count)
                        if atom == "HY3":
                            component_count += 1
                    elif atom in self.params.AA['ACE'].atoms[1]:
                        self.component.append('ACE2T')
                        self.component_id.append(component_count)
                        if atom == "OY":
                            component_count += 1
                    elif atom in self.params.AA['NTER'].atoms[0]:
                        self.component.append('NTER')
                        self.component_id.append(component_count)
                        if atom == "HA":
                            component_count += 1
                    elif atom in self.params.AA['CTER'].atoms[0]:
                        self.component.append('CTER')
                        self.component_id.append(component_count)
                        if atom == "OT2":
                            carbonyl_Carbon_CTER_loc.append(j-2)
                            component_count += 1
                    else:
                        self.component.append('SIDE')
                        self.component_id.append(component_count)
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
                if atom in self.params.AA['TIP3'].atoms[0]:
                    self.component.append('WATR'+str(water_count))
                    self.component_id.append(component_count)
                    if atom == 'H2':
                        component_count += 1
                        water_count += 1
            elif resi in ion:
                if atom in self.params.AA['SOD'].atoms[0]:
                    self.component.append('IONS'+str(ion_count))
                    self.component_id.append(component_count)
                    component_count += 1
                    ion_count += 1
                if atom in self.params.AA['CLA'].atoms[0]:
                    self.component.append('IONS'+str(ion_count))
                    self.component_id.append(component_count)
                    component_count += 1        
                    ion_count += 1
            else:
                print('ERROR: Amino Acid '+resi+' not found in parameters.')
                print('Exiting now.')
                sys.exit(1)
        for i in carbonyl_Carbon_CTER_loc:
            #print(i,"Converting ",self.component[i]," to ",self.component[i+1])
            self.component[i] = "CTER"
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
        for i in range(len(self.crd_file.aa)):
            if i == 0:
                comp_indx.append(i)
                comp_chain_id.append(self.crd_file.chain_id[i])
                comp_entity_id.append(self.entity_id[i])
                comp_aaid.append(self.crd_file.aaid[i])
                comp_aa.append(self.crd_file.aa[i])
                comp.append(self.component[i])
                current_comp = self.component[i]
                comp_cnt = 0
            else:
                if self.component[i] != current_comp:
                    comp_indx.append(i)
                    comp_chain_id.append(self.crd_file.chain_id[i])
                    comp_entity_id.append(self.entity_id[i])
                    comp_aaid.append(self.crd_file.aaid[i])
                    comp_aa.append(self.crd_file.aa[i])
                    comp.append(self.component[i])
                    current_comp = self.component[i]
                    comp_end.append(comp_cnt)
                    comp_cnt += 1
                else:
                    comp_cnt += 1
        comp_end.append(comp_cnt)
        for i in range(len(comp_indx)):
            self.comp_indx.append((comp_chain_id[i],comp_entity_id[i],comp_aaid[i],
                                   comp_aa[i],comp[i],comp_indx[i],comp_end[i]))

    def mmgbsa_CA_bindingMatrix(self,options):
        self.params = CP.read_charmm_FF()
        self.load_data(options)
        self.component_indexing()
        # Gen A vector of displacement.
        xyz_d = []
        move_chains = options.mov.split(',')
        for i in range(len(self.crd_file.xyz)):
            if self.crd_file.chain_id[i] in move_chains:
                xyz_d.append(self.crd_file.xyz[i]+np.array((0,0,options.trans)))
            else:
                xyz_d.append(self.crd_file.xyz[i])
        EE_K = -327.90
        EE_K2 = -163.95
        #Surface tension coefficient     (sgamma) =   0.010 [kcal/mol/Angs**2]
        SA_K = 0.010
        l = len(self.comp_indx)
        idxFile = open(self.dirpath+'/'+self.out_file+'_MMGBSA.idx', 'w')
        outFile = open(self.dirpath+'/'+self.out_file+'_MMGBSA.out', 'w')
        tempGB = 0
        tempGB_Z = 0
        tempSA = 0
        tempSA_Z = 0
        tempEE = 0
        tempEE_Z = 0
        tempVDW = 0
        tempVDW_Z = 0
        pairGB = 0
        pairGB_Z = 0
        DGB = 0.0
        DGB_Z = 0.0
        DEE = 0.0
        DEE_Z = 0.0
        DVW = 0.0
        DVW_Z = 0.0
        DSA = 0.0
        DSA_Z = 0.0
        for i in range(0,l):
            idxFile.write(str(i)+" "+str(self.comp_indx[i][2])+" "+\
                          self.comp_indx[i][0]+" "+str(self.comp_indx[i][3])+\
                          " "+self.comp_indx[i][4]+" "+str(self.comp_indx[i][5])+\
                          " "+str(self.comp_indx[i][6])+"\n")
            for k in range(self.comp_indx[i][5],self.comp_indx[i][6]):
                temp5 = EE_K2*self.psf_file.charge[k]*self.psf_file.charge[k]/self.gb_z0[k]
                tempGB += temp5
                temp5_Z = EE_K2*self.psf_file.charge[k]*self.psf_file.charge[k]/self.gb_z500[k]
                tempGB_Z += temp5_Z
                tempSA += SA_K*self.sa_z0[k]
                tempSA_Z += SA_K*self.sa_z500[k]
                for kk in range(k+1,self.comp_indx[i][6]+1):
                    last = kk
                    d1 = np.linalg.norm(self.crd_file.xyz[k]-self.crd_file.xyz[kk])
                    r2 = d1*d1
                    d1_Z = np.linalg.norm(self.crd_file.xyz[k]-self.crd_file.xyz[kk])
                    r2_Z = d1_Z*d1_Z
                    temp3 = self.gb_z0[k]*self.gb_z0[kk]
                    temp4 = (-1*r2)/(4*temp3)
                    temp3_Z = self.gb_z500[k]*self.gb_z500[kk]
                    temp4_Z = (-1*r2_Z)/(4*temp3_Z)
                    denm = np.sqrt(r2 + (temp3*np.exp(temp4)))
                    denm_Z = np.sqrt(r2_Z + (temp3_Z*np.exp(temp4_Z)))
                    Kqq = EE_K*self.psf_file.charge[k]*self.psf_file.charge[kk]
                    temp5 = Kqq/denm
                    temp5_Z = Kqq/denm_Z
                    tempGB += temp5
                    tempGB_Z += temp5_Z
            tempSA += SA_K*self.sa_z0[last]
            tempSA_Z += SA_K*self.sa_z500[last]
            temp5 = EE_K2*self.psf_file.charge[last]*self.psf_file.charge[last]/self.gb_z0[last]
            tempGB += temp5
            temp5_Z = EE_K2*self.psf_file.charge[last]*self.psf_file.charge[last]/self.gb_z500[last]
            tempGB_Z += temp5_Z
            DGB += tempGB
            DGB_Z += tempGB_Z
            DSA += tempSA
            DSA_Z += tempSA_Z
            outFile.write(str(i)+" "+str(i)+" "+str(tempGB-tempGB_Z)+" "+str(tempSA-tempSA_Z)+"\n")
            tempSA = 0.0
            tempSA_Z = 0.0
            for j in range(i+1,l):
                for m in range(self.comp_indx[i][5],self.comp_indx[i][6]+1):
                    for n in range(self.comp_indx[j][5],self.comp_indx[j][6]+1):
                        r = np.linalg.norm(self.crd_file.xyz[m]-self.crd_file.xyz[n])
                        r2 = r*r
                        r_Z = np.linalg.norm(xyz_d[m]-xyz_d[n])
                        r2_Z = r_Z*r_Z
                        temp = (332.06*self.psf_file.charge[m]*self.psf_file.charge[n])/r
                        temp_Z = (332.06*self.psf_file.charge[m]*self.psf_file.charge[n])/r_Z
                        tempEE += temp
                        tempEE_Z += temp_Z
                        Eps = np.sqrt(self.epsilon[m]*self.epsilon[n])
                        Rmin = self.rmin_half[m] + self.rmin_half[n]
                        A = Rmin/r
                        A_Z = Rmin/r_Z
                        A2 = A*A
                        A2_Z = A_Z*A_Z
                        A6 = A2*A2*A2
                        A6_Z = A2_Z*A2_Z*A2_Z
                        A12 = A6*A6
                        A12_Z = A6_Z*A6_Z
                        tempVDW += Eps*(A12-(2*A6))
                        tempVDW_Z += Eps*(A12_Z-(2*A6_Z))
                        temp3 = self.gb_z0[m]*self.gb_z0[n]
                        temp4 = (-1*r2)/(4*temp3)
                        temp3_Z = self.gb_z500[m]*self.gb_z500[n]
                        temp4_Z = (-1*r2_Z)/(4*temp3_Z)
                        denm = np.sqrt(r2 + (temp3*np.exp(temp4)))
                        denm_Z = np.sqrt(r2_Z + (temp3_Z*np.exp(temp4_Z)))
                        temp5 = (EE_K*self.psf_file.charge[m]*self.psf_file.charge[n])/denm
                        temp5_Z = (EE_K*self.psf_file.charge[m]*self.psf_file.charge[n])/denm_Z
                        pairGB += temp5
                        pairGB_Z += temp5_Z
                DGB += pairGB
                DGB_Z += pairGB_Z
                DEE += tempEE
                DEE_Z += tempEE_Z
                DVW += tempVDW
                DVW_Z += tempVDW_Z
                outFile.write(str(i)+" "+str(j)+" "+str(tempEE-tempEE_Z)+" "+str(tempVDW-tempVDW_Z)+" "+\
                              str(pairGB-pairGB_Z)+"\n")
                pairGB = 0
                pairGB_Z = 0
                tempEE = 0
                tempEE_Z = 0
                tempVDW = 0
                tempVDW_Z = 0
            tempGB = 0
            tempGB_Z = 0
        os.system("rm gbsa_z* sasa_z* gb_z* temp_xplo.psf")
        
    def mmgbsa_CA_bindingMatrix_14exclusions(self,options):
        # TODO: Not working. Some problems with 1-4exclusions. Some missing exclusion? not enough testing.
        #ef = EF.Energy_Functions(insu, idx_ss)
        #ef.get_exclusions_1_4()
        #count = 0
        #exclu = []
        #for i in ef.exclusions_1_4:
        #    for j in ef.exclusions_1_4[i]:
        #        if not len(ef.exclusions_1_4[i][j]) == 0:
        #            for k in ef.exclusions_1_4[i][j]:
        #                exclu.append(k)
        #                count += 1
        ELECB = 0.0
        ELECU = 0.0
        count = 0
        #for k in range(ch_AB.shape[0]-1):
        #    for kk in range(k+1,ch_AB.shape[0]):
        #        if not (k,kk) in exclu:
        #            Eqq = 332.06*ch_AB[k][0]*ch_AB[kk][0]
        #            d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
        #            ELECB += Eqq/d1
        #            d1_z = np.linalg.norm(np.array((ch_AB[k][3],ch_AB[k][4],ch_AB[k][8]))-\
        #                                  np.array((ch_AB[kk][3],ch_AB[kk][4],ch_AB[kk][8])))
        #            ELECU += Eqq/d1_z
        #            count += 1
        print("Pair calculations=",count)
        print("ELECB=",ELECB)
        print("ELECU=",ELECU)
    def fix_cols_components(self):
        for i in self.df_self.columns:
            self.df_self.rename(columns={i:i.replace(' ','')}, inplace=True)
        for i in self.df_pair.columns:
            self.df_pair.rename(columns={i:i.replace(' ','')}, inplace=True)
        for i in self.df_self.index:
            self.components[self.df_self.loc[i,'i']] = self.df_self.loc[i,'Chain_entity_aa_aaid_comp']
