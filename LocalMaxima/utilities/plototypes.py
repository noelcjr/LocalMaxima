#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 18:23:50 2021

@author: noel
"""
import sys
import numpy as np
import Bio.PDB as PDB
import struct
import matplotlib.pyplot as plt
import os
import pandas as pd

class plototypes(object):
    def __init__(self,options,plotType):
        self.multiPlotDIM = {1:(1,1),2:(1,2),3:(2,2),4:(2,2),5:(2,3),6:(2,3),
                             7:(3,3),8:(3,3),9:(3,3),10:(3,4),11:(3,4),12:(3,4),
                             30:(5,6),36:(6,6)} 
        self.plot_out = options.out
        if plotType == 'zmplot':
            self.plot_title = os.path.basename(options.input).split('.')[0]
            #if options.out:
            #    self.plot_out = os.path.dirname(options.input)+options.out
            #else:
            #    self.plot_out = os.path.dirname(options.input)+self.plot_title+".png"
            self.DF = pd.read_csv(options.input,index_col=0)
            #self.xy_axis = options.xy
            print(self.plot_title)
            print("Output:",self.plot_out)
            print(self.DF.columns)
            print(self.DF.shape, self.DF.shape[1])
        elif plotType == 'simout':
            pass
            #if options.out:
            #    self.plot_out = options.out
            #else:
            #    self.plot_out = os.getcwd().split('/')[-2]+"_"+os.getcwd().split('/')[-1]+"_out.png"
    def plot_sim_outputs(self,options):
        clms = []
        for i in options.columns.split(','):
            clms.append(i)
        #count_input_groups = 0
        selected_cols = ['PEvdw', 'PEee','P']
        DF_list = []
        self.plot_title = os.getcwd().split('/')[-2]+" "+os.getcwd().split('/')[-1]
        subplots = []
        for i in options.input:
            subplots.append(i.split('/')[0]+" "+i.split('/')[1].split("_")[0])
            out_path = os.path.dirname(i)
            #out_box = out_path.split('/')[-1]
            individual_files = os.path.basename(i).split(",")
            out_files = [out_path+"/"+jj for jj in individual_files]
            #outID = self.check_file_names(individual_files)
            df_list = []
            for i in out_files:
                tmp = pd.read_csv(i, delim_whitespace=True, header=None)
                df_list.append(tmp)
            dfFinal = pd.concat(df_list)
            #dfFinal.drop([0], inplace=True) 
            dfFinal.columns = clms
            dfFinal.set_index(clms[0], inplace=True)
            dfFinal.index = (dfFinal.index)*0.001
            DF_list.append(dfFinal[selected_cols])
        num_subplots = len(DF_list)
        num_rows = self.multiPlotDIM[num_subplots][0]
        num_cols = self.multiPlotDIM[num_subplots][1]
        print(num_subplots,num_rows,num_cols)
        print("Columns",clms,len(DF_list))
        print(DF_list[0].head())
        count = 0
        f, axarr = plt.subplots(num_rows*3, num_cols, sharex=True)
        #f.add_gridspec(nrows=3, ncols=3, left=0.15, right=0.48, wspace=0.05)
        # 78.42 KE -> Temp factor
        min_X, max_X, min_Y, max_Y = self.get_min_max_df_list_by_cols(DF_list)
        #print(min_X,max_X,min_Y,max_Y)
        subCounts = 0
        for i in range(0,num_rows*3,3):
            for j in range(num_cols):
                #print(i,j,subplots[subCounts])
                if i != 8:
                    plt.setp([a.get_xticklabels() for a in axarr[i,:]], visible=False)
                else:
                    axarr[i, j].set_xlabel('ns')
                if j != 0:
                    plt.setp([a.get_yticklabels() for a in axarr[:,j]], visible=False)
                else:
                    axarr[i, j].set_ylabel('KJ/mol/kpal')
                if count < num_subplots:
                    axarr[i, j].plot(DF_list[count],linewidth=0.8)
                    axarr[i, j].set_title(subplots[subCounts],fontsize=10)
                    subCounts += 1
                    axarr[i, j].set_xlim([min_X, max_X])
                    axarr[i, j].set_ylim([max_Y-12, max_Y])
                    axarr[i, j].spines['bottom'].set_visible(False)
                    axarr[i, j].tick_params(bottom=False, labelbottom=False)
                    d = .015  # how big to make the diagonal lines in axes coordinates
                    # arguments to pass to plot, just so we don't keep repeating them
                    kwargs = dict(transform=axarr[i, j].transAxes, color='k', clip_on=False)
                    axarr[i, j].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
                    axarr[i, j].plot((1 - d, 1 + d), (-d, +d), **kwargs)
                    ###########################################################
                    axarr[i+1, j].plot(DF_list[count],linewidth=0.8)
                    axarr[i+1, j].set_ylim([max_Y-17, max_Y-13])
                    axarr[i+1, j].spines['bottom'].set_visible(False)
                    axarr[i+1, j].spines['top'].set_visible(False)
                    axarr[i+1, j].tick_params(top=False,bottom=False, labelbottom=False)
                    d = .015
                    kwargs = dict(transform=axarr[i+1, j].transAxes, color='k', clip_on=False)
                    axarr[i+1, j].plot((-d, +d), (-d, +d), **kwargs)
                    axarr[i+1, j].plot((1 - d, 1 + d), (-d, +d), **kwargs)                    
                    kwargs.update(transform=axarr[i+1, j].transAxes)  # switch to the bottom axes
                    axarr[i+1, j].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
                    axarr[i+1, j].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
                    ###########################################################
                    axarr[i+2, j].plot(DF_list[count],linewidth=0.8)
                    axarr[i+2, j].set_ylim([min_Y, min_Y+8])
                    axarr[i+2, j].spines['top'].set_visible(False)
                    kwargs.update(transform=axarr[i+2, j].transAxes)  # switch to the bottom axes
                    axarr[i+2, j].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
                    axarr[i+2, j].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
                    count += 1
        # TODO Fix bottom lable 
        for i in range(num_cols):
            axarr[num_rows, i].set_xlabel('ns')
        #axarr[8, 0].set_xlabel('ns')
        #axarr[8, 1].set_xlabel('ns')
        #axarr[8, 2].set_xlabel('ns')
        #axarr[8, 3].set_xlabel('ns')
        plt.suptitle(self.plot_title)
        f.set_size_inches(18.5, 10.5)
        plt.legend(['PE(vdw)','PE(ee)','Pressure'],loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.subplots_adjust(left=None, bottom=0.1, right=None, top=None, wspace=None, hspace=0.5)
        plt.savefig(self.plot_out, dpi=f.dpi)
        plt.show()
                
    def get_min_max_df_list_by_cols(self,DF_list):
        min_X = 100000000
        max_X = -100000000
        min_Y = 100000000
        max_Y = -100000000
        for i in DF_list:
            if min(i.index) < min_X:
                min_X = min(i.index)
            if max(i.index) > max_X:
                max_X = max(i.index)
            if i.to_numpy().min() < min_Y:
                min_Y = i.to_numpy().min()+0.1
            if i.to_numpy().min() > max_Y:
                max_Y = i.to_numpy().max()+0.2
        return min_X, max_X, min_Y, max_Y
        
    def OH_LJ_Average(self):
        num_rows = self.multiPlotDIM[self.DF.shape[1]/4][0]
        num_cols = self.multiPlotDIM[self.DF.shape[1]/4][1]
        count = 0
        f, axarr = plt.subplots(num_rows, num_cols)
        min_X = min(self.DF.index)
        max_X = max(self.DF.index)
        min_Y = 0
        max_Y = 0
        for i in list(self.DF.columns):
            new_Y = np.float(i.split('_')[2][:-1])
            if new_Y > max_Y:
                max_Y = new_Y
        print("minXY",min_X,max_X,min_Y,max_Y)
        print("DF.index: ",self.DF.index)
        for i in range(num_rows):
            for j in range(num_cols):
                if i != (num_rows-1):
                    plt.setp([a.get_xticklabels() for a in axarr[i,:]], visible=False)
                else:
                    axarr[i, j].set_xlabel('ns')
                if j != 0:
                    plt.setp([a.get_yticklabels() for a in axarr[:,j]], visible=False)
                else:
                    axarr[i, j].set_ylabel('nm')
                if count < (len(self.DF.columns)):
                    subplot_title = self.DF.columns[count].split("_")
                    coefOH = np.polyfit(self.DF.index,self.DF[self.DF.columns[count]],1)
                    poly1d_fnOH = np.poly1d(coefOH)
                    axarr[i, j].text(0.2,max_Y-0.25,"y="+'{:.4f}'.format(coefOH[0])+"x+"+'{:.4f}'.format(coefOH[1]),c="red")
                    axarr[i, j].plot(self.DF.index,poly1d_fnOH(self.DF.index),
                                     linewidth=0.2,label=subplot_title[3],c="red")
                    axarr[i, j].plot(self.DF[self.DF.columns[count]],
                                     linewidth=0.5,label=subplot_title[3],c="red")
                    count += 1
                    coefOH = np.polyfit(self.DF.index,self.DF[self.DF.columns[count]],1)
                    poly1d_fnOH = np.poly1d(coefOH)
                    axarr[i, j].text(0.2,0.25,"y="+'{:.4f}'.format(coefOH[0])+"x+"+'{:.4f}'.format(coefOH[1]),c="green")
                    axarr[i, j].plot(self.DF.index,poly1d_fnOH(self.DF.index),
                                     linewidth=0.2,label=subplot_title[3],c="green")
                    subplot_title = self.DF.columns[count].split("_")
                    axarr[i, j].plot(self.DF[self.DF.columns[count]],
                                     linewidth=0.5,label=subplot_title[3],c="green")
                    #x = range(len(LJ_mZ))
                    #coefOH = np.polyfit(x,OH_mZ,1)
                    #poly1d_fnOH = np.poly1d(coefOH)
                    #coefLJ = np.polyfit(x,LJ_mZ,1)
                    #poly1d_fnLJ = np.poly1d(coefLJ)
                    count += 1
                    count += 1
                    axarr[i, j].set_title("BOX"+str(subplot_title[0:3])+" "+str(subplot_title[-1]),fontsize=10)
                    axarr[i, j].set_xlim([min_X, max_X])
                    axarr[i, j].set_ylim([min_Y, max_Y])
                    count += 1
        plt.suptitle(self.plot_title,fontsize=20)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        f.set_size_inches(18.5, 10.5)
        plt.savefig(self.plot_out, dpi=f.dpi)
        plt.show()
        # To add a mean line for each line
        #axarr[1, 2].axhline(mean_mean, color='blue')
        #axarr[1, 2].axhline(y=df[simulations[5]]['mean'], color='red')
        # If I want to draw a line fir to the plots do this below.

        #plt.title("Water and Argon, Zcoord Mean distribution v. Time")
        #plt.plot(x,OH_mZ,'b',x,poly1d_fnOH(x),'b--',LJ_mZ,'g',x,poly1d_fnLJ(x),'g--')
        #plt.show()        
    def OH_LJ_Zhistog(self):
        num_rows = self.multiPlotDIM[self.DF.shape[1]/4][0]
        num_cols = self.multiPlotDIM[self.DF.shape[1]/4][1]
        count = 0
        f, axarr = plt.subplots(num_rows, num_cols)
        min_X = min(self.DF.index)
        max_X = max(self.DF.index)
        min_Y = 0
        max_Y = 0
        for i in list(self.DF.columns):
            new_Y = np.float('2.0X_2.0Y_2.4Z_LJ_F'.split('_')[2][:-1])
        if new_Y > max_Y:
            max_Y = new_Y
        for i in range(num_rows):
            for j in range(num_cols):
                if i != (num_rows-1):
                    plt.setp([a.get_xticklabels() for a in axarr[i,:]], visible=False)
                else:
                    axarr[i, j].set_xlabel('ns')
                if j != 0:
                    plt.setp([a.get_yticklabels() for a in axarr[:,j]], visible=False)
                else:
                    axarr[i, j].set_ylabel('nm')
                if count < (len(self.DF.columns)):
                    count += 1
                    count += 1
                    subplot_title = self.DF.columns[count].split("_")
                    for nn in self.DF.index:
                        tmp = [float(n) for n in self.DF[self.DF.columns[count]].loc[nn].strip('][').split(', ')]
                        axarr[i, j].plot([nn]*len(tmp),tmp,linewidth=0.3,marker='x',c="red")
                    count += 1
                    subplot_title = self.DF.columns[count].split("_")
                    tmp = [float(n) for n in self.DF[self.DF.columns[count]].loc[0.02].strip('][').split(', ')]
                    axarr[i, j].plot([0]*len(tmp),tmp,
                                     linewidth=0.3,marker='+',c="green")                    
                    axarr[i, j].set_title("BOX"+str(subplot_title[0:3])+" "+str(subplot_title[-1]),fontsize=10)
                    axarr[i, j].set_xlim([min_X, max_X])
                    axarr[i, j].set_ylim([min_Y, max_Y])
                    count += 1
        #plt.suptitle(self.plot_title,fontsize=20)
        #plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        #f.set_size_inches(18.5, 10.5)
        #plt.savefig(self.plot_out, dpi=f.dpi)
        plt.show()        
    def check_file_names(self,individual_files):
        # This same class also exist in structure/XYZ_Formats.
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
                    A naming conventions s        # To add a mean line for each line
        #axarr[1, 2].axhline(mean_mean, color='blue')
        #axarr[1, 2].axhline(y=df[simulations[5]]['mean'], color='red')
        # If I want to draw a line fir to the plots do this below.
        #x = range(len(LJ_mZ))
        #coefOH = np.polyfit(x,OH_mZ,1)
        #poly1d_fnOH = np.poly1d(coefOH)
        #coefLJ = np.polyfit(x,LJ_mZ,1)
        #poly1d_fnLJ = np.poly1d(coefLJ)
        #plt.title("Water and Argon, Zcoord Mean distribution v. Time")
        #plt.plot(x,OH_mZ,'b',x,poly1d_fnOH(x),'b--',LJ_mZ,'g',x,poly1d_fnLJ(x),'g--')
        #plt.show()uch as A_1.dcd, A_2.dcd ... will 
                    represent two DCDs belonging to the same contiguous 
                    simulation. If this lables are different, the first one
                    will be used to identified the concatenated trajectories.""")
        return prefix_label
        
        