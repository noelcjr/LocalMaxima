import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
import mpl_toolkits.axes_grid1
import matplotlib.widgets
import pandas as pd
import argparse
import sys
import os
#import LocalMaxima.structure.XYZ_Formats as xyzf

class Player(FuncAnimation):
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
                 save_count=None, mini=0, maxi=100, pos=(0.125, 0.92), **kwargs):
        self.i = 0
        self.min=mini
        self.max=maxi
        self.runs = True
        self.forwards = True
        self.fig = fig
        self.func = func
        self.setup(pos)
        FuncAnimation.__init__(self,self.fig, self.func, frames=self.play(), 
                                           init_func=init_func, fargs=fargs,
                                           save_count=save_count, **kwargs )    
    def play(self):
        while self.runs:
            self.i = self.i+self.forwards-(not self.forwards)
            if self.i > self.min and self.i < self.max:
                yield self.i
            else:
                self.stop()
                yield self.i
    def start(self):
        self.runs=True
        self.event_source.start()
    def stop(self, event=None):
        self.runs = False
        self.event_source.stop()
    def forward(self, event=None):
        self.forwards = True
        self.start()
    def backward(self, event=None):
        self.forwards = False
        self.start()
    def oneforward(self, event=None):
        self.forwards = True
        self.onestep()
    def onebackward(self, event=None):
        self.forwards = False
        self.onestep()
    def onestep(self):
        if self.i > self.min and self.i < self.max:
            self.i = self.i+self.forwards-(not self.forwards)
        elif self.i == self.min and self.forwards:
            self.i+=1
        elif self.i == self.max and not self.forwards:
            self.i-=1
        self.func(self.i)
        self.fig.canvas.draw_idle()
    def setup(self, pos):
        playerax = self.fig.add_axes([pos[0],pos[1], 0.22, 0.04])
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(playerax)
        bax = divider.append_axes("right", size="80%", pad=0.05)
        sax = divider.append_axes("right", size="80%", pad=0.05)
        fax = divider.append_axes("right", size="80%", pad=0.05)
        ofax = divider.append_axes("right", size="100%", pad=0.05)
        self.button_oneback = matplotlib.widgets.Button(playerax, label=u'$\u29CF$')
        self.button_back = matplotlib.widgets.Button(bax, label=u'$\u25C0$')
        self.button_stop = matplotlib.widgets.Button(sax, label=u'$\u25A0$')
        self.button_forward = matplotlib.widgets.Button(fax, label=u'$\u25B6$')
        self.button_oneforward = matplotlib.widgets.Button(ofax, label=u'$\u29D0$')
        self.button_oneback.on_clicked(self.onebackward)
        self.button_back.on_clicked(self.backward)
        self.button_stop.on_clicked(self.stop)
        self.button_forward.on_clicked(self.forward)
        self.button_oneforward.on_clicked(self.oneforward)

def register_parser(subparsers):
    parser = subparsers.add_parser('zhist', usage=usage(), \
                                   description=description())
    add_arguments(parser)

def add_arguments(parser):
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument("-c","--csv", metavar="FILE",
                              help="A comma separated file with lists of z-coordinate \
                                    distributions of selected atoms.", required=True)
    requiredArgs.add_argument("-o","--outmp4", metavar="FILE",
                              help="File name for output video in .mp4 format.",
                              required=True)
    optionalArgs = parser.add_argument_group("Optional arguments")
    optionalArgs.add_argument("-s","--single", type=str, 
                               help="A String with only two columns names for \
                                   more detailed analysis of a moving histogram.\
                                   If present, only one moving histogram will be\
                                   output and displayed in a window.",
                               required=False)
    parser.set_defaults(func=run)
### using this class is as easy as using FuncAnimation:            
# https://stackoverflow.com/questions/44985966/managing-dynamic-plotting-in-matplotlib-animation-module
#dd='/home/noel/Datasets_nobackup/Phase_Diag/t3_216_wt0_108_150EE_T2/288/t3_216_wt0_108_150EE_T2_288_zmean_1.csv'
#dd='/home/noel/Datasets_nobackup/Phase_Diag/t5_216_wt0_108_150EE_T2/298/t5_216_wt0_108_150EE_T2_298_zmean_1.csv'
#DF = pd.read_csv(dd,index_col=0)
#optionssingle = '2.0X_2.0Y_2.0Z_OH_A_L,2.0X_2.0Y_2.0Z_LJ_A_L'
def run(options):
    if not os.path.exists(options.csv):
        print("ERROR: CSV input file does not exist.")
        sys.exit(1)
    DF = pd.read_csv(options.csv,index_col=0)
    n = DF.shape[0]
    print(DF.shape)
    data = []
    dataT = []
    maxZ = 0
    if options.single:
        col_list = []
        min_max = {'minX':0,'maxX':40,'minY':0,'maxY':-1000000}
        for ii in options.single.split(','):
            array = []
            for nn in DF.index: 
                array.append([float(n) for n in DF[ii].loc[nn].strip('][').split(', ')])
            data.append(np.array(array[:n]))
            dataT.append(ii)
            if maxZ < float(ii.split('_')[2][:-1])*10:
                maxZ = float(ii.split('_')[2][:-1])*10
        min_max['maxY'] = maxZ/10
        bns = [float(i/10) for i in range(0,int(maxZ)+2)]
        fig, ax = plt.subplots(1)
        def update(i):
            ax.cla()
            ax.set_xlim([min_max['minX'],min_max['maxX']])
            ax.set_ylim([min_max['minY'],min_max['maxY']+0.1])
            ax.hist([data[0][i],data[1][i]],bins=bns,orientation='horizontal')
            subtitle = "".join(dataT[0].split('_')[0:3])+" "+dataT[0].split('_')[4]
            ax.set_title(subtitle,fontsize=10)
            ax.text(35,2.2,"n="+str(i))
            ax.set_xlabel('Count')
            ax.set_ylabel('Z(nm)')
        plt.suptitle("Water (Blue) and Argon (Orange), Zcoord distribution")
        ani = Player(fig, update, maxi=n-1, save_count=n)
        plt.gca().invert_yaxis()
        fig.set_size_inches(18.5, 10.5)
        Writer = writers['ffmpeg']
        writer = Writer(fps=5, bitrate=1800)
        plt.show()
        #if os.path.basename(options.outmp4).split('.')[1] == 'mp4':
        #    ani.save(options.outmp4, writer)
        #else:
        #    outfile = options.outmp4+".mp4"
        #    print("WARNING: .mp4 suffix added to output file.")
    else:
        print("Else")
        col_list = []
        min_max = {'minX':0,'maxX':40,'minY':0,'maxY':-1000000}
        for i in DF.columns:
            if i.split('_')[-1] == 'L':
                tmp = float(i.split('_')[2][:-1])
                if tmp > min_max['maxY']:
                    min_max['maxY'] = tmp
                col_list.append(i)
        multiPlotDIM = {1:(1,1),2:(1,2),3:(2,2),4:(2,2),5:(2,3),6:(2,3),
                        7:(3,3),8:(3,3),9:(3,3),10:(3,4),11:(3,4),12:(3,4),
                        30:(5,6),36:(6,6)} 
        num_rows = multiPlotDIM[DF.shape[1]/4][0]
        num_cols = multiPlotDIM[DF.shape[1]/4][1]
        print(num_rows,num_cols)
        for ii in range(2,len(DF.columns),4):
            array = []
            isnan = pd.isna(DF[DF.columns[ii]])
            for nn in DF.index:
                if isnan.loc[nn]:
                    array.append([])
                else:
                    array.append([float(n) for n in DF[DF.columns[ii]].loc[nn].strip('][').split(', ')])
            data.append(np.array(array[:n]))
            dataT.append(DF.columns[ii])
            if maxZ < float(DF.columns[ii].split('_')[2][:-1])*10:
                maxZ = float(DF.columns[ii].split('_')[2][:-1])*10
            array = []
            for nn in DF.index:
                if isnan.loc[nn]:
                    array.append([])
                else:
                    array.append([float(n) for n in DF[DF.columns[ii+1]].loc[nn].strip('][').split(', ')])
            data.append(np.array(array[:n]))
            dataT.append(DF.columns[ii+1])
            if maxZ < float(DF.columns[ii+1].split('_')[2][:-1])*10:
                maxZ = float(DF.columns[ii+1].split('_')[2][:-1])*10
        print(data)
        #print(dataT)
        bns = [float(i/10) for i in range(0,int(maxZ)+2)]
        print(f"bns {bns}")
        fig, ax = plt.subplots(num_rows,num_cols)
        def update(i):
            count = 0
            for ii in range(num_rows):
                for j in range(num_cols):
                    print(f"num_rows = {ii}  num_col = {j}")
                    if count < (len(DF.columns)):
                        ax[ii,j].cla()
                        ax[ii,j].set_xlim([min_max['minX'],min_max['maxX']])
                        ax[ii,j].set_ylim([min_max['minY'],min_max['maxY']+0.1])
                        ax[ii,j].hist([data[count][i],data[count+1][i]],bins=bns,orientation='horizontal')
                        subtitle = "".join(dataT[count].split('_')[0:3])+" "+dataT[count].split('_')[4]
                        ax[ii,j].set_title(subtitle,fontsize=10)
                        ax[0,0].text((min_max['maxX']*(num_cols-1)),((maxZ/10)+0.5),"n="+str(i))
                        count += 2
                    if ii != (num_rows-1):
                        plt.setp([a.get_xticklabels() for a in ax[ii,:]], visible=False)
                    else:
                        ax[ii, j].set_xlabel('Count')
                    if j != 0:
                        plt.setp([a.get_yticklabels() for a in ax[:,j]], visible=False)
                    else:
                        ax[ii, j].set_ylabel('Z(nm)')
        plt.suptitle("Water (Blue) and Argon (Orange), Zcoord distribution")
        ani = Player(fig, update, maxi=n-1, save_count=n)
        plt.gca().invert_yaxis()
        fig.set_size_inches(18.5, 10.5)
        Writer = writers['ffmpeg']
        writer = Writer(fps=5, bitrate=1800)
        if os.path.basename(options.outmp4).split('.')[1] == 'mp4':
            ani.save(options.outmp4, writer)
        else:
            outfile = options.outmp4+".mp4"
            print("WARNING: .mp4 suffix added to output file.")
    #plt.show()
def description():
    return '''Gets the mean z-axis location and a list of z-axis distributions for
              selected atoms of one or more DCD files. The output of this function
              canb be used to build moving z-distributions histograms with, for example:
              ntral zhist --csv A_0.csv --outmp4 A_0.mp4
           '''
def usage():
    return '''\nntraj zhist --csv out.csv --outmp4 out.mp4'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)

