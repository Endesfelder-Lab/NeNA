# -*- coding: cp1252 -*-
#Coordinate Based Precision Estimator
#please contact
#<Sebastian Malkusch> <malkusch@chemie.uni-frankfurt.de>
### BEGIN LICENSE
# Copyright (C) 2013 <Sebastian Malkusch> <malkusch@chemie.uni-frankfurt.de>
# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License version 3, as published 
# by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranties of 
# MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR 
# PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with this program.  If not, see <http://www.gnu.org/licenses/>.
### END LICENSE

from Tkinter import *
import tkFileDialog
from numpy import*
from scipy import*
from scipy.optimize import curve_fit
import numpy as np
import scipy as sp
import gc


def print_menu():
    print('1. define roi')
    print('2. print roi')
    print('3. coordinate based precision estimation')
    print('4. quit')
    print ' '

def def_roi(roi_frame):
    roi_frame=zeros([5,2])
    roi_frame[0,0]=input('Xmin [um]: ')*1000
    roi_frame[0,1]=input('Xmax [um]: ')*1000
    roi_frame[1,0]=input('Ymin [um]: ')*1000
    roi_frame[1,1]=input('Ymax [um]: ')*1000
    roi_frame[2,0]=input('Zmin [um]: ')*1000
    roi_frame[2,1]=input('Zmax [um]: ')*1000
    roi_frame[3,0]=input('Tmin [frame]: ')
    roi_frame[3,1]=input('Tmax [frame]: ')
    roi_frame[4,0]=input('Imin [counts]: ')
    roi_frame[4,1]=input('Imax [counts]: ')
    return roi_frame
    
def print_roi(roi_frame):
    print('Xmin [nm]:     '+ str(int(roi_frame[0,0])))
    print('Xmax [nm]:     '+ str(int(roi_frame[0,1])))
    print('Ymin [nm]:     '+ str(int(roi_frame[1,0])))
    print('Ymax [nm]:     '+ str(int(roi_frame[1,1])))
    print('Zmin [nm]:     '+ str(int(roi_frame[2,0])))
    print('Zmax [nm]:     '+ str(int(roi_frame[2,1])))
    print('Tmin [frame]:  '+ str(int(roi_frame[3,0])))
    print('Tmax [frame]:  '+ str(int(roi_frame[3,1])))
    print('Imin [counts]: '+ str(int(roi_frame[4,0])))
    print('Imax [counts]: '+ str(int(roi_frame[4,1])))

def filename(Name):
    master = Tk()
    Name = tkFileDialog.askopenfilename(title="Open File (rapidSTORM Malk format)", filetypes=[('MALK File','*.txt')])
    master.quit()
    master.destroy()
    return Name

def loadfile(Locs, Name):
    locs = np.loadtxt(Name, skiprows = 1)
    l=locs.shape
    Mlocs=zeros([l[0],5])
    if l[1]==4:
        Mlocs[:,0]=locs[:,0]
        Mlocs[:,1]=locs[:,1]
        Mlocs[:,3]=locs[:,2]
        Mlocs[:,4]=locs[:,3]
    else:
        Mlocs=locs
    return Mlocs

def ROI(locs,roi_frame):
    idx1=locs[:,0]>=roi_frame[0,0]
    idx2=locs[:,0]<=roi_frame[0,1]
    idy1=locs[:,1]>=roi_frame[1,0]
    idy2=locs[:,1]<=roi_frame[1,1]
    idz1=locs[:,2]>=roi_frame[2,0]
    idz2=locs[:,2]<=roi_frame[2,1]
    idt1=locs[:,3]>=roi_frame[3,0]
    idt2=locs[:,3]<=roi_frame[3,1]
    idi1=locs[:,4]>=roi_frame[4,0]
    idi2=locs[:,4]<=roi_frame[4,1]
    locs1 = (locs[(idx1&idx2&idy1&idy2&idz1&idz2&idt1&idt2&idi1&idi2), :])
    return locs1

def save_locs(roi_locs,Name):
    outfilename=Name[0:(len(Name)-4)]+'-ROI'+'.txt'
    out_file = open(outfilename, "w")
    out_file.write('# localizations whitin roi'+"\n")
    out_file.close
    out_file = open(outfilename, "a")
    for i in range(1, (len(roi_locs))):
        out_file.write(str(roi_locs[i,0])+' '+str(roi_locs[i,1])+' '+str(roi_locs[i,2])+' '+str(int(roi_locs[i,3]))+' '+str(roi_locs[i,4])+"\n")
    out_file.close
    print "done!"

def save_frame(roi_frame,Name):
    outfilename=Name[0:(len(Name)-4)]+'-frame'+'.txt'
    out_file = open(outfilename, "w")
    out_file.write('# roi dimensions'+"\n")
    out_file.close
    out_file = open(outfilename, "a")
    out_file.write("Xmin [nm]   " + str(int(roi_frame[0,0]))+"\n")
    out_file.write("Xmax [nm]   " + str(int(roi_frame[0,1]))+"\n")
    out_file.write("Ymin [nm]   " + str(int(roi_frame[1,0]))+"\n")
    out_file.write("Ymax [nm]   " + str(int(roi_frame[1,1]))+"\n")
    out_file.write("Zmin [nm]   " + str(int(roi_frame[2,0]))+"\n")
    out_file.write("Zmax [nm]   " + str(int(roi_frame[2,1]))+"\n")
    out_file.write("Tmin [nm]   " + str(int(roi_frame[3,0]))+"\n")
    out_file.write("Tmax [nm]   " + str(int(roi_frame[3,1]))+"\n")
    out_file.write("Imin [nm]   " + str(int(roi_frame[4,0]))+"\n")
    out_file.write("Imax [nm]   " + str(int(roi_frame[4,1])))    
    out_file.close
    print "done!"

def N_N(array,value):
    idx=(np.abs(array-value)).argmin()
    return array[idx]

def Acc_Calculator(Nlocs,MaxFrame):
    print 'estimating coordinate based localization precision...'
    #sort Nlocs
    leftMax=np.searchsorted(Nlocs[:,3],MaxFrame,side='left')
    lengthlocs=leftMax-1
    NearNeigh= zeros ([lengthlocs,1])    
    for i in range(0,lengthlocs):
        loc=Nlocs[i,:]
        j=1+loc[3]
        if i==0:
            left=np.searchsorted(Nlocs[:,3],j,side='left')
            right=np.searchsorted(Nlocs[:,3], j, side = 'right')
            Frlocs=Nlocs[left:right,:]
        elif loc[3]>Nlocs[(i-1),3]:
            left=np.searchsorted(Nlocs[:,3], j, side='left')
            right=np.searchsorted(Nlocs[:,3], j, side='right')
            Frlocs=Nlocs[left:right,:]
        if left==right:
           NearNeigh[i]=0
        else:
           if j>(max(Nlocs[:,3])):
              NearNeigh[i]=200
           else:
              Dist=sqrt(((loc[0]-Frlocs[:,0])*(loc[0]-Frlocs[:,0]))+((loc[1]-Frlocs[:,1])*(loc[1]-Frlocs[:,1]))+((loc[2]-Frlocs[:,2])*(loc[2]-Frlocs[:,2])))
              NearNeigh[i]=N_N(Dist,0)
              if NearNeigh[i]>200:
                 NearNeigh[i]=200
                    
    print 'Done!'
    return NearNeigh    

def Area(r,y):
    Areaf=abs(np.trapz(y, r))
    return Areaf

def CFunc2dCorr(r,a,rc,w,F,A,O):
    y=(r/(2*a*a))*exp((-1)*r*r/(4*a*a))*A+(F/(w*sqrt(pi*2)))*exp(-0.5*((r-rc)/w)*((r-rc)/w))+O*r

    return y

def CFit_resultsCorr(r,y):
    A=Area(r,y)
    print A
    p0 = np.array([10.0,201,100,(A/2),(A/2),((y[98]/200))])
    popt, pcov = curve_fit(CFunc2dCorr,r,y,p0)
    print 'estimated localization precision at:'
    print popt
    return popt, pcov

def save_results(ar,ay,ayf,aF,aFerr, Name):
    outfilename=Name[0:(len(Name)-4)]+'-precision'+'.txt'
    out_file = open(outfilename, "w")
    out_file.write('# localization precision whitin roi at: ' +str(aF[0])+'+/-'+str(aFerr[0,0])+' [nm]'+"\n")
    out_file.write('# r[nm] hist[a.u.]  fit[a.u.]'+"\n")
    out_file.close
    out_file = open(outfilename, "a")
    for i in range(1, (len(ar))):
        out_file.write(str(int(ar[i]))+'    '+str(ay[i])+'  '+str(ayf[i])+"\n")
    out_file.close
    print "done!"



def Loc_Acc(roi_frame):
    Name=[]
    Locs=[]
    Name=filename(Name)
    save_frame(roi_frame,Name)
    locs=loadfile(Locs,Name)
    roi_locs=ROI(locs,roi_frame)
    save_locs(roi_locs,Name)
    NearNeigh=Acc_Calculator(roi_locs,roi_frame[3,1])
    ahist=np.histogram(NearNeigh,bins=99,range=(1,199),normed=False)
    ar=ahist[1][1:len(ahist[1])]-1
    ay=ahist[0]
    aF,aFerr=CFit_resultsCorr(ar,ay)
    ayf=CFunc2dCorr(ar,aF[0],aF[1],aF[2],aF[3],aF[4],aF[5])
    save_results(ar,ay,ayf,aF,aFerr,Name)
    check()

#Garbage Collector
def check():
    gc.collect()
    print 'Job done!'
    
# Main Program
roi_frame = []
while True:
    print
    print ('coordinate based localization precision estimator:')
    print 
    print_menu()
    menu_choice = int(input('type in a number (1-4): '))
    if menu_choice==1:
        print ' '
        print 'define roi:'
        print ' '
        roi_frame=def_roi(roi_frame)
    elif menu_choice == 2:
        if roi_frame == []:
            print ' '
            print 'define roi first!'
            print ' '
        else:
            print ' '
            print 'defined roi:'
            print ' '
            print_roi(roi_frame)
    elif menu_choice == 3:
        if roi_frame == []:
            print ' '
            print 'define roi first!'
            print ' '
        else:
            Loc_Acc(roi_frame)
            
    elif menu_choice == 4:
        break
    else:
        print_menu()
