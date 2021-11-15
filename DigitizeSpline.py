#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 13:53:23 2021

@author: arun
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy as sp
from scipy import interpolate as ip
from itertools import islice
import ars_support_functions as sf
# def range_list(k,n=0):
#     if n==0:
#         n=k
#     a= list(range(k+1))
#     b=list(item/k*n for item in a)
#     return (b)

class DigitizeSpline():
    """ Converts an input file into a spline representation Example input file below (copy and paste and you can run it)
    # 2-D boundary description
    # "#" is the comment line signal
    # *SIZE= How many loops do we have in the boundary : go ccw it externl and CW for internal boundaries 
    #*LOOP is to list the segents in a loops Dont put comments in the loop section 
    #each segement has Number of points, Order of spline, points x,y,x,y,x,y etc 
    #*END is the end of the file 
    *SIZE
    1
    *LOOP
    10,1, 0,0, 1,0
    10,1, 1,0, 1,1
    10,1, 1,1, 0,1
    10,1, 0,1, 0,0
    *END
    """
    
    def __init__(self,fname ):
        self.filename=fname
        # this is where the commands in the input file are converted into functions in a safe way. Add new if you find 
        self.function_switch={}
        self.function_switch["*SIZE\n"]=self._read_size
        self.function_switch["*LOOP\n"]=self._read_loop
        self.read_input()
        plt.close('all')
        return
    
    def read_input(self):
        print("\n Hello, Welcome to SPLINECONTOURS, A way to generate CONTOURS (2-D REGION BOUNDARIES)  using SPLINE sections ")
        while True:
            try:
                self.fp=open(self.filename,'r')
            except IOError:
                print("\n Sorry I couldnt find this file, please try again\n")
                continue
            else:
                break
        
        self.loops=[]
        self.total_nodes=0
        
        sf.read_files(self.fp,self.function_switch)
        
        self.fp.close()
        return
    
    def _read_size(self):
        """ reads the number of loops in the file and find the next token or sends an end token"""
        print("read_size")
        self.NLOOPS= int(self.fp.readline())
        line=self.fp.readline()
        token= line if(line) else "*END"
        return token
        
    def _read_loop(self):
       token,segments=sf.read_number_list(self.fp,[int,int,int,float])
       nodes_in_loop=0
       for seg in segments:
           nodes_in_loop+=seg[0]-1
           spline_k=int(seg[1])
           periodic=int(seg[2])
           length=len(seg[3::2]) # first one is the number of points , second is the spline  
           t=list(range(length))#np.arange(0,length,1)
           # convert to spline representation
           tckx = ip.splrep(t, seg[3::2], k=spline_k,s=0,per=periodic)
           tcky = ip.splrep(t, seg[4::2], k=spline_k,s=0,per=periodic)
           seg.append(tckx)
           seg.append(tcky)
           
       segments.append(nodes_in_loop)    
       self.total_nodes+=nodes_in_loop
       self.loops.append(segments)
       return token
      
        
    def distribute_nodes(self,Ncycles):
        number_of_loops = self.NLOOPS
        alpha=0.25 # this is the force power 
        eta=0.01
        loc2=self.filename.rfind(".")
        loc1=self.filename.rfind("_")
        if(loc1==-1):
            loc1=loc2
        print("Assigning Nodes....\n"+"...opening output file...\n")
        outfile=self.filename[0:loc1]+"_bdrypts"+self.filename[loc2:]
        if os.path.exists(outfile):
            os.remove(outfile) 
        #___________________#________
        with open(outfile,'a+') as f:
            f.write("# number of loops and nodes per loop\n")
            f.write('# dont forget to include injection points and refniment regions\n')
            f.write(str(self.NLOOPS))
            for loop in self.loops:
                f.write(' '+str(loop[-1]))
            f.write('\n')
            f.write('# boundary node locations\n')
            kk=1 # this is the node type
            for loop in self.loops:
                for seg in loop[:-1]:
                    n_nodes=seg[0]
                    length=len(seg[3:-2:2])
                    t=np.linspace(0,length-1,n_nodes)
                    tckx=seg[-2]
                    tcky=seg[-1]
                    Forcev=np.zeros(n_nodes)
                    #____________#________adjust the force s
                    for j in range(Ncycles):
                        Xv=ip.splev(t, tckx, der=0)
                        Yv=ip.splev(t, tcky, der=0)
                        Forcev[:]=0
                        for i in range(1,n_nodes-1):
                            # distance based attractive forces
                            Forcev[i]=-((Xv[i-1]-Xv[i])**2+(Yv[i-1]-Yv[i])**2)**alpha+((Yv[i+1]-Yv[i])**2+(Xv[i+1]-Xv[i])**2)**alpha
                        Forcev[0]=0
                        Forcev[-1]=0
                        # force on the last two nodes are zero 
                        t+=eta*Forcev # adjust the t values to take care of the forces 
                    plt.plot(Xv,Yv,'.')
                    plt.axis('equal')
                    for i in range(0,n_nodes-1):
                          f.write(str(Xv[i])+' '+str(Yv[i])+' '+str(kk)+'\n')
                    kk+=1 # next set
        f.closed
        print("Output written in "+outfile)
        return
    
    def set_injection_points(self,injvals,refvals=[]):
        loc2=self.filename.rfind(".")
        loc1=self.filename.rfind("_")
        outfile=self.filename[0:loc1]+"_bdrypts"+self.filename[loc2:]
        number=len(injvals)
        with open(outfile,'a+') as f:
            f.write("# number of injection loctions\n")
            f.write(str(number)+"\n")
            f.write('# injection point locations\n')
            for pt in injvals:
                f.write(str(pt[0])+' '+str(pt[1])+'\n')
                plt.plot(pt[0],pt[1],'+',color='r')
            f.write('# number of refinement regions\n')
            f.write(str(0))
        f.closed
        print("Injection points  added to "+outfile)
        return
                    
    def set_refinement_regions(self,refvals=[]):
        loc2=self.filename.rfind(".")
        loc1=self.filename.rfind("_")
        outfile=self.filename[0:loc1]+"_bdrypts"+self.filename[loc2:]
        with open(outfile,'a+') as f:
            ref_n=len(refvals)
            f.write("# number of refinement areas\n")
            f.write(str(ref_n)+"\n")
            if(ref_n >0):
                f.write('# refinement region corner locations\n')
                for pts in refvals:
                    for val in pts[:-1]:
                        f.write(str(val)+' ')
                        f.write(str(pts[-1]))
                    f.write('\n')
             
        f.closed
        print("refinement areas added to "+outfile)
        return             
    
                 