#!python
#cython: language_level=3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:45:19 2021

@author: Arun
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import ctypes
from numpy.ctypeslib import ndpointer
import  ars_delaunay_new as ars_fns
from itertools import islice

def f_mylines(inOpenedFile, inN):
    "function to extract N lines at a time from a csv file"
    return [x.strip() for x in islice(inOpenedFile, inN)]

class FlowMesherPython():
    def __init__(self,fname,delimiter=" ",nodetype=-1):
        self.fname=fname
        loc2=fname.rfind(".")
        loc1=fname.rfind("_")
        if(loc1==-1):
            loc1=loc2
        self.filename=fname[:loc1]+"_data"+fname[loc2:]
       
        #THere are no checks right now for format errors, I am still learning how to do this correctly
# this is a 2-D program
        self.read_file(delimiter)
        self.initialize_geometry()
        self.nodetype=nodetype
        self.delim=delimiter
#_____________________________#________________________________________________ 
    def read_file(self,delim):
        """INPUT FILE 
        # number of loops and nodes per loop
        2,36,19
        # boundary node locations (36+19), last column is the boundary type, useful for assigning boundary conditions 
        -1.0,-1.0,1
        -0.7777777777777777,-1.0,1
        -0.5555555555555556,-1.0,1
        .
        .
        .
        #Number Injection points
        4
        0.5, 0.5
        -0.5,0.5
        0.5,-0.5
        -0.5,-0.5
        #Refinement regions polygon corners 
        3, -1, -1, 1,1 1,-1 
        """

        print("\n Hello, Welcome to FLOWMESHER, A SIMPLE 2_D MeshGenerator\n Given the bounday nodes\n")
        while True:
            try:
                #self.filename= input(" Please type-in the path to your input file and press 'Enter': ")
                inputFile=open(self.fname,'r')
            except IOError:
                print("\n Sorry I couldnt find this file, please try again\n")
                continue
            else:
                break
        data = np.loadtxt(inputFile, delimiter=delim, dtype=np.intc,max_rows=1)
        self.N_SEGMENT_NODES=np.cumsum(data[1:])
        self.N_BDRY_SEGMENTS=data[0]
        self.N_BNODES=sum(data[1:])

        self.BRm=np.loadtxt(inputFile, delimiter=delim, dtype=np.single,max_rows=self.N_BNODES,usecols=[0,1]) 
        #now find the injection locations 
        self.N_INJ=np.loadtxt(inputFile, delimiter=delim, dtype=np.intc,max_rows=1).item() 
        
        self.INJ_DATA=np.loadtxt(inputFile, delimiter=delim, dtype=np.single,max_rows=self.N_INJ) 
        if self.N_INJ==1:
            self.INJ_DATA=np.array([self.INJ_DATA.tolist()],dtype=np.single)  
        # now if you want to add some refinements here is the place to do it 
        self.N_REFINEMENTS=np.loadtxt(inputFile, delimiter=delim, dtype=np.intc,max_rows=1).item(0)
        self.REFINEMENT_DATA=[]
        if (self.N_REFINEMENTS>0):
            for i in range(self.N_REFINEMENTS):
                self.REFINEMENT_DATA.append(np.loadtxt(inputFile, delimiter=delim, dtype=np.single,max_rows=1))
        inputFile.close()
        self.constrained_edges=[]# needed for Boundary Condition Specification
        
    
        return
#_________________________
    def initialize_geometry(self):
        # create tangent vectors 
        self.TGTm=np.zeros_like(self.BRm)
        start=0 
        for i in self.N_SEGMENT_NODES:
            self.TGTm[start:i,:]=np.diff(self.BRm[start:i,:], axis=0,append=[self.BRm[start,:]])
            start=i
        refinedarea=self.__find_refined_area()
       
        # #Area= xdy-ydx
        AREA=(sum(self.BRm[:,0]*self.TGTm[:,1]-self.BRm[:,1]*self.TGTm[:,0]))/2
        # Convert to unit tgt vector 
        BDRY_DISTv=np.sqrt(self.TGTm[:,0]**2+self.TGTm[:,1]**2)
        self.TGTm/=BDRY_DISTv[:,None]
        self.tgt_v=np.reshape(self.TGTm,2*self.N_BNODES)
        # max distance in the bdry nodes 
        self.R0= np.max(BDRY_DISTv)
        
        # Interior nodes:
        self.N_INODES=int((AREA+refinedarea)/(np.pi*(self.R0/2)**2))+20



        # Now assign the nodes and their properties 
        self.N_TOTAL=self.N_INODES+len(self.BRm) 
        print("interior nodes= ", self.N_INODES, "boundary nodes= ", self.N_BNODES," total nodes= ", self.N_TOTAL)
        ## better to store the nodes as x, y in a single array to save a lot of pain
        #_slices of data from storage________________________________________________
        self.data=np.zeros((6*self.N_TOTAL),dtype=np.single) # this contains all stored data in a contiguous list 
        self.data6=self.data.reshape(6,self.N_TOTAL)
        self.xs,self.ys,self.vxs, self.vys, self.fxs, self.fys=self.data6
        self.rs=self.data6[:2,:]
        self.xs[:self.N_BNODES],self.ys[:self.N_BNODES]=self.BRm[:,0],self.BRm[:,1]
        self.NODE_TYPES=np.zeros(self.N_TOTAL,dtype=np.intc)
        # assign the injected nodes 
        
        r_0= np.random.rand(self.N_INODES)*self.R0/2
        theta=np.random.rand(self.N_INODES)*np.pi*2
        cs=np.cos(theta)
        sine=np.sin(theta)
        j_start=self.N_BNODES
        self.xs[j_start:],  self.ys[j_start:] =r_0*cs, r_0*sine
        self.vxs[j_start:],  self.vys[j_start:] =0.5*cs, 0.5*sine
     
        j=0

        for i in range(self.N_INODES):
            j=j%self.N_INJ
            self.xs[i+j_start]+=self.INJ_DATA[j,0]
            self.ys[i+j_start]+=self.INJ_DATA[j,1]
            j=j+1
        # #_________________________________    
        self.MAX=np.max(self.BRm,axis=0)+2*self.R0
        self.MIN=np.min(self.BRm,axis=0)-2*self.R0
        self.MIN_X,self.MAX_X=self.MIN[0],self.MAX[0]
        self.MIN_Y,self.MAX_Y=self.MIN[1],self.MAX[1]
        return
        #__________________________________________        
    
    def showgrid(self):
        plt.close('all')
        self.fig=plt.figure()
        
        ax = plt.axes(xlim=(self.MIN_X, self.MAX_X), ylim=(self.MIN_Y, self.MAX_Y))
        plt.axis('equal')
        self.line, = ax.plot([], [], linestyle='none', marker='.', color='k')
        self.line.set_data(self.xs[:self.N_BNODES],self.ys[:self.N_BNODES])
        plt.show()
#
    
    def triangulate(self):
        
        self.tri=ars_fns.f_ars_triangulate(self.data6[:2,:],self.NODE_TYPES,self.TGTm,self.filename)
        # this writes output to the file also
# #____________________________#________________________#_________________#

    def __find_refined_area(self):
        refinedarea=0
        if (self.N_REFINEMENTS>0):
            for refi in self.REFINEMENT_DATA:
                N=len(refi)-2 # first two contatin ht and sharpness data 
                ht=refi[0]
                X=refi[2:].reshape(int(N/2),2)
                tgt=np.diff(X,axis=0)
                refinedarea+=((sum(X[:-1,0]*tgt[:,1]-X[:-1,1]*tgt[:,0]))/2)*(1+ht)**2
        return refinedarea
# #cdef functions for doing the cals go here 

#___________________________________________________________________________________________
    def _move_node_c(self, MAX_ITER):
        if (self.nodetype==-1):
            inputFile=open(self.fname,'r')
            self.NODE_TYPES[:self.N_BNODES]= np.loadtxt(inputFile, delimiter=self.delim, dtype=np.intc,max_rows=self.N_BNODES,usecols=[2],skiprows=3) 
            inputFile.close()
        else:
            self.NODE_TYPES[:self.N_BNODES]=self.nodetype 
        
        ref_lengths=np.array([0]+[len(x) for x in self.REFINEMENT_DATA],dtype=np.intc)
        ref_v=np.array([x for array in self.REFINEMENT_DATA for x in array],dtype=np.single)
        print(ref_lengths,ref_v)
        mydir=os.getcwd()+"/move_nodes_ctypes.so"
        #"/media/arun/ARUN_PROG/PROGRAMMING/MeshGenerator/move_nodes_ctypes.so"
        my_functions=ctypes.cdll.LoadLibrary(mydir)
        fun =my_functions.f_moveNodes
        fun.restype = ctypes.c_double
        
        self.REFINEMENTS=0
        theDIMS=np.array([self.R0,self.MIN_X,self.MAX_X,self.MIN_Y,self.MAX_Y],dtype=np.single)
        theNSIZES=np.array([self.N_TOTAL,self.N_INODES,self.N_REFINEMENTS],dtype=np.intc) # therefinements
        success=0
        
        fun.argtypes = [ndpointer(ctypes.c_float, ndim=1, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_int,ndim=1, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_float,ndim=1, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_int, ndim=1,flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_float,ndim=1, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_int, ndim=1,flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_float,ndim=1, flags="C_CONTIGUOUS"),
                        ctypes.c_int]
        # fun.argtypes=[ndpointer(ctypes.c_float, ndim=1, flags="C_CONTIGUOUS")]
        # fun(self.data)
        success=fun(self.data,self.NODE_TYPES,self.tgt_v,theNSIZES,theDIMS,ref_lengths,ref_v,MAX_ITER)


    
                           
    def set_bc(self, bdry, dof,value):
        self.constrained_edges.append([bdry,dof,value])
        
    def make_bcfile(self,tdofs):
      
        c_dofs=[]
        cvals=[]
        loc=self.filename.rfind(".")
        for i in self.constrained_edges:
            bnodes=np.where(self.NODE_TYPES==i[0])[0] 
            c_dofs+=(bnodes*tdofs+i[1]-1).tolist()# extract from tuple
            cvals+=[i[2]]*len(bnodes)
        rows=list(range(len(c_dofs)))
        with open(self.filename,"a+") as f:
            f.write("#"+self.filename[:loc]+":CONSTRAINT="+str(len(c_dofs))+"\n")
            for i in rows:
                f.write(str(i)+" ")
            f.write("\n")
            for i in c_dofs:
                f.write(str(i)+" ")
            f.write("\n")
            for i in c_dofs:
                f.write("1 ")
            f.write("\n")
            for i in cvals:
                f.write(str(i)+" ")
        f2=self.filename[:loc]+"_T"+self.filename[loc:]    
        with open(f2,"a+") as f:
            f.write("#"+self.filename[:loc]+":CONSTRAINT="+str(len(c_dofs))+"\n")
            for i in rows:
                f.write(str(i)+" ")
            f.write("\n")
            for i in c_dofs:
                f.write(str(i)+" ")
            f.write("\n")
            for i in c_dofs:
                f.write("1 ")
            f.write("\n")
            for i in cvals:
                f.write(str(i)+" ")
       
            
            
            
            
            