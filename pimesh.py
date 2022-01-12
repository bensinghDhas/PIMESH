import numpy as np
import ctypes
import os
from scipy.spatial import Delaunay
from numpy.ctypeslib import ndpointer
import matplotlib.pylab as plt


def initializeGeometry(BRm,N_SEGMENT_NODES,N_BNODES,N_INJ,INJ_DATA, noRefine, RefineData):
    # create tangent vectors 
    TGTm=np.zeros_like(BRm)
    start=0 
    for i in N_SEGMENT_NODES:
        TGTm[start:i,:]=np.diff(BRm[start:i,:], axis=0,append=[BRm[start,:]])
        start=i
    refinedarea=find_refined_area(RefineData,noRefine)
       
    # #Area= xdy-ydx
    AREA=(sum(BRm[:,0]*TGTm[:,1]-BRm[:,1]*TGTm[:,0]))/2
    # Convert to unit tgt vector 
    BDRY_DISTv=np.sqrt(TGTm[:,0]**2+TGTm[:,1]**2)
    TGTm/=BDRY_DISTv[:,None]
    tgt_v=np.reshape(TGTm,2*N_BNODES)
    # max distance in the bdry nodes 
    R0= np.max(BDRY_DISTv)
        
    # Interior nodes:
    N_INODES=int((AREA+refinedarea)/(np.pi*(R0/2)**2))+20

    # Now assign the nodes and their properties 
    N_TOTAL=N_INODES+len(BRm) 
    print("interior nodes= ", N_INODES, "boundary nodes= ", N_BNODES," total nodes= ", N_TOTAL)
    ## better to store the nodes as x, y in a single array to save a lot of pain
    #_slices of data from storage________________________________________________
    data=np.zeros((6*N_TOTAL),dtype=np.single) # this contains all stored data in a contiguous list 
    data6=data.reshape(6,N_TOTAL)
    xs,ys,vxs, vys, fxs, fys=data6
    rs=data6[:2,:]
    xs[:N_BNODES],ys[:N_BNODES]=BRm[:,0],BRm[:,1]
    NODE_TYPES=np.zeros(N_TOTAL,dtype=np.intc)
    # assign the injected nodes 
        
    r_0= np.random.rand(N_INODES)*R0/2
    theta=np.random.rand(N_INODES)*np.pi*2
    cs=np.cos(theta)
    sine=np.sin(theta)
    j_start=N_BNODES
    xs[j_start:],  ys[j_start:] =r_0*cs, r_0*sine
    vxs[j_start:], vys[j_start:] =0.5*cs, 0.5*sine
     
    j=0
    for i in range(N_INODES):
        j=j%N_INJ
        xs[i+j_start]+=INJ_DATA[j,0]
        ys[i+j_start]+=INJ_DATA[j,1]
        j=j+1

    MAX=np.max(BRm,axis=0)+2*R0
    MIN=np.min(BRm,axis=0)-2*R0
    MIN_X,MAX_X=MIN[0],MAX[0]
    MIN_Y,MAX_Y=MIN[1],MAX[1]
    return data, tgt_v, NODE_TYPES, N_INODES, N_TOTAL, R0, MIN_X, MIN_Y, MAX_X, MAX_Y

def find_refined_area(REFINEMENT_DATA,N_REFINEMENTS):
    refinedarea=0
    if (N_REFINEMENTS>0):
        for refi in REFINEMENT_DATA:
            N=len(refi)-2 # first two contatin ht and sharpness data 
            ht=refi[0]
            X=refi[2:].reshape(int(N/2),2)
            tgt=np.diff(X,axis=0)
            refinedarea+=((sum(X[:-1,0]*tgt[:,1]-X[:-1,1]*tgt[:,0]))/2)*(1+ht)**2
    return refinedarea

def moveNodes(data, tgt_v, N_BNODES, NODE_TYPES, nodetype, REFINEMENT_DATA, N_TOTAL, N_INODES, N_REFINEMENTS,R0, MIN_X, MAX_X, MIN_Y, MAX_Y, fname, MAX_ITER):
    if (nodetype==-1):
        inputFile=open(fname,'r')
        NODE_TYPES[:N_BNODES]= np.loadtxt(inputFile, delimiter=' ', dtype=np.intc,max_rows=N_BNODES,usecols=[2],skiprows=3) 
        inputFile.close()
    else:
        NODE_TYPES[:N_BNODES]=nodetype 
        
    ref_lengths=np.array([0]+[len(x) for x in REFINEMENT_DATA],dtype=np.intc)
    ref_v=np.array([x for array in REFINEMENT_DATA for x in array],dtype=np.single)
    print(ref_lengths,ref_v)

    #C function call
    mydir=os.getcwd()+"/move_nodes_ctypes.so"
    my_functions=ctypes.cdll.LoadLibrary(mydir)
    fun =my_functions.f_moveNodes
    fun.restype = ctypes.c_double
        
    REFINEMENTS=0
    theDIMS=np.array([R0, MIN_X, MAX_X, MIN_Y, MAX_Y],dtype=np.single)
    theNSIZES=np.array([N_TOTAL, N_INODES, N_REFINEMENTS],dtype=np.intc) # therefinements
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
    success=fun(data, NODE_TYPES, tgt_v, theNSIZES,theDIMS,ref_lengths,ref_v,MAX_ITER)
    plt.plot(data[0:N_TOTAL],data[(N_TOTAL):(2*N_TOTAL)],'o')
    pos=np.column_stack((data[0:N_TOTAL],data[N_TOTAL:(2*N_TOTAL)]))
    # data=np.reshape(data[0:(2*N_TOTAL)],(N_TOTAL,2))
    # plt.plot(pos[:,0],pos[:,1],'o')
    # plt.show()
    tri=Delaunay(pos)
    return tri