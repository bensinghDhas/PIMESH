import numpy as np
import ctypes
import os
from scipy.spatial import Delaunay

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

    # MAX=np.max(self.BRm,axis=0)+2*self.R0
    # self.MIN=np.min(self.BRm,axis=0)-2*self.R0
    # self.MIN_X,self.MAX_X=self.MIN[0],self.MAX[0]
    # self.MIN_Y,self.MAX_Y=self.MIN[1],self.MAX[1]
    return xs, ys

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

def moveNodes(self, N_BNODES,NODE_TYPES, nodetype,REFINEMENT_DATA, fname, MAX_ITER):
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
    return tri

def triangulate(rs, inNodetypes,tgtm,filename):
    xs=rs[0,:]
    ys=rs[1,:]
    n_nodes=np.size(xs)
    n_bnodes=sum(inNodetypes>0) # this is the number of boundary nodes
    n_inodes=n_nodes-n_bnodes
    tri=Delaunay(rs.T)
    Elements=tri.simplices# arrange the node numbers in ascending order.
    N_TRIANGLES=len(Elements)
    element_types=np.zeros(N_TRIANGLES, dtype=np.int16)
    x=np.zeros(4)
    y=np.zeros(4)
    # for i in (Elements):
    #     r01=pts[i[1],:]-pts[i[0],:]
    #     r02=pts[i[2],:]-pts[i[0],:]
    #     A=r01[0]*r02[1]-r01[1]*r02[0]
    #     if (A<0):
    #         temp=i[2]
    #         i[2]=i[1]
    #         i[1]=temp
    centroid=np.zeros(2,dtype=np.single)    
    for index, i in enumerate(Elements):#{
        if(sum(inNodetypes[i]>0)>2): # no nodes outside so all three are boundary nodes 
            centroid[:]=[sum(xs[i])/3,sum(ys[i]/3)]
            x[0:3]=xs[i]
            y[0:3]=ys[i]
            for j in range(3):
                ny=tgtm[i[j],0]
                nx=tgtm[i[j],1]
                ny=-ny
                r=centroid-[x[j],y[j]]
                if(nx*r[0]+ny*r[1]>0):#{ # if centroid is outside the boundary
                        element_types[index]=-1 
    # j=0
    # # PLOTTING HERE_________________
    # # N_bdrys=np.max(inNodetypes)
    # # colors = cm.rainbow(np.linspace(0, 1, N_bdrys+1))
    # for idx, i in zip(element_types,Elements):#{
    #     if(idx==0):#{
    #         j=j+1
    #         x[0:3]=xs[i]
    #         y[0:3]=ys[i]
    #         x[3],y[3]=xs[i[0]],ys[i[0]]
    #         pyp.plot(x,y,color='k')
    
    # pyp.scatter(x=xs, y=ys, c=inNodetypes, cmap="rainbow")
    # # for idx,x,y in zip(inNodetypes,xs,ys):
    # #     pyp.plot(x,y,'.',color=colors[idx])
    # # pyp.
    # pyp.colorbar(label="Node Types", orientation="horizontal")
    
    # #WRITING TO FILE___________________________________
    # loc=filename.rfind('.')
    # with open(filename,"w+") as f:
    #     f.write("#"+filename[:loc]+":NODES="+str(n_nodes)+', BOUNDARY_NODES='+str(n_bnodes)+"\n")
    #     for x,y  in zip(xs,ys):
    #         f.write(str(x)+' '+str(y)+'\n')
    # #felem=filename[:loc]+'_elements'+filename[loc:]
    # #with open(felem,"w+") as f:
    #     f.write("#"+filename[:loc]+":ELEMENTS="+str(j)+'\n')
    #     for ithelem,ithtri in zip(Elements, element_types):
    #         if(ithtri==0): # if it is a legit node
    #             f.write("1 "+str(ithelem[0])+' '+str(ithelem[1])+' '+str(ithelem[2])+'\n')
    
    # f2=filename[:loc]+"_T"+filename[loc:]
    # with open(f2,"w+") as f:
    #     f.write("#"+filename[:loc]+":NODES="+str(n_nodes)+', BOUNDARY_NODES='+str(n_bnodes)+"\n")
    #     for r  in rs:
    #         for i in r:
    #             f.write(str(i)+" ")
    #         f.write("\n")
    
    #     f.write("#"+filename[:loc]+":ELEMENTS="+str(j)+'\n')
    #     for ithelem,ithtri in zip(Elements, element_types):
    #         if(ithtri==0): # if it is a legit node
    #             f.write("1 "+str(ithelem[0])+' '+str(ithelem[1])+' '+str(ithelem[2])+'\n')
                
    
    # with open(filename,"w+") as f:
    #     f.write('# Number of Nodes, Number of Boundary Nodes and Number of Elements\n')
    #     f.write(str(n_nodes)+','+str(n_bnodes)+','+str(j)+'\n')
    #     f.write("#nodal location and type: (type useful for different degrees of freedom or boundary conditions or forces etc) ")
    #     for ithpt,ithnode  in zip(pts,inNodetypes):
    #         f.write(str(ithpt[0])+','+str(ithpt[1])+','+str(ithnode)+'\n')
    #     f.write("#Elements connectivity")
    #     for ithelem,ithtri in zip(Elements, element_types):
    #         if(ithtri==0): # if it is a legit node
    #             f.write(str(ithelem[0])+','+str(ithelem[1])+','+str(ithelem[2])+'\n')
             
         
         
    return tri