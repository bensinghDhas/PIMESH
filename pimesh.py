import numpy as np

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
