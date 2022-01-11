import numpy as np

def readFile(fName):
    while True:
        try:
            inputFile=open(fName,'r')
        except IOError:
            print("\n File not found.\n")
            continue
        else:
            break
    data = np.loadtxt(inputFile, delimiter=' ', dtype=np.intc,max_rows=1)
    N_SEGMENT_NODES=np.cumsum(data[1:])
    N_BDRY_SEGMENTS=data[0]
    N_BNODES=int(sum(data[1:]))

    BRm=np.loadtxt(inputFile, delimiter=' ', dtype=np.single, max_rows=N_BNODES, usecols=[0,1]) 
        #now find the injection locations 
    N_INJ=np.loadtxt(inputFile, delimiter=' ', dtype=np.intc,max_rows=1).item() 
        
    INJ_DATA=np.loadtxt(inputFile, delimiter=' ', dtype=np.single,max_rows=N_INJ) 
    if N_INJ==1:
        INJ_DATA=np.array([INJ_DATA.tolist()],dtype=np.single)  

        # now if you want to add some refinements here is the place to do it 
    N_REFINEMENTS=np.loadtxt(inputFile, delimiter=' ', dtype=np.intc,max_rows=1).item(0)
    REFINEMENT_DATA=[]
    if (N_REFINEMENTS>0):
        for i in range(N_REFINEMENTS):
            REFINEMENT_DATA.append(np.loadtxt(inputFile, delimiter=' ', dtype=np.single,max_rows=1))
    inputFile.close()
    constrained_edges=[]# needed for Boundary Condition Specification
    return BRm, N_SEGMENT_NODES, N_BNODES, N_INJ, INJ_DATA, N_REFINEMENTS, REFINEMENT_DATA