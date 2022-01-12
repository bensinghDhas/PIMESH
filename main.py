import util as util
import pimesh as pm
import matplotlib.pylab as plt

MAX_ITER=3000
nodetype=-1
fName="./examples/hole_bdrypts.txt"
bndryNodes,noSegNodes,noBndryNodes,N_INJ,INJ_DATA, N_REFINEMENTS, REFINEMENT_DATA=util.readFile(fName)
data,tgt_v, NODE_TYPES, N_INODES, N_TOTAL, R0, MIN_X, MIN_Y, MAX_X, MAX_Y=pm.initializeGeometry(bndryNodes,noSegNodes,noBndryNodes,N_INJ,INJ_DATA,N_REFINEMENTS, REFINEMENT_DATA)
tri=pm.moveNodes(data, tgt_v, noBndryNodes, NODE_TYPES, nodetype, REFINEMENT_DATA, N_TOTAL, N_INODES, N_REFINEMENTS,R0, MIN_X, MAX_X, MIN_Y, MAX_Y, fName, MAX_ITER)
plt.plot(tri.points[:,0],tri.points[:,1],'o')
plt.show()