import util as util
import pimesh as pm
import matplotlib.pylab as plt


fName="./examples/hole_bdrypts.txt"
bndryNodes,noSegNodes,noBndryNodes,N_INJ,INJ_DATA, N_REFINEMENTS, REFINEMENT_DATA=util.readFile(fName)
data=pm.initializeGeometry(bndryNodes,noSegNodes,noBndryNodes,N_INJ,INJ_DATA,N_REFINEMENTS, REFINEMENT_DATA)
pm.moveNodes(nodes,)

plt.plot(bndryNodes[:,0],bndryNodes[:,1],'o')
plt.show()