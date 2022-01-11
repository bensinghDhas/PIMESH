# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 17:38:19 2021

@author: Arun
"""

import FlowMesherPythonBC as FM



wm=FM.FlowMesherPython("./examples/hole_bdrypts.txt")
wm.showgrid()

wm._move_node_c(3000)
# #wedgemesh._move_node_function(100)
wm.triangulate()
# fin=time.time()
wm.set_bc(bdry=1,dof=2,value=0)
wm.set_bc(bdry=1,dof=1,value=0)
wm.set_bc(bdry=3,dof=2,value=1)
wm.make_bcfile(tdofs=2)