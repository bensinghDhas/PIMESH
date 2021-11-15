# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:32:28 2018

@author: arun.srinivasa
"""

""" 
NAMING CONVENTIONS
    prescripts 
        f_=myfunction
        o_=myobject instance
        c_=myclass
        in=input variable
        out=output variable
        the=global variable
        a_=a list iterator 


    postscripts
        v=np array vector
        m=np array matrix
        l=list 
        _c=iterator or counter
        s=list (plural)

    Variable names
        global variables will be descriptive camel notation starting with "the"
        localvariables all lower case with underscore where needed 
        all caps=gobal constant for that particular execution file ie the program will not alter it in any way
        if there is a name clash i will start my variables and functions with ars (arun r srinivasa)
        variables to be used as indexes to vectors will conform to the mathetics standard of i, j,k A,B,C etc for better readablity 
"""

# in and out (perhaps abbreviated to i and o? represent input and output variables in functions 

import numpy as np
import matplotlib.pyplot as pyp
from scipy.spatial import Delaunay
import matplotlib.cm as cm

def f_ars_triangulate(rs, inNodetypes,tgtm,filename):
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
    j=0
    # PLOTTING HERE_________________
    # N_bdrys=np.max(inNodetypes)
    # colors = cm.rainbow(np.linspace(0, 1, N_bdrys+1))
    for idx, i in zip(element_types,Elements):#{
        if(idx==0):#{
            j=j+1
            x[0:3]=xs[i]
            y[0:3]=ys[i]
            x[3],y[3]=xs[i[0]],ys[i[0]]
            pyp.plot(x,y,color='k')
    
    pyp.scatter(x=xs, y=ys, c=inNodetypes, cmap="rainbow")
    # for idx,x,y in zip(inNodetypes,xs,ys):
    #     pyp.plot(x,y,'.',color=colors[idx])
    # pyp.
    pyp.colorbar(label="Node Types", orientation="horizontal")
    
    #WRITING TO FILE___________________________________
    loc=filename.rfind('.')
    with open(filename,"w+") as f:
        f.write("#"+filename[:loc]+":NODES="+str(n_nodes)+', BOUNDARY_NODES='+str(n_bnodes)+"\n")
        for x,y  in zip(xs,ys):
            f.write(str(x)+' '+str(y)+'\n')
    #felem=filename[:loc]+'_elements'+filename[loc:]
    #with open(felem,"w+") as f:
        f.write("#"+filename[:loc]+":ELEMENTS="+str(j)+'\n')
        for ithelem,ithtri in zip(Elements, element_types):
            if(ithtri==0): # if it is a legit node
                f.write("1 "+str(ithelem[0])+' '+str(ithelem[1])+' '+str(ithelem[2])+'\n')
    
    f2=filename[:loc]+"_T"+filename[loc:]
    with open(f2,"w+") as f:
        f.write("#"+filename[:loc]+":NODES="+str(n_nodes)+', BOUNDARY_NODES='+str(n_bnodes)+"\n")
        for r  in rs:
            for i in r:
                f.write(str(i)+" ")
            f.write("\n")
    
        f.write("#"+filename[:loc]+":ELEMENTS="+str(j)+'\n')
        for ithelem,ithtri in zip(Elements, element_types):
            if(ithtri==0): # if it is a legit node
                f.write("1 "+str(ithelem[0])+' '+str(ithelem[1])+' '+str(ithelem[2])+'\n')
                
    
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
    #pyp.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
    
# def f_ars_truss(inNodes, inNodetypes,tgt_v,filename):
#     length=np.size(inNodes)
#     n_nodes=length//2
#     n_bnodes=sum(inNodetypes<0) # this is the number of boundary nodes
#     n_inodes=n_nodes-n_bnodes
#     pts=inNodes.reshape(n_nodes,2) # this is the data

#     #pts=inNodes.reshape(n_nodes,2) # this is the data
#     tri=Delaunay(pts)
#     Elements=tri.simplices# arrange the node numbers in ascending order.
#     #print(Elements[10:20,:])
#     #Elements.sort()
#     #print(Elements[10:20,:])
#     #Elements.sort()
#     N_TRIANGLES=len(Elements)
#     element_types=np.zeros(N_TRIANGLES, dtype=np.int16)
#     x=np.zeros(2)
#     y=np.zeros(2)
#     # for i in (Elements):
#     #     r01=pts[i[1],:]-pts[i[0],:]
#     #     r02=pts[i[2],:]-pts[i[0],:]
#     #     A=r01[0]*r02[1]-r01[1]*r02[0]
#     #     if (A<0):
#     #         temp=i[2]
#     #         i[2]=i[1]
#     #         i[1]=temp
        
#     for index, i in enumerate(Elements):#{
#         if(sum(inNodetypes[i]<0)>2): # no nodes outside so all three are boundary nodes 
#             centroid=sum(pts[i,:])/3
#             x[0:3]=pts[i,0]
#             y[0:3]=pts[i,1]
#             x[3],y[3]=pts[i[0],:]
#             for j in range(3):
#                 xij=x[j+1]-x[j]
#                 yij=y[j+1]-y[j]
#                 ny=tgt_v[2*i[j]]
#                 nx=tgt_v[2*i[j]+1]
#                 ny=-ny
#                 r=centroid-[x[j],y[j]]
#                 if(nx*r[0]+ny*r[1]>0):#{ # if centroid is outside the boundary
#                         element_types[index]=-1 
#     j=0
#     # PLOTTING HERE_________________
#     for idx, i in zip(element_types,Elements):#{
#         if(idx==0):#{
#             for jj in range(3):
#                 kk=(jj+1)%3
#                 s=i[jj]
#                 e=i[kk]
#                 bdry=(inNodetypes[s]!=0) and (inNodetypes[e]!=0)
#                 if (s<e  or bdry==True):
#                     j=j+1
#             # x[0:3]=pts[i,0]
#             # y[0:3]=pts[i,1]
#             # x[3],y[3]=pts[i[0],:]
#             #pyp.plot(x,y,'b')
#            # if (sum(inNodetypes[i]<0)==2):
#     #________________            
#             #}
#         #}
#     # write these to files
#     #Elements2=Elements.transpose()
#     # with open(filename,"w+") as f:
#     #     f.write("*TRUSS2D\n*SIZE\n")
#     #     f.write("DIM=2,ELEMENTS="+str(j)+",DOFS=2,NODES="+str(n_nodes)+"\n")
#     #     f.write("*NODES #nodal location and type: (type useful for different degrees of freedom or boundary conditions or forces etc) \n")
#     #     for ithpt,ithnode  in zip(pts,inNodetypes):
#     #         f.write(str(ithpt[0])+','+str(ithpt[1])+','+str(ithnode)+'\n')
#     #     f.write("*ELEMENTS\n")
#     #     for ithelem,ithtri in zip(Elements, element_types):
#     #         if(ithtri==0): # if it is a legit element
#     #             for jj in range(3):
#     #                 kk=(jj+1)%3
#     #                 s=ithelem[jj]
#     #                 e=ithelem[kk]
#     #                 bdry=(inNodetypes[s]!=0) and (inNodetypes[e]!=0)
#     #                 if (s<e  or bdry==True):
#     #                     f.write(str(s)+','+str(e)+',0'+'\n')
#     #                     x[0:2]=pts[[s,e],0]
#     #                     y[0:2]=pts[[s,e],1]
#     #                     pyp.plot(x,y,'b')
                        
         
         
    return tri