# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:36:44 2021

@author: Arun

"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# def fig_close_handler(evt,in_fig_d,in_name): # thiere is no way that this can be in any object, it is the fucntion to close figures etc in tkinter
#     if (in_name in in_fig_d):
#         del(in_fig_d[in_name])
#     print("closed reference pic")

def f_ars_setlims(ax,in_minvals,in_maxvals):
    ax.set_xlim(in_minvals[0],in_maxvals[0]),   ax.set_ylim(in_minvals[1],in_maxvals[1])
    return

def f_ars_setaxes(in_fig,in_minvals, in_maxvals,in_ratio=1):
    DIM=len(in_minvals)
    xloc=-0.015
    yloc=-0.015
    if (DIM==2):
            ax=in_fig.gca()
            ax.set_xlim(in_minvals[0],in_maxvals[0]),   ax.set_ylim(in_minvals[1],in_maxvals[1])
            ax.set_xlabel('X(mm)', fontsize=8),   ax.set_ylabel('Y(mm)',fontsize=8)
            ax.xaxis.set_label_coords(0.95,yloc), ax.yaxis.set_label_coords(in_ratio*xloc, 0.95)
            ax.tick_params(axis="y",direction="in", pad=-22, which='major',labelsize=8)
            ax.tick_params(axis="x",direction="in", pad=-15, which='major',labelsize=8)
    else:
            #self.ax = self.Fig.gca(projection='3d')
            ax=Axes3D(in_fig)
            ax.set_xlim3d(in_minvals[0],in_maxvals[0]),   ax.set_ylim3d(in_minvals[1],in_maxvals[1]),   ax.set_zlim3d(in_minvals[2],in_maxvals[2])
            ax.set_xlabel('X(mm)',fontsize=6, fontweight='bold'),   ax.set_ylabel('Y(mm)',fontsize=6, fontweight='bold'), ax.set_zlabel('Z(mm)',fontsize=6, fontweight='bold')
            ax.xaxis.labelpad=-15
            ax.yaxis.labelpad=-15
            ax.zaxis.labelpad=-15
            #ax.xaxis.set_label_coords(0.95,yloc,0), ax.yaxis.set_label_coords(in_ratio*xloc, 0.95,0)
            ax.tick_params(axis="y",direction="in", pad=-15, which='major',labelsize=6)
            ax.tick_params(axis="x",direction="in", pad=-15, which='major',labelsize=6)
            ax.tick_params(axis="z",direction="in", pad=-15, which='major',labelsize=6)
    return ax

def f_ars_graph_loc(alpha, in_minvals,in_maxvals):
    DIM=len(alpha)
    loc=list((1-alpha[k])*in_minvals[k]+alpha[k]*in_maxvals[k] for k in range(DIM))
    return loc

def f_ars_round(val,n):# a beautifully clever way to round to n digits
    if  np.size(val)==1:
        ans= "0" if np.allclose(val,0) else "{:g}".format(float('{:.{p}g}'.format(val,p=n)))
    else:
        s=", "
        ans="{"+s.join([("0" if np.allclose(i,0) else'{:g}'.format(float('{:.1g}'.format(i)))) for i in val])+"}"
    return ans

def f_ars_myArrow(x,y, c,t, in_zorder, ax):
    xi=x[0]
    yi=y[0]
    dx=x[1]-xi
    dy=y[1]-yi
    return ax.arrow(xi,yi,dx, dy, head_width=6*t, head_length=10*t, fc=c, ec=c, width=t,length_includes_head=True,zorder=in_zorder)


class f_ars_Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def read_files(fp, switcher):
    # reads a file and splits into tokens
    # switcher={}
    # for i in token_list:
    #     switcher[i[0]]=i[1]
    # print(switcher)
    token=fp.readline()
    while(token!="*END"):
        if(token =="\n" or token[0]=="#"): # skip blank and comment lines
            token=fp.readline()# just get the next line 
            print(token)
        else:
            print("Reading",token)
            func=switcher.get(token, lambda :"*END")
            token=func()
    return

def read_number_list(fp,format_list=[float],array_type="list",rows=0,cols=0):
    token="*END"
    if (array_type=="list" or rows==0):
        storage_list=[]
        for my_line in  fp:
            indicator=my_line[0] # test if we have reached end of line
            if indicator=='*' or indicator=='\n' or indicator=='#':
                token=my_line
                break
            else:
                num_list=my_line.split(',')
                length=len(format_list)
                if (length==1):# if only one format is given 
                    num_list=list(map(format_list[0],num_list)) #apply to all elements 
                else:
                    for idx,numtype in enumerate(format_list[:-1]): # apply to the first list 
                        num_list[idx]=numtype(num_list[idx])
                num_list[length-1:]= list(map(format_list[-1],num_list[length-1:]))  # the final format is applied to the remaining elements  
                storage_list.append(num_list) # then add this to the list 
    elif (array_type=="array"):
        storage_list=np.loadtxt(fp, delimiter=",", dtype=format_list[0],max_rows=rows)
        if (rows==1):
            storage_list=np.array([storage_list.tolist()])
        line=fp.readline()
        token= line if(line) else "*END"
        
    elif(array_type=="array_list"): 
        storage_list=[]
        # create the list 
        for i,j in zip(cols,format_list):
            storage_list.append(np.zeros([rows,len(i)],dtype=j))
        # reading the lises
        idx=0 # this is the line number 
        for my_line in  fp:
            indicator=my_line[0] # test if we have reached end of line
            if indicator=='*' or indicator=='\n' or indicator=='#':
                token=my_line
                break
            else:
                num_list=my_line.split(',')
                for jdx,(fmt,icol) in enumerate(zip(format_list,cols)):
                    storage_list[jdx][idx,:]=[fmt(num_list[i]) for i in icol]
                idx+=1
                
    else:        
        raise Exception('support_function: read_number_list: type is neither a list nor array')
    return token,storage_list