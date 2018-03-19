import numpy as np
import h5py
from .NDTree import Node


def grid_lines(self,i1=0,i2=1):
    """
        Return the lines which split the node
    """
    if self.leaf:
        return None

    dx = 2.**(-self.global_index[0])
    indx = self.global_index[1:]
    dx /= 2
    i = indx[i1]
    j = indx[i2]
    i_line = [ (dx*(2*j+1),dx*(2*i)),(dx*(2*j+1),dx*(2*(i+1)))]
    j_line = [ (dx*(2*j),dx*(2*i+1)),(dx*(2*(j+1)),dx*(2*i+1))]
    return [i_line,j_line]



def generate_grid(self,i1=0,i2=1,max_level=np.infty,save=None,xmin=None,xmax=None):
    if xmin is None:
        xmin = [0]*self.dim
    if xmax is None:
        xmax = [1]*self.dim


    lines = []
    self.walk(node_func=lambda x: lines.extend(x.grid_lines(i1=i1,i2=i2) if x.global_index[0]<max_level else [None,None]))
    xscale = (xmax[i1]-xmin[i1])
    yscale = xmax[i2]-xmin[i2]
    xstart = xmin[i1]
    ystart = xmin[i2]
    grid = []
    for line in lines:
        if line is not None:
            grid.append( [
                (line[0][0]*xscale + xstart, line[0][1]*yscale+ystart),
                (line[1][0]*xscale + xstart,line[1][1]*yscale+ystart)])



    if save is not None:
        np.array(grid).tofile(save)

    return grid

def grid_plot(self,i1=0,i2=1,max_level=np.infty,save=None,xmin=None,xmax=None,savefig=None,
             lw=1,colors='k',figsize=(6,6)):
    import matplotlib.collections as mc
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(figsize=figsize)

    if xmin is None:
        xmin = [0]*self.dim
    if xmax is None:
        xmax = [1]*self.dim


    grid = self.generate_grid(i1=i1,i2=i2,max_level=max_level,save=save,xmin=xmin,xmax=xmax)
    lc = mc.LineCollection(grid,colors=colors,lw=lw)

    ax.add_collection(lc)


    ax.set_xlim((xmin[i1],xmax[i1]))
    ax.set_ylim((xmin[i2],xmax[i2]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(i1+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(i2+1),fontsize=20)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    return fig,ax