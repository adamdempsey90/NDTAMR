import numpy as np
import h5py
import matplotlib.pyplot as plt
#from .NDTree import Node


def grid_lines(node,slice=[0,1]):
    """
        Return the lines which split the node
    """
    if node.leaf:
        return None

    dx = 2.**(-node.global_index[0]-1)
    indx = np.array(node.global_index[1:])
    
    i,j = indx[slice]
    i_line = [ (dx*(2*i+1),dx*(2*j)),(dx*(2*i+1),dx*(2*(j+1)))]
    j_line = [ (dx*(2*i),dx*(2*j+1)),(dx*(2*(i+1)),dx*(2*j+1))]
    return [i_line,j_line]



def generate_grid(node,slice=[0,1],max_level=np.infty,save=None,xmin=None,xmax=None):
    if xmin is None:
        xmin = [0]*node.dim
    if xmax is None:
        xmax = [1]*node.dim

        
    xmin = np.array(xmin)
    xmax = np.array(xmax)

    lines = []
    node.walk(node_func=lambda x: lines.extend(grid_lines(x,slice=slice) if x.global_index[0]<max_level else [None,None]))
    
    xscale = xmax[slice]-xmin[slice]
    xstart = xmin[slice]
    grid = []
    for line in lines:
        if line is not None:
            grid.append( [
                (line[0][0]*xscale[0] + xstart[0], line[0][1]*xscale[1]+xstart[1]),
                (line[1][0]*xscale[0] + xstart[0],line[1][1]*xscale[1]+xstart[1])])



    if save is not None:
        np.array(grid).tofile(save)

    return grid

def grid_plot(node,slice=[0,1],max_level=np.infty,save=None,xmin=None,xmax=None,savefig=None,
             fig=None,ax=None,lw=1,colors='k',figsize=(6,6),**kargs):
    import matplotlib.collections as mc
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    if xmin is None:
        xmin = [0]*node.dim
    if xmax is None:
        xmax = [1]*node.dim

    xmin = np.array(xmin)
    xmax = np.array(xmax)

    grid = generate_grid(node,slice=slice,max_level=max_level,save=save,xmin=xmin,xmax=xmax)
    lc = mc.LineCollection(grid,colors=colors,lw=lw)

    ax.add_collection(lc)


    
    xmin = xmin[slice]
    xmax = xmax[slice]
    ax.set_xlim((xmin[0],xmax[0]))
    ax.set_ylim((xmin[1],xmax[1]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(slice[0]+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(slice[1]+1),fontsize=20)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    return fig,ax

def convert_to_uniform(tree,slice=[0,1],q=None,func=lambda x: x,**kargs):
    """Convert the tree to a numpy array for fast (and versitile) plotting"""
    lmax = tree.depth()
    
    res = np.zeros((2**lmax,2**lmax))
    
    for n in tree.list_leaves(attr='self'):
        lvl = n.global_index[0]
        i,j = np.array(n.global_index[1:])[slice]
        x,y = np.array(n.get_coords())[slice] 
        d = func(getattr(n.data,q))
        
        if lvl == lmax:
            res[i,j] = d
        else:
            fac = 2**(lmax-lvl)
            res[fac*i:fac*(i+1),fac*j:fac*(j+1)] = d
    return res
        
def plot(tree,slice=[0,1],q=None,cmap='viridis',rflag=False,func=lambda x: x,grid=False,figsize=(6,6),fig=None,ax=None,xmin=None,xmax=None,**kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    if xmin is None:
        xmin = [0]*tree.dim
    if xmax is None:
        xmax = [1]*tree.dim

    xmin = np.array(xmin)
    xmax = np.array(xmax)
    xmin1 = xmin[slice]
    xmax1 = xmax[slice]
    
    res = convert_to_uniform(tree,slice=slice,q=q,func=func)
    
    ax.imshow(res.T,extent=(xmin1[0],xmax1[0],xmin1[1],xmax1[1]),origin='lower',interpolation='none',cmap=cmap)

    if rflag:
        coords = []
        cfunc = lambda x: coords.append(x.get_coords(xmin=xmin,xmax=xmax,shift=True) if x.rflag else None)
        tree.walk(leaf_func=cfunc)
        for c in coords:
            if c is not None:
                ax.plot(c[slice[0]],c[slice[1]],'r.')
        
        

    
    
    ax.set_xlim((xmin1[0],xmax1[0]))
    ax.set_ylim((xmin1[1],xmax1[1]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(slice[0]+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(slice[1]+1),fontsize=20)
    ax.tick_params(labelsize=16)
    
    if grid:
        grid_plot(tree,slice=slice,fig=fig,xmin=xmin,xmax=xmax,ax=ax,**kargs)
    fig.tight_layout()
    return fig,ax
    
