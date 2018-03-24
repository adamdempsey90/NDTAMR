import numpy as np
import h5py
import matplotlib.pyplot as plt
#from .NDTree import Node


def grid_lines(node,dims=[0,1],slice=None):
    """
        Return the lines which split the node
    """
    if node.leaf:
        return None

    dx = 2.**(-node.global_index[0]-1)
    indx = np.array(node.global_index[1:])
    
    coords = np.array(node.coords)
    ndx = np.array(node.dx)
    if slice is not None:
        for s in slice:
            if not (s[1] >= coords[s[0]] and s[1] < coords[s[0]]+ndx[s[0]]):
                return [None,None]
            
    i,j = indx[dims]
    idx,jdx = np.array(node.dx)[dims] / 2
    istart,jstart = np.array(node.xmin)[dims]
    i_line = [ (istart + idx*(2*i+1),jstart+jdx*(2*j)),(istart+idx*(2*i+1),jstart+jdx*(2*(j+1)))]
    j_line = [ (istart + idx*(2*i),jstart+jdx*(2*j+1)),(istart+idx*(2*(i+1)),jstart+jdx*(2*j+1))]
    return [i_line,j_line]



def generate_grid(node,dims=[0,1],slice=None,max_level=np.infty,save=None,):
    lines = []
    node.walk(node_func=lambda x: lines.extend(grid_lines(x,dims=dims,slice=slice) if x.global_index[0]<max_level else [None,None]))
    
    grid = []
    for line in lines:
        if line is not None:
            grid.append( [
                (line[0][0], line[0][1]),
                (line[1][0],line[1][1])])



    if save is not None:
        np.array(grid).tofile(save)

    return grid

def grid_plot(node,dims=[0,1],slice=None,max_level=np.infty,save=None,savefig=None,
             fig=None,ax=None,lw=1,colors='k',figsize=(6,6),**kargs):
    import matplotlib.collections as mc
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)


    xmin = np.array(node.xmin)
    xmax = np.array(node.xmax)

    grid = generate_grid(node,dims=dims,slice=slice,max_level=max_level,save=save)
    lc = mc.LineCollection(grid,colors=colors,lw=lw)

    ax.add_collection(lc)


    
    xmin = xmin[dims]
    xmax = xmax[dims]
    ax.set_xlim((xmin[0],xmax[0]))
    ax.set_ylim((xmin[1],xmax[1]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(dims[0]+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(dims[1]+1),fontsize=20)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    return fig,ax

def convert_to_uniform(tree,dims=[0,1],slice=None,q=None,func=lambda x: x,**kargs):
    """Convert the tree to a numpy array for fast (and versitile) plotting.
        slice = [(dimension,value)] 
    """
    
   
    if tree.dim > len(dims) and slice is None:
        slice = [(-1,0)]
    xmin = np.array(tree.xmin)
    xmax = np.array(tree.xmax)
    
    lmax = tree.depth()
    
    res = np.zeros((2**lmax,2**lmax))
    leaves = tree.list_leaves(attr='self')
    for n in leaves:
        lvl = n.global_index[0]
        indices = np.array(n.global_index[1:])
        dx = np.array(n.dx)
        coords = np.array(n.coords)
        i,j = indices[dims]
        
        if slice is None:
            d = func(getattr(n.data,q))
            if lvl == lmax:
                res[i,j] = d
            else:
                fac = 2**(lmax-lvl)
                res[fac*i:fac*(i+1),fac*j:fac*(j+1)] = d
        else:
            good = all([s[1] >= coords[s[0]] and s[1] < coords[s[0]]+dx[s[0]] for s in slice])
            if good:
                d = func(getattr(n.data,q))
                if lvl == lmax:
                    res[i,j] = d
                else:
                    fac = 2**(lmax-lvl)
                    res[fac*i:fac*(i+1),fac*j:fac*(j+1)] = d
    return res

def get_slice(tree,dim,q,func,slice):
    if tree.dim > 1 and slice is None:
        slice = [(-1,0)]
    def _lfunc(n,dim,slice,q,func):
        lvl = n.global_index[0]
        indices = np.array(n.global_index[1:])
        dx = np.array(n.dx)
        coords = np.array(n.coords)
        
        good = all([s[1] >= coords[s[0]] and s[1] < coords[s[0]]+dx[s[0]] for s in slice])
        if good:
            return coords[dim],func(getattr(n.data,q))
        return None,None
    vals = []
    tree.walk(leaf_func = lambda x: vals.append(_lfunc(x,dim,slice,q,func)))
    return list(filter(lambda x: not None in x,vals))

def line_plot(tree,dim=0,slice=None,q=None,func=lambda x: x,fig=None,ax=None,**kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=(8,6))
        
    if tree.dim == 1:
        vals=[]
        tree.walk(leaf_func = lambda x: vals.append([x.coords[0], func(getattr(x.data,q))]))    
        vals = np.array(vals)
    else:
        vals = np.array(get_slice(tree,dim,q,func,slice))

        

    ax.plot(vals[:,0],vals[:,1],**kargs)
    
    ax.set_xlabel('$x$')
    ax.set_ylabel(q)
    ax.minorticks_on()
    fig.tight_layout()
    return fig,ax
def plot(tree,dims=[0,1],slice=None,q=None,cmap='viridis',rflag=False,func=lambda x: x,grid=False,figsize=(6,6),fig=None,ax=None,**kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)


    xmin = np.array(tree.xmin)
    xmax = np.array(tree.xmax)
    xmin1 = xmin[dims]
    xmax1 = xmax[dims]
    
    res = convert_to_uniform(tree,dims=dims,slice=slice,q=q,func=func)
    
    ax.imshow(res.T,extent=(xmin1[0],xmax1[0],xmin1[1],xmax1[1]),origin='lower',interpolation='none',cmap=cmap)
    
    _create_colorbar(ax,vmin=res.min(),vmax=res.max(),cmap=cmap)

    if rflag:
        coords = []
        cfunc = lambda x: coords.append(x.get_coords(shift=True) if x.rflag else None)
        tree.walk(leaf_func=cfunc)
        for c in coords:
            if c is not None:
                ax.plot(c[dims[0]],c[dims[1]],'r.')
        
        

    
    
    ax.set_xlim((xmin1[0],xmax1[0]))
    ax.set_ylim((xmin1[1],xmax1[1]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(dims[0]+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(dims[1]+1),fontsize=20)
    ax.tick_params(labelsize=16)
    
    if grid:
        grid_plot(tree,dims=dims,slice=slice,fig=fig,ax=ax,**kargs)
    fig.tight_layout()
    return fig,ax
def contour(tree,Nconts=20,dims=[0,1],slice=None,q=None,cmap='viridis',rflag=False,func=lambda x: x,grid=False,colors='grey',figsize=(6,6),fig=None,ax=None,**kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)


    xmin = np.array(tree.xmin)
    xmax = np.array(tree.xmax)
    xmin1 = xmin[dims]
    xmax1 = xmax[dims]
    
    res = convert_to_uniform(tree,dims=dims,slice=slice,q=q,func=func)
    
    vmin = res.min()
    vmax = res.max()
    
    ax.contour(res.T,Nconts,extent=(xmin1[0],xmax1[0],xmin1[1],xmax1[1]),origin='lower',
               vmin=vmin,vmax=vmax,cmap=cmap,**kargs)
    
    _create_colorbar(ax,vmin=vmin,vmax=vmax,cmap=cmap)

  

    
    
    ax.set_xlim((xmin1[0],xmax1[0]))
    ax.set_ylim((xmin1[1],xmax1[1]))

    ax.minorticks_on()
    ax.set_xlabel('$x_{:d}$'.format(dims[0]+1),fontsize=20)
    ax.set_ylabel('$x_{:d}$'.format(dims[1]+1),fontsize=20)
    ax.tick_params(labelsize=16)
    
    if grid:
        grid_plot(tree,dims=dims,colors=colors,slice=slice,fig=fig,ax=ax,**kargs)
    fig.tight_layout()
    return fig,ax 
def _create_colorbar(ax,vmin,vmax,log=False,cmap='viridis',**kargs):
    import matplotlib
    import matplotlib.cm
    import matplotlib.colors as colors
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('top',size='3%',pad=.05)
    if log:
        norm = colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = matplotlib.cm.get_cmap(cmap)
    cb = matplotlib.colorbar.ColorbarBase(ax=cax,cmap=cmap,norm=norm,orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=12)
    return cb
