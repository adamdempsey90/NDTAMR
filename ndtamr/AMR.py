"""
    The AMR module contains functions which adaptively refine the domain.
"""
from __future__ import print_function, division
import numpy as np
from .Vis import plot

def compression(tree):
    """
    Print some statistics about the efficiency of the AMR grid.
    """
    depth = tree.depth()
    nleaves = len(tree.list_leaves())
    
    tot = (2**depth)**tree.dim
    print('{:d} points out of {:d}^{:d} = {:d} for full grid'.format(nleaves,2**depth,tree.dim,tot))
    print('You have saved a factor of {:.0f}'.format(tot/nleaves))
    print('With a compression factor of {:.1f}%'.format((1-nleaves/tot)*100))
    
def clear_refine(tree):
    """
    Set all refinemnet flags to False
    """
    tree.walk(leaf_func=lambda x: setattr(x,'rflag',False))

def start_refine(tree):
    """
    Look through leaves and split if flagged for refinement.
    """
    
    def _do_split(node,count):
        """
        Split the node if it is flagged for refinement.

        Parameters
        ----------
        node : NDtree.Node
            Node we are evaluating
        count : list
            The list of nodes which have been refined
        """
        if node.rflag:
            node.rflag = False
            count.append(node.name)
            node.split()
            
    total = []
    tree.walk(leaf_func = lambda x: _do_split(x,total))
    return len(total)

def start_derefine(tree):
    """
    Look through leaves and derefine if needed.
    """
    
    def _do_unsplit(x,count):
        """
        Unsplit the node if it is flagged for derefinement.

        Parameters
        ----------
        node : NDtree.Node
            Node we are evaluating
        count : list
            The list of nodes which have been derefined
        """
        if x.rflag:
            count.append(x.name)
            x.remove()
            x.rflag = False
    total = []
    tree.walk(leaf_func = lambda x: _do_unsplit(x,total))
    return len(total)
def refine(tree,tol=.8,eps=.01,show=False,**kargs):
    """
    The main AMR routine which evaluates and refines 
    each node in the tree.

    Parameters
    ----------
    tree : NDTree.node
        The tree we want to refine
    tol : float
        The tolerance level for refinement
    eps : float
        Helps with smoothing out flucuations in the 
        refinement variable
    show : bool
        If True plot the domain showing which cells 
        will be refined
    **kargs :
        Keyword arguments passed to the refinement_check function
    
    Returns
    -------
    total : int
        The total number of refined cells

    """
    depth = tree.depth()
    for lvl in range(depth+1)[::-1]:
        tree.walk(target_level=lvl,
                  leaf_func = lambda x: refinement_check(x,tol=tol,eps=eps,**kargs))
    for lvl in range(depth+1)[::-1]:
        tree.walk(target_level=lvl,leaf_func = neighbor_check)
    if show:
        plot(tree,q='d',grid=True,rflag=True)
    
    total = start_refine(tree)
    return total


def neighbor_check(node):
    """
    Check that if a coarser neighbor refined we also refine.
    This enforces that the maximum discrepency in neighbor levels
    is one.
    """
    if not node.rflag:
        return
    _,_,_,neighbors = node.find_neighbors()
    
    for n in neighbors:
        if n is not None:
            if n.leaf:
                n.rflag = True
    
    
   
def refinement_lohner(leaf,nodes,tol=.8,eps=.01,
                      min_value=1e-5,reverse=False,corners=True,**kargs):
    """
    The refinement criteria of L\"{o}hner (1987).
    This function does not evaulate neighbors which are on a finer level, as
    they should have already been evaulated.
    The user has the option of including the cross derivative terms with the 
    corners keyword argument.

    Parameters
    ----------
    leaf : NDTree.node
        The leaf node we are evaluating
    nodes : list
        List of neighbor leaves.
    tol : float
        The tolerance level for refinement
    eps : float
        Helps with smoothing out flucuations in the 
        refinement variable
    min_value : float
        If the second derivative values are below this value then
        we do not refine.
    reverse : bool
        If True then we flag the cell if it does not satisfy the 
        refinement criteria
    corners : bool
        If True we include the "corner" cells in the evauluation
    Returns
    -------
    res: bool
        If True we refine this cell.
    value: float
        The numerical value for the refinement criteria
        
    """

    total_neighbors = 3**leaf.dim

    u = np.zeros((total_neighbors,))


    u1 = leaf.data.get_refinement_data()

    for i,node in enumerate(nodes):
        if node is None:
            u[i] = u1
        else:
            if not node.leaf:
                d = node.restrict().get_refinement_data()
            else:
                d = node.data.get_refinement_data()
#            try:
#                d = node.data.get_refinement_data()
#            except:
#                print('Node',node,'failed on get_refinement_data()')
#                raise

            u[i] = d

    numerator = 0
    denominator = 0


    ifunc = lambda x: sum([j * 3**(leaf.dim-1-k) for k,j in enumerate(x)])

    iC = ifunc([1]*leaf.dim)

    for i in range(leaf.dim):
        iL = [1]*leaf.dim
        iR = [1]*leaf.dim
        iL[i] += 1
        iR[i] -= 1
        
        iL = ifunc(iL)
        iR = ifunc(iR)
        numerator += (u[iR] - 2*u[iC] + u[iL])**2
        if corners:
            for j in range(i+1,leaf.dim):
                jRR = [1]*leaf.dim
                jLL = [1]*leaf.dim
                jLR = [1]*leaf.dim
                jRL = [1]*leaf.dim
                jRR[i] += 1
                jRR[j] += 1
                jLL[i] -= 1
                jLL[j] -= 1
                jLR[i] -= 1
                jLR[j] += 1
                jRL[i] += 1
                jRL[j] -= 1
                jRR = ifunc(jRR)
                jLL = ifunc(jLL)
                jRL = ifunc(jRL)
                jLR = ifunc(jLR)
                numerator += (.5*( u[jRR] + u[jLL] - u[jRL] - u[jLR]))**2
            
        denominator += (abs(u[iR]-u[iC]) + abs(u[iL]-u[iC]) + eps*(abs(u[iL]) -2*abs(u[iC]) + abs(u[iR])))**2

    if abs(denominator) < min_value or abs(numerator) < min_value:
        value = 0.
    else:
        value = np.sqrt(numerator/denominator)
    

    res =  value >= tol
    if reverse:
        res = not res
        

    return res,value

def refinement_check(leaf,criteria=refinement_lohner,**kargs):
    """
    Deterimine neighbors and see if this node should be refined.
    If the node satisfies the criteria, then we also flag all of its
    leaf neighbors.

    Parameters
    ----------
    leaf : NDTree.node
        The leaf node we are evaluating
        
    criteria : function
        The function which evaluates the refinement criteria.
    **kargs :
        Keyword arguments which are passed to the criteria function.

    Returns
    -------
    res: bool
        If True we refine this cell.
    value: float
        The numerical value for the refinement criteria

    """

    # Get neighbors
    total_neighbors = 3**leaf.dim
    offsets, neighbor_indices,neighbors, upper_neighbors = leaf.find_neighbors()

    # Even if already tagged, still need to check new neighbors
    final_list = [None]*total_neighbors

    for i in range(total_neighbors):

        if upper_neighbors[i] is not None:
            node = upper_neighbors[i]
            if not node.leaf:
                node = neighbors[i]
            final_list[i] = node


    res,value = criteria(leaf,final_list,**kargs)
    
    for node in final_list:
        if node is not None:
            if node.leaf:
                node.rflag  |= res
            

    return res,value
