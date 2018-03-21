import numpy as np
from .Vis import plot

def compression(tree):
    """Print some statistics about the efficiency of the AMR grid."""
    depth = tree.depth()
    nleaves = len(tree.list_leaves())
    
    tot = (2**depth)**tree.dim
    print('{:d} points out of {:d}^{:d} = {:d} for full grid'.format(nleaves,2**depth,tree.dim,tot))
    print('You have saved a factor of {:.0f}'.format(tot/nleaves))
    print('With a compression factor of {:.1f}%'.format((1-nleaves/tot)*100))
def clear_refine(tree):
    """Set all refinemnet flags to False"""
    tree.walk(leaf_func=lambda x: setattr(x,'rflag',False))

def start_refine(tree):
    """Look through leaves and split if flagged for refinement."""
    
    def do_split(x,count):
        if x.rflag:
            x.rflag = False
            count.append(x.name)
            x.split()
            
    total = []
    tree.walk(leaf_func = lambda x: do_split(x,total))
    return len(total)

def start_derefine(tree):
    """Look through leaves and derefine if needed."""
    
    def do_split(x,count):
        if x.rflag:
            count.append(x.name)
            x.remove()
            x.rflag = False
    total = []
    tree.walk(leaf_func = lambda x: do_split(x,total))
    return len(total)
def refine(tree,tol=.8,eps=.01,show=False,**kargs):
    depth = tree.depth()
    for lvl in range(depth+1)[::-1]:
        tree.walk(target_level=lvl,
                  leaf_func = lambda x: refinement_check(x,tol=tol,eps=eps,**kargs))
    print("Enforcing neighbors")
    for lvl in range(depth+1)[::-1]:
        tree.walk(target_level=lvl,leaf_func=lambda x: refinement_check(x,criteria=neighbor_check))
    if show:
        plot(tree,q='d',grid=True,rflag=True)
    
    total = start_refine(tree)
    return total


def neighbor_check(leaf,neighbors,**kargs):
    res = False
    for n in neighbors:  
        if n is not None:
            if leaf.global_index[0] < n.global_index[0]:    
                res |= n.rflag
    return res,0





def refinement_lohner(leaf,nodes,tol=.8,eps=.01,
                      min_value=1e-8,reverse=False,**kargs):

    total_neighbors = 3**leaf.dim
    ans = [False]*total_neighbors

    u = np.zeros((total_neighbors,))


    u1 = leaf.data.get_refinement_data()
#        u1 = clean_data(root.data.get_refinement_data())
#        unst = root.data.Bools.unst
#        inres = root.data.Bools.inres

#        if unst:
#            u1 = 180.
    for i,node in enumerate(nodes):
        if node is None:
            u[i] = u1
        else:
            try:
                d = node.data.get_refinement_data()
            except:
                print('Node',n,'failed on get_refinement_data()')
                print(d)
                raise
            #inres |= nodes[i][j].data.Bools.inres
            #if nodes[i][j].data.Bools.unst:
            #    d = u1

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
        denominator += (abs(u[iR]-u[iC]) + abs(u[iL]-u[iC]) + eps*abs(u[iL] + 2*u[iC] + u[iR]))**2
    #if corners:
    #numerator += (.5*abs( u[2,2] + u[0,0] - u[0,2] - u[2,0]))**2


    resx = np.sqrt(numerator/denominator)
    if abs(numerator) < min_value and abs(denominator) < min_value:
        resx = 0.
    if abs(denominator) < min_value:
        resx = 0.

    ans =  resx >= tol
    if reverse:
        ans = not ans
        
#        ans[iC] = True
#        for i in range(leaf.dim):
#            iL = [1]*leaf.dim
#            iR = [1]*leaf.dim
#            iL[i] += 1
#            iR[i] -= 1
#            ans[ifunc(iL)] = True
#            ans[ifunc(iR)] = True

    return ans,resx

def refinement_check(leaf,criteria=refinement_lohner,**kargs):
    """
        Check neighbors to see if this node should
        be refined.
    """

    # Get neighbors
  #  print(leaf.global_index[0],leaf.name,leaf.rflag)
    total_neighbors = 3**leaf.dim
    offsets, neighbor_indices,neighbors, upper_neighbors = leaf.find_neighbors()

    # Even if already tagged, still need to check new neighbors
    final_list = [None]*total_neighbors

    for i in range(total_neighbors):
#        ind = sum([j * 3**(2-k) for k,j in enumerate(i)])

        if upper_neighbors[i] is not None:
            node = upper_neighbors[i]
#            node = leaf.find(upper_neighbors[i])
            if not node.leaf:
                node = neighbors[i]
#                node = node.find(neighbors[i])
            final_list[i] = node


    res,value = criteria(leaf,final_list,**kargs)
    
    for node in final_list:
        if node is not None:
            node.rflag  |= res

    return res,value