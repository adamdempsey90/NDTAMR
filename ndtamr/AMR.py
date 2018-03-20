import numpy as np

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
            count.append(x.name)
            x.split()
            x.rflag = False
    total = []
    tree.walk(leaf_func = lambda x: do_split(x,total))
    return len(total)
def refine(tree,tol=.8,eps=.01,maxdiff=1,**kargs):
    depth = tree.depth()
    for lvl in range(depth+1)[::-1]:
        tree.walk(target_level=lvl,
                  leaf_func = lambda x: refinement_check(x,tol=tol,eps=eps,**kargs))
    if maxdiff > 0:
        tree.walk(leaf_func = lambda x: maxlevel_check(x,maxdiff,**kargs))
    total = start_refine(tree)
    return total
    
def refinement_check(leaf,refine_all=False,
                     corners=False,**kargs):
    """
        Check neighbors to see if this node should
        be refined.
    """

    # Get neighbors

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


    res,num,den,result = refinement_lohner(leaf,final_list,**kargs)
    

    for i in range(total_neighbors):
        if final_list[i] is not None:
            final_list[i].rflag |= res

    return num,den,result

def maxlevel_check(leaf,refine_all=False,
                     corners=False,**kargs):
    """
        Check for a neighbor that is being refined
    """

    # Ignore if already flagged for refinement
    if leaf.rflag:
        return
    
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


    res,num,den,result = refinement_lohner(leaf,final_list,**kargs)
    

    for i in range(total_neighbors):
        if final_list[i] is not None:
            final_list[i].rflag |= res[i]

    return num,den,result

def refinement_lohner(leaf,nodes,tol=.8,eps=.01,
                      min_value=1e-8,**kargs):

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
#        ans[iC] = True
#        for i in range(leaf.dim):
#            iL = [1]*leaf.dim
#            iR = [1]*leaf.dim
#            iL[i] += 1
#            iR[i] -= 1
#            ans[ifunc(iL)] = True
#            ans[ifunc(iR)] = True

    return ans,np.sqrt(numerator),np.sqrt(denominator),resx