#import ..AMR as amr
from ..AMR import *
from ..NDTree import *
from ..Data import *

def make_example_tree(with_data=False):
    if with_data:
        t = Node(dim=2,data_class=SimpleTest2D,prolongate_func=prolongate_datafunc,
                restrict_func=restrict_datafunc)
    else:
        t = Node(dim=2)
    t.split()
    for c in t.child:
        c.split()
    t.child[0].child[0].split()
    t.child[0].child[3].split()
    t.child[1].child[0].split()
    t.child[1].child[1].split()
    t.child[3].child[0].split()
    return t
def make_example_tree2(with_data=False):
    if with_data:
        t = Node(dim=2,data_class=SimpleTest2D,prolongate_func=prolongate_datafunc,
                restrict_func=restrict_datafunc)
    else:
        t = Node(dim=2)
    t.split()
    for c in t.child:
        c.split()
    t.child[0].child[0].split()
    t.child[0].child[3].split()
    t.child[1].child[0].split()
    t.child[1].child[1].split()
    t.child[3].child[0].split()
    t.child[1].child[2].split()
    t.child[2].child[1].split()
    return t
def func(xc,yc):
    """Function which sets the data value"""
    cx = .65 + .5* 2.**(-8)
    cy = .65 +  .5* 2.**(-8)
    s = .1
    res = np.exp(-((xc-cx)**2+(yc-cy)**2)/(2*s**2))
    cx = .3 + .5* 2.**(-8)
    cy = .3 + .5* 2.**(-8)
    s = .1
    res += np.exp(-((xc-cx)**2+(yc-cy)**2)/(2*s**2))
    return res
class TestAMR():
#
#        t = Node(dim=2,prolongate_func=data.prolongate_datafunc,
#                restrict_func=data.restrict_datafunc,
#                 data = data.SimpleTest2D)
    def test_clear_refine(self):
        t = Node(dim=2)
        t.split()
        for n in t.list_leaves():
            n.rflag = True
        clear_refine(t)
        
        assert all([n.rflag == False for n in t.list_leaves()]) 
        
    def test_compression(self):
        import io
        import sys
        t = make_example_tree()
        capturedOutput = io.StringIO()                  # Create StringIO object
        sys.stdout = capturedOutput                     #  and redirect stdout.
        compression(t)                                     # Call function.
        sys.stdout = sys.__stdout__                     # Reset redirect.
        ans = '31 points out of 8^2 = 64 for full grid\nYou have saved a factor of 2.06\nWith a compression factor of 51.56%\n'
        assert ans == capturedOutput.getvalue()

    def test_neighbor_check(self):
        t = make_example_tree()
        t.child[3].child[0].child[0].rflag = True
        t.child[1].child[2].rflag= False
        t.child[2].child[1].rflag= False
        neighbor_check(t.child[3].child[0].child[0])
        assert t.child[1].child[2].rflag
        assert t.child[2].child[1].rflag

#    def test_refine(self):

#    def test_refinement_check(self):
#
    def test_refinement_lohner(self):
        
        # Neighbors are all on the same level
        t = make_example_tree2(with_data=True)
        leaf = t.child[3].child[0].child[0]
        total_neighbors = 3**leaf.dim
        offsets, neighbor_indices,neighbors, upper_neighbors = leaf.find_neighbors()
        final_list = [None]*total_neighbors

        for i in range(total_neighbors):

            if upper_neighbors[i] is not None:
                node = upper_neighbors[i]
                if not node.leaf:
                    node = neighbors[i]
                final_list[i] = node
        
   
        u = np.zeros((3,3))
        u[1,1] = func(.5,.5)
        u[0,1] = func(.375,.5)
        u[2,1] = func(.625,.5)
        u[1,0] = func(.5,.375)
        u[1,2] = func(.5,.625)
        u[0,0] = func(.375,.375)
        u[0,2] = func(.375,.625)
        u[2,0] = func(.625,.375)
        u[2,2] = func(.625,.625)
        au = abs(u)
        
        
        eps = 0 
        num = (u[0,1]-2*u[1,1]+u[2,1])**2 
        num += (u[1,0]-2*u[1,1]+u[1,2])**2
        num += (.25*(u[2,2]+u[0,0] - u[0,2]-u[2,0]))**2
        num += (.25*(u[2,2]+u[0,0] - u[0,2]-u[2,0]))**2
        den = (abs(u[1,1]-u[0,1])+abs(u[1,1]-u[2,1])+eps*(au[0,1]+2*au[1,1]+au[2,1]))**2
        den += (abs(u[1,1]-u[0,1])+abs(u[1,1]-u[2,1])+eps*(au[2,2]+au[0,0]+au[0,2]+au[2,0]))**2
        den += (abs(u[1,1]-u[1,0])+abs(u[1,1]-u[1,2])+eps*(au[1,0]+2*au[1,1]+au[1,2]))**2
        den += (abs(u[1,1]-u[1,0])+abs(u[1,1]-u[1,2])+eps*(au[2,2]+au[0,0]+au[0,2]+au[2,0]))**2
        ans = np.sqrt(num/den)
        res,value = refinement_lohner(leaf,final_list,tol=0,eps=0)
        return u,ans,value
        assert ans == value
        assert res == (ans>=0)
        res,value = refinement_lohner(leaf,final_list,tol=0.5,eps=0)
        assert res == (ans>=.5)
        eps = 0.01 
        num = (u[0,1]-2*u[1,1]+u[2,1])**2 
        num += (u[1,0]-2*u[1,1]+u[1,2])**2
        num += (.25*(u[2,2]+u[0,0] - u[0,2]-u[2,0]))**2
        num += (.25*(u[2,2]+u[0,0] - u[0,2]-u[2,0]))**2
        den = (abs(u[1,1]-u[0,1])+abs(u[1,1]-u[2,1])+eps*(au[0,1]+2*au[1,1]+au[2,1]))**2
        den += (abs(u[1,1]-u[0,1])+abs(u[1,1]-u[2,1])+eps*.25*(au[2,2]+au[0,0]+au[1,2]+au[2,1]))**2
        den += (abs(u[1,1]-u[1,0])+abs(u[1,1]-u[1,2])+eps*(au[1,0]+2*au[1,1]+au[1,2]))**2
        den += (abs(u[1,1]-u[1,0])+abs(u[1,1]-u[1,2])+eps*.25*(au[2,2]+au[0,0]+au[1,2]+au[2,1]))**2
        ans = np.sqrt(num/den)
        assert ans == value
        assert res == (ans>=.5)
        
        # 2 neighbors are on a coarser level
        
        
#        e0_num = 0 
#        e0_den = 0
#        eps = 0
#        u = np.zeros((3,3))
#        u[1,1] = func(.5,.5)
#        u[0,1] = func(.25,.5)
#        u[2,1] = func(.625,.5)
#        u[1,0] = func(.5,.25)
#        u[1,2] = func(.5,.625)
#        u[0,0] = func(.375,.375)
#        u[0,2] = func(.25,.5)
#        u[2,0] = func(.5,.25)
#        u[2,2] = func(.625,.625)
#
#    def test_start_derefine(self):
#
#    def test_start_refine(self):