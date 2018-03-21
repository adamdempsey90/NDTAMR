import numpy as np
import h5py
#from .DataClass import Sim as Data


class Node():
    """
        Nodes either point to their children or they have no
        children and instead hold some data.
    """
    def __init__(self,name='0x0',dim=2,file=None,parent=None,data=None):
        """
            name is the hexadecimal chain of children indices that
            traces the node back to the root
            For example, the 5th child of the root would be called
            0x00x5, and the 7th child of the 3rd child of the root
            would be 0x00x30x7
        """
        self.dim = dim
        self.fmt = '0{:d}b'.format(dim)

        self.name = name
        self.global_index = (0,) + tuple([0]*dim)
        self.parent = parent
        self.leaf = True
        self.rflag = False
        self.nchildren = 2**dim
        self.child = [None]*self.nchildren
        self.file = file
        self.data = data


    #   Take binary number and return the index relative to parent
        self.index_from_bin = lambda x: tuple(map(int,x))

    #  Take the child index and convert it to binary
    #  Ex: 2 --> '010'
    #      6 --> '110'
        self.tobin = lambda x: format(x,self.fmt)

    #  Take a binary number and convert it to an integer
        self.frombin = lambda x: int(x,base=2)

        self.child_index ={i:self.index_from_bin(self.tobin(i)) for i in range(self.nchildren)}

        self.global_index = self.get_global_index(name)
#    def save(self,file):
#        """
#            Write this node to the hdf5 group/file.
#        """
#        grp = file.create_group(self.indx[-1])
#        if self.leaf:
#            # We are a leaf, so we should dump our data
#            dset = grp.create_group('Data')
#            self.data.save(dset)
#        else:
#            # We are not a group, so call the children
#            for c in self.child:
#                c.save(grp)
#        return
#    def build(self,f):
#        """
#            Look in the hdf5 group f for child cells
#        """
#        for i in range(self.nchildren):
#            try:
#                grp = f[hex(i)]
#                self.leaf = False
#                self.child[i] = Node(self.name+hex(i),parent=self,file=grp)
#                self.child[i].build(grp)
#            except KeyError:
#                self.leaf = True
#                self.datastr = self.name + '/' + hex(i) + '/Data'
#                self.data=Data(fname=f['Data'],node=self)
#        return
    def get_local_index(self,name):
        """
            Get the local index relative to the parent from the name
        """
        try:
            indx = int(name,base=16)
        except ValueError:
            indx = int( name.split('0x')[-1],base=16)
        binary_name = self.tobin(indx)
        return self.index_from_bin(binary_name)

    def get_global_index(self,name):
        """
            Calculate the global index from the name
        """
        glindx = [0]*self.dim

        names = name.split('0x')[1:]

        level = len(names)-1
        #print('Global',glindx)
        for n in name.split('0x')[1:]:
            lindx = self.get_local_index(n)
           # print('Local',lindx)
            glindx = [2*g+i for g,i in zip(glindx,lindx)]
           # print('New global',glindx)
        return (level,) + tuple(glindx)
    def move_index_up(self,indx):
        """
            Move an index up a level, returning its name and parent index
        """
        pindx = [i//2 for i in indx]
        name = ''.join(map(str,[i%2 for i in indx]))
        name = hex( self.frombin(name) )

        return pindx, name

    def get_level(self,name):
        """
            Get the level from a name
        """
        return len(name.split('0x')[1:])-1

    def get_name(self,indx):
        """
            Get the name of the cell from the global index
        """
        name = []

        level = indx[0]
        if level == 0:
            return hex(0)


        glindx = indx[1:]

        while level > 0:
            #print(glindx)
            glindx, n = self.move_index_up(glindx)
            #print(glindx,n)
            name.append(n)
            level -= 1
        return hex(0) + ''.join(name[::-1])
    def copy(self):
        """
            Return a copy of this node
        """
        return Node(name=self.name,dim=self.dim,data=self.data)
    def deepcopy(self):
        """
            Deep copy of this node
        """
        import copy
        return copy.deepcopy(self)
    def remove(self):
        """
            Remove the tree below this node
        """
        if self.child[0].data is not None:
            self.data = self.child[0].data.copy()
        self.child = [None]*self.nchildren
        self.leaf = True
    def pop(self):
        """
            Same as remove(), but also returns the tree below
        """
        new_tree =  self.deepcopy()
        self.remove()
        return new_tree
    def insert(self,name): 
        """
            Insert a new point in the tree.
            This is the same as find, but will
            grow the tree to accommodate the new point.
        """
        
        return self.find(name,insert=True)


    
    def split(self,return_children=False):
        """
            Split the node into 2^dim children, and pass the data to the
            first born.
        """
        self.leaf=False
        if self.data is None:
            mydata = None
        else:
            mydata = self.data.copy()
        self.child[0] = Node(self.name+hex(0),dim=self.dim,parent=self,data=mydata)
        for i in range(1,self.nchildren):
            self.child[i] = Node(self.name+hex(i),dim=self.dim,parent=self,data=None)
        if return_children:
            return self.child
        
    def up(self):
        """
            Move up the tree
        """
        return self.parent
    def down(self,i=0):
        """
            Move down the tree to child i
        """
        return self.child[i]
    def walk(self,printname=False,leaf_func=None,node_func=None, target_level=None,
             maxlevel=np.infty):
        """
            Recursively walk the tree.
            If the node is a leaf apply the leaf_func,
            if it is not a leaf apply the node_func.
        """
        if target_level is not None:
            if self.global_index[0] > target_level:
                return
        if self.leaf:
            if target_level is None:
                if leaf_func is not None:
                    leaf_func(self)
                if printname:
                    print(self)
            else:
                if self.global_index[0] == target_level:
                    if leaf_func is not None:
                        leaf_func(self)
            return 
        if node_func is not None:
            node_func(self)
        for c in self.child:
            c.walk(printname=printname,
                    leaf_func=leaf_func,node_func=node_func,
                   target_level=target_level,maxlevel=maxlevel)
    def depth(self):
        """
            Find the depth of the tree
        """
        res = [self.global_index[0]]
        func = lambda x: res.append(x.global_index[0])
        self.walk(leaf_func=func)
        return max(res)

    def find(self,name,insert=False):
        """
           Find the next step towards the node given
           by name.
           If insert is True then the tree will grow to 
           accomidate the new point
        """
        my_level = self.global_index[0]

        names = name.split('0x')[1:]
        target_level = len(names)-1


        if self.name == name:
            # Found it!
            #print('Node ', self.name, 'FOUND ',name)
            return self
        if my_level < target_level:
            # It's below us so we need to determine which direction to go
            new_name = '0x'+'0x'.join(names[:my_level+1])
            if self.name == new_name:
                # It's one of our descendents
                child = names[my_level+1:][0]
                #print('Node ', self.name, ' is going to child ',child)
                if self.leaf:
                    # Point doesn't exist currently  
                    if not insert:
                        #print(name,'is below this LEAF',self.name)
                        return self
                    #print('Node', self.name, 'SPLITS') 
                    self.split()
                    #print(self.child)
                #print('Node',self.name,'MOVES DOWN to ',child)
                return self.down(int(child,base=16)).find(name,insert=insert)
        # It's not below us, so move up
        #print('Node ', self.name, 'MOVES UP',name)
        if self.parent is None:
            #print('{} is not in the tree!'.format(name))
            return None
        return self.up().find(name,insert=insert)
    def find_neighbors(self):
        """
            Find the neighbors and their parents.
        """
        import itertools
        level = self.global_index[0]
        indx = self.global_index[1:]

        total_neighbors = 3**self.dim
        offsets = list(itertools.product([-1,0,1],repeat=self.dim))
        

        neighbor_indices = [(level,)+tuple([x+j for j,x in zip(i,indx)]) for i in offsets]

        # Default to None
        neighbors = [None]*total_neighbors
        upper_neighbors = [None]*total_neighbors
        for i,ind in enumerate(neighbor_indices):
            # Check that the point is inside the domain
            if all(j>=0 for j in ind): 
                n = self.get_name(ind)
                node = self.find(n)
                if node.name == n: 
                    # Node exists at this level
                    neighbors[i] = node
                    upper_neighbors[i] = node.parent
                else: 
                    # Node doesn't exist at this level, we have its parent
                    upper_neighbors[i] = node


        return offsets, neighbor_indices,neighbors, upper_neighbors

    def get_dx(self,xmin=0,xmax=1):
        """
            Return the spacing for this level
        """
        dx = 2.**(-self.global_index[0])
        return dx*(xmax-xmin)
        
        
        
    def get_coords(self,xmin=None,xmax=None,shift=False):
        """
            Return the coordinates for this node.
            xmin and xmax are the extent of the entire domain 
            and default to [0,1]
        """
        if xmin is None:
            xmin=[0]*self.dim
        if xmax is None:
            xmax=[1]*self.dim
        dx = 2.**(-self.global_index[0])
        indx = self.global_index[1:]
        shift = .5 if shift else 0
        return [(i+shift)*dx*(xo-xi) + xi for i,xi,xo in zip(indx,xmin,xmax)]
    def list_leaves(self,attr='name',func=None):
        """
            Return a list of leaves below this node.
            What is returned is controlled by attr and func
            attr returns the given attribute of the leaf
            func is a function applied to the leaf.
        """
        leaves = []
        if func is None:
            if attr is 'self' or attr is 'obj':
                func = lambda i: i
            else:
                func = lambda i: getattr(i,attr)
        
        self.walk(leaf_func=lambda i: leaves.append(func(i)))
        return leaves    
    def __repr__(self):
        return self.name
    def __str__(self):
        return self.name
