#import numpy as np
#import h5py
#from .DataClass import Sim as Data

class Node():
    """
        Nodes either point to their children or they have no
        children and instead hold some data.
    """
    def __init__(self,name,dim=2,file=None,parent=None,data=None):
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
        #self.name = '/'.join(self.indx)
        #self.level = len(indx)-1
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
    def split(self):
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
        return
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
    def walk(self,printname=False,leaf_func=None,node_func=None):
        """
            Recursively walk the tree.
            If the node is a leaf apply the leaf_func,
            if it is not a leaf apply the node_func.
        """
        if self.leaf:
            if leaf_func is not None:
                leaf_func(self)
            if printname:
                print(self)
            return self
        if node_func is not None:
            node_func(self)
        for c in self.child:
            c.walk(printname=printname,
                    leaf_func=leaf_func,node_func=node_func)
    def depth(self):
        """
            Find the depth of the tree
        """
        res = [self.global_index[0]]
        func = lambda x: res.append(x.global_index[0])
        self.walk(leaf_func=func)
        return max(res)

    def find(self,name):
        """
           Find the next step towards the node given
           by name.
        """
        my_level = self.global_index[0]

        names = name.split('0x')[1:]
        target_level = len(names)-1


        if self.name == name:
            # Found it!
            print('Node ', self.name, ' found ',name)
            return self
        if my_level < target_level:
            # It's below us so we need to determine which direction to go
            new_name = '0x'+'0x'.join(names[:my_level+1])
            if self.name == new_name:
                # It's one of our descendents
                child = names[my_level+1:][0]
                print('Node ', self.name, ' is going to child ',child)
                if self.leaf:
                    return self
                else:
                    return self.down(int(child,base=16)).find(name)
        # It's not below us, so move up
        print('Node ', self.name, ' is going up to find ',name)
        if self.parent is None:
            print('{} is not in the tree!'.format(name))
            return None
        return self.up().find(name)
    def find_neighbors(self):
        """
            Find the neighbors and their parents.
        """
        import itertools
        level = self.global_index[0]
        indx = self.global_index[1:]

        offsets = list(itertools.product([-1,0,1],repeat=self.dim))

        neighbor_indices = [(level,)+tuple([x+j for j,x in zip(i,indx)]) for i in offsets]

        neighbors = []
        upper_neighbors = []
        for ind in neighbor_indices:
            if not all(i>=0 for i in ind):
                neighbors.append(None)
                upper_neighbors.append(None)
            else:
                n = self.get_name(ind)
                node = self.find(n)
                if node.name == n: # Node exists
                    neighbors.append(node)
                    upper_neighbors.append(node.parent)
                else: # Node doesn't exist, we have its parent
                    neighbors.append(None)
                    upper_neighbors.append(node)


        return offsets, neighbor_indices,neighbors, upper_neighbors


    def __repr__(self):
        return self.name
    def __str__(self):
        return self.name
