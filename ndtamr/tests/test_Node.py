from ..NDTree import Node
class TestNode():

#    def test_build(self):

    def test_copy(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            t1 = t.child[1]
            t2 = t1.copy()

            assert t1 == t2

    def test_deepcopy(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            for c in t.child:
                c.split()
            t2 = t.deepcopy()

            assert t.list_leaves() == t2.list_leaves()

    def test_depth(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            t.child[0].split()
            t.child[1].split()
            t.child[1].child[0].split()
            t.child[1].child[0].child[1].split()
            
            res = [x.depth() for x in t.list_leaves()]
            
            assert res == [4]*len(res)
        

    def test_down(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            t.child[0].split()
            t.child[1].split()
            assert t.down() == t.child[0]
            assert t.child[1].down() == t.child[1].child[0]
    def test_up(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            t.child[0].split()
            t.child[1].split()
            assert t.up() == None
            assert t.child[0].up() == t
            assert t.child[1].child[1].up() == t.child[1]
    
    def test_find(self):
        for i in range(1,4):
            t = Node(dim=i)
            t.split()
            t.child[0].split()
            t.child[1].split()
            target = t.child[1].child[0]
            assert target == t.child[0].find(target.name)
            assert target.find('0x0') == t
            
    def test_find_neighbors(self):
        t = Node(dim=2)
        t.split()
        for c in t.child:
            c.split()
        t.child[0].child[0].split()
        t.child[1].child[0].split()
        t.child[1].child[1].split()
        t.child[3].child[0].split()
        t.child[0].child[3].split()

        p1 = t.child[0].child[0].child[0]
        p2 = t.child[1].child[1].child[2]
        p3 = t.child[3].child[0].child[0]
        p4 = t.child[2].child[2]
        p1_ans = [(None, None), (None, None), (None, None), (None, None), ('0x00x00x00x0', '0x00x00x0'), ('0x00x00x00x1', '0x00x00x0'), (None, None), ('0x00x00x00x2', '0x00x00x0'), ('0x00x00x00x3', '0x00x00x0')]
        p2_ans = [('0x00x10x00x1', '0x00x10x0'), ('0x00x10x10x0', '0x00x10x1'), ('0x00x10x10x1', '0x00x10x1'), ('0x00x10x00x3', '0x00x10x0'), ('0x00x10x10x2', '0x00x10x1'), ('0x00x10x10x3', '0x00x10x1'), (None, '0x00x10x2'), (None, '0x00x10x3'), (None, '0x00x10x3')]
        p3_ans = [('0x00x00x30x3', '0x00x00x3'), (None, '0x00x10x2'), (None, '0x00x10x2'), (None, '0x00x20x1'), ('0x00x30x00x0', '0x00x30x0'), ('0x00x30x00x1', '0x00x30x0'), (None, '0x00x20x1'), ('0x00x30x00x2', '0x00x30x0'), ('0x00x30x00x3', '0x00x30x0')]
        p4_ans = [(None, None), ('0x00x20x0', '0x00x2'), ('0x00x20x1', '0x00x2'), (None, None), ('0x00x20x2', '0x00x2'), ('0x00x20x3', '0x00x2'), (None, None), (None, None), (None, None)]
        
        _,_,neighs, uneighs = p1.find_neighbors()
        neighs = [None if n is None else n.name for n in neighs]
        uneighs = [None if n is None else n.name for n in uneighs]
        
        assert [(n,nu) for n,nu in zip(neighs,uneighs)] == p1_ans
        _,_,neighs, uneighs = p2.find_neighbors()
        neighs = [None if n is None else n.name for n in neighs]
        uneighs = [None if n is None else n.name for n in uneighs]
        
        assert [(n,nu) for n,nu in zip(neighs,uneighs)] == p2_ans
        _,_,neighs, uneighs = p3.find_neighbors()
        neighs = [None if n is None else n.name for n in neighs]
        uneighs = [None if n is None else n.name for n in uneighs]
        
        assert [(n,nu) for n,nu in zip(neighs,uneighs)] == p3_ans
        _,_,neighs, uneighs = p4.find_neighbors()
        neighs = [None if n is None else n.name for n in neighs]
        uneighs = [None if n is None else n.name for n in uneighs]
        assert [(n,nu) for n,nu in zip(neighs,uneighs)] == p4_ans
    def test_index_from_bin(self):
        t = Node()
        assert (0,1) == t.index_from_bin(bin(1))
        assert (1,0,1) == t.index_from_bin(bin(5))
        assert (0,1) == t.index_from_bin('01')
        assert (1,0,1) == t.index_from_bin('101')
        assert (0,1) == t.index_from_bin(t.tobin(1))
        assert (1,0,1) == t.index_from_bin(t.tobin(5))
        
    def test_tobin(self):
        t = Node()
        assert '10' == t.tobin(2)
        assert '110' == t.tobin(6)
        
    def test_frombin(self):
        t = Node()
        assert 2 == t.frombin('10')
        assert 6 == t.frombin('110')
        assert 8 == t.frombin(bin(8))
        
    def test_get_coords(self):
        
        t = Node(dim=2)
        assert [0,0] == t.get_coords()
        assert [0.5,0.5] == t.get_coords(shift=True)
        t = Node(dim=2,xmin=(-1,-1),xmax=(1,1))
        t.split()
        assert [0,0] == t.child[3].get_coords()
        assert [0.5,0.5] == t.child[3].get_coords(shift=True)

    def test_get_global_index(self):
        t = Node(dim=2)
        
        p1 = t.get_global_index('0x00x00x00x0')
        p2 = t.get_global_index('0x00x10x10x2')
        p3 = t.get_global_index('0x00x30x00x0')
        p4 = t.get_global_index('0x00x20x2')
        assert p1 == (3, 0, 0)
        assert p2 == (3, 1, 6)
        assert p3 == (3, 4, 4)
        assert p4 == (2, 3, 0)

    def test_get_level(self):
        t = Node(dim=2)
        p0 = t.get_level('0x0')
        p1 = t.get_level('0x00x30x00x0')
        p2 = t.get_level('0x00x20x2')
        
        assert p0 == 0 
        assert p1 == 3 
        assert p2 == 2 
        
    def test_get_local_index(self):
        t = Node(dim=2)
        p0 = t.get_local_index('0x0')
        p1 = t.get_local_index('0x00x30x00x1')
        p2 = t.get_local_index('0x00x20x2')
        assert p0 == (0,0)
        assert p1 == (0,1)
        assert p2 == (1,0)
        
    def test_get_name(self):
        t = Node(dim=2)
        p1 = t.get_global_index('0x00x00x00x0')
        p2 = t.get_global_index('0x00x10x10x2')
        p3 = t.get_global_index('0x00x30x00x0')
        p4 = t.get_global_index('0x00x20x2')
        
        assert '0x0' == t.get_name((0,0,0))
        assert '0x00x00x00x0' == t.get_name((3,0,0))
        assert '0x00x10x10x2' == t.get_name((3,1,6))
        assert '0x00x20x2' == t.get_name((2,3,0))
        
        
    def test_insert(self):
        t = Node(dim=2)
        
        p1 = t.insert('0x00x20x2')
        p2 = t.insert('0x00x10x10x2')
        
        assert t.child[2].child[2].leaf 
        assert t.child[1].child[1].child[2].leaf 
        assert t.child[2].child[2] == p1 
        assert t.child[1].child[1].child[2] == p2 
        assert t.find('0x00x20x2') == p1
        assert t.find('0x00x10x10x2') == p2
       
    def test_list_leaves(self):
        t = Node(dim=2)
        t.split()
        t.child[0].split()
        l1 = ['0x00x1','0x00x2','0x00x3'] 
        l2 = ['0x00x0{}'.format(hex(i)) for i in range(4)]
        ans = l1 + l2
        
        leaves = t.list_leaves(attr='name')
        assert all([a in leaves for a in ans])
        leaves = [x.name for x in t.list_leaves(attr='self')]
        assert all([a in leaves for a in ans])
        leaves = [x.name for x in t.list_leaves(attr='self',criteria=lambda x: x.global_index[0]>1)]
        assert all([a in leaves for a in l2])
        

    def test_move_index_up(self):
        t = Node()
        assert ([2,3],'0x0') == t.move_index_up((4,6))
        assert ([32,16],'0x3') == t.move_index_up((65,33))
    
    def test_pop(self):
        t = Node(dim=2)
        t.split()
        t.child[1].split()
        p1 = t.child[1].deepcopy()
        p2 = t.child[1].pop()
        assert p1 == p2
        assert t.child[1].leaf
        
    def test_unsplit(self):
        t = Node(dim=2)
        t.split()
        t.child[1].split()
        t.child[1].unsplit()
        assert t.child[1].leaf
        assert t.child[1].child == [None]*4
    def test_split(self):
        for i in range(1,4):
            nchildren = 2**i
            t = Node(dim=i)
            t.split()
            assert [x.name for x in t.child] == ['0x0'+hex(j) for j in range(2**i)]
        
    def test_query(self):
        t = Node(dim=2)
        t.split()
        for c in t.child:
            c.split()
        n = t.query((0.6,0.3))
        assert n.name == '0x00x20x1'
        n = t.query((0.5,0.5))
        assert n.name == '0x00x30x0'
#    def test_walk(self):
        
#
#    def test_save(self):
#
#
#
#