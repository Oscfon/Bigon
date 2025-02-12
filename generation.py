r"""
Generation of FatGraph starting from small examples. We ensure small operations to augment the size of the FatGraph.
"""

from surface_dynamics.misc.permutation import perm_invert
from surface_dynamics import FatGraph
from tree_generation import classical_generation

def ed_perm(k):
    if k%2==1:
        return k-1
    else:
        return k+1

def flip(vp, fp, h): #take the two permutations of the FatGraph and a semi-edge and "flip" the vertex that end the edge
    if len(vp)-1<h:
        raise ValueError("This edge doesn't belong to this graph")
    if fp[h] != h:
        h2 = ed_perm(h)
        i1 = vp[h]
        i2 = vp[i1]
        i3 = vp[i2]
        j1 = vp[h2]
        j2 = vp[j1]
        j3 = vp[j2]
    
        vp[h] = j3
        vp[j3] = i1
        vp[i2] = h
    
        vp[h2] = i3
        vp[i3] = j1
        vp[j2] = h2

        fp[h2] = i2
        fp[ed_perm(i3)]= h2
        fp[ed_perm(j1)] = i3

        fp[h] = j2
        fp[ed_perm(j3)] = h
        fp[ed_perm(i1)] = j3


def three_three(vp,fp,h):
    if len(vp)-1<h:
        raise ValueError("This edge doesn't belong to this graph")
    if fp[h] == h:
        raise ValueError("This edge form a loop")

    h2 = fp[h]
    h3 = fp[h2]
    if fp[h3] != h:
        raise ValueError("This face is not triangular")
    j = ed_perm(h)
    j2 = ed_perm(h2)
    j3 = ed_perm(h3)

    a1 = vp[j]
    a2 = vp[a1]
    b1 = vp[j2]
    b2 = vp[b1]
    c1 = vp[j3]
    c2 = vp[c1]

    vp[c2] = a1
    vp[a1] = h
    vp[j3] = c2
    vp[j] = a2
    vp[a2] = b1
    vp[b1] = h2
    vp[j2] = b2
    vp[b2] = c1
    vp[c1] = h3

    fp[ed_perm(c1)] = b2
    fp[ed_perm(c2)] = j3
    fp[j3] = c1
    fp[ed_perm(a1)] = c2
    fp[ed_perm(a2)] = j
    fp[j] = a1
    fp[ed_perm(b1)] = a2
    fp[ed_perm(b2)] = j2
    fp[j2] = b1
    

def g_n_init(g,n): #create a g genus map with n vertex
    if n < 2:
        raise ValueError("Not implemented yet, sorry !")
    if n < 2*g:
        raise ValueError("There is no such map (not enough vertex)")
    vp = [None for i in range(4*n)]
    fp = [None for i in range(4*n)]
    r = n-2*g

    for i in range(g):
        vp[4*i] = 4*(g+i)
        vp[4*i+1] = 4*(g+i)+2
        vp[4*i+2] = 4*i+1
        vp[4*i+3] = 4*i
        vp[4*(g+i)+2] = 4*i+3
        vp[4*(g+i)+1] = 4*(g+i)+3
        if i!=g-1:
            vp[4*(g+i)+3] = 4*(g+i)+5
        else: 
            vp[4*(g+i)+3] = 4*g+1
        if i!=0:
            vp[4*(g+i)] = 4*i-2
        else:
            vp[4*(g+i)] = 4*g-2

        fp[4*i] = 4*i+2
        fp[4*i+2] = 4*(g+i)+2
        fp[4*(g+i)+2] = 4*(g+i)+1
        fp[4*(g+i)+1] = 4*i

        fp[4*i+1] = 4*i+3
        fp[4*(g+i)+3] = 4*i+1
        if i!=g-1:
            fp[4*i+3] = 4*(g+i+1)
            fp[4*(g+i+1)] = 4*(g+i)+3 
        else:
            fp[4*i+3] = 4*g
            fp[4*g] = 4*(g+i)+3
    
    if g == 0:
        for j in range(r):
            vp[4*j+1] = 4*j+1
            vp[4*j] = 4*j+2
            if j!=0:
                vp[4*j+2] = 4*(j-1)
            else:
                vp[2] = 4*(r-1)
                vp[4*r-1] = 3
                fp[4*r-3] = 2
                fp[2] = 4*r-1
            
            fp[4*j] = 4*j+1
            fp[4*j+3] = 4*j
            if j != r-1:
                fp[4*j+1] = 4*j+6
                fp[4*j+6] = 4*j+3
                vp[4*j+3] = 4*j+7
            else:
                vp[4*r-1] = 3
    else:
        for j in range(r):
            vp[8*g+4*j+1] = 8*g+4*j+1
            vp[8*g+4*j] = 8*g+4*j+2
            if j!=0:
                vp[8*g+4*j+2] = 8*g+4*(j-1)
            else:
                vp[8*g+4*j+2] = 4*g-2
                fp[4*g-1] = 8*g+2
                fp[8*g+2] = 8*g-1
                vp[4*g] = 8*g+4*(r-1)
                vp[8*g-1] = 8*g+3
            
            fp[8*g+4*j] = 8*g+4*j+1
            fp[8*g+4*j+3] = 8*g+4*j
            if j != r-1:
                fp[8*g+4*j+1] = 8*g+4*j+6
                fp[8*g+4*j+6] = 8*g+4*j+3
                vp[8*g+4*j+3] = 8*g+4*j+7
            else: 
                fp[8*g+4*j+1] = 4*g
                fp[4*g] = 8*g+4*j+3
                vp[8*g+4*j+3] = 4*g+1
            
    return perm_invert(fp),perm_invert(vp)
        

r"""
Sampling of eulerian maps using Fusy-Guitter bijection (http://igm.univ-mlv.fr/~fusy/Articles/AllGenus.pdf)
"""

def incr_l(l,node):
    l[0] += 1

def incr_l_leaf(l,node):
    if not(bool(node)):
        l[0] += 1


class Decorated_blossoming_tree:

    def __init__(self, data): 
        r"""
        Return a tree and the list of the decoration on the leaves corresponding to the opening, closing and normal leaves.
        The closing leaves also have an offset.
        INPUT:
        - data -- Can be any one of the following:
            - a tuple containing the tree (as an ordered tree), the list of decorations as a DyckWord with None corresponding to unused leaves and the offset associated to each closing leaves
            - a tuple containing the tree (as an binary tree), the list of the position of the closing leaves and the offset associated to each closing leaves

        EXAMPLES::

            sage: A = Decorated_blossoming_tree((OrderedTree([[],[]]),[1,0,None],[0]))
            sage: A._tree
            [[],[]]
            sage: A._word
            [1,0,None]
            sage: A._offset
            [0]
        
            sage: B = Decorated_blossoming_tree((BinaryTree([[],[]]),[0,1,0],[0,0,1,0]))
            sage: B._tree
            [[], [[], [], []], [[], [], []]]
            sage: B._word
            [1,1,0,0,1,1,0,0]
            sage: B._offset
            [0,0,1,0]
        """

        if isinstance(data, tuple) and len(data) == 3:
            if isinstance(data[0], (OrderedTree)) and isinstance(data[1], (list, DyckWord)) and isinstance(data[2], list):
                self._tree = data[0]
                self._word = list(data[1])
                self._offset = data[2]

            elif isinstance(data[0], BinaryTree) and isinstance(data[1], list) and isinstance(data[2], list):
                self._offset = data[2]
                w = [1]
                current = data[0]
                sub_tree = []
                corner = [0]
                left_tree = []
                node = 0
                out = []
                while corner[-1] != 2 or current != data[0]:
                    if current.is_empty():
                        corner.pop()
                        corner[-1] += 1
                        current = sub_tree.pop()
                        w.append(0)
                        left_tree.append(OrderedTree([]))
                    elif corner[-1] == 0:
                        out.append(data[1][node])
                        node += 1
                        if out[-1] == 0:
                            w.append(1)
                            left_tree.append(OrderedTree([]))
                        sub_tree.append(current)
                        current = current[0]
                        corner.append(0)
                    elif corner[-1] == 1:
                        if out[-1] == 1:
                            w.append(1)
                            left_tree.append(OrderedTree([]))
                        sub_tree.append(current)
                        current = current[1]
                        corner.append(0)
                    else:
                        if out[-1] == 2:
                            w.append(1)
                            left_tree.append(OrderedTree([]))
                        t3 = left_tree.pop()
                        t2 = left_tree.pop()
                        t1 = left_tree.pop()
                        left_tree.append(OrderedTree([t1,t2,t3]))
                        out.pop()
                        corner.pop()
                        corner[-1] += 1
                        current = sub_tree.pop()
                if out[-1] == 2:
                    w.append(1)
                    left_tree.append(OrderedTree([]))
                t3 = left_tree.pop()
                t2 = left_tree.pop()
                t1 = left_tree.pop()
                self._tree = OrderedTree([t1,t2,t3])
                self._word = w
        else:
            raise ValueError("The data has not a good type.")


    def _check(self):
        closing_leaves = 0
        opening_leaves = 0
        for elt in self._word:
            if elt == 0:
                closing_leaves += 1
            elif elt == 1:
                opening_leaves += 1
        l = [0]
        self._tree.iterative_post_order_traversal(lambda node: incr_l_leaf(l,node))
        leaves = l[0]+1
        if leaves != len(self._word):
            raise ValueError("Not the same value of leaves and decorations on the leaves.")
        elif opening_leaves != closing_leaves:
            raise ValueError("Not the same value of opening and closing leaves.")
        elif opening_leaves != len(self._offset):
            raise ValueError("Not the same value of closing leaves and offset on this leaves.")

        height = 0
        index = 0
        for elt in self._word:
            if elt==0:
                if height<=0:
                    raise ValueError("The decoration on the leaves is not a Dyck word.")
                elif self._offset[index]<0 or self._offset[index]>=height:
                    raise ValueError("The offset on the {} closing leaf is bigger of the height or smaller than 0.".format(index))
                index += 1
                height -= 1
            elif elt==1:
                height += 1


    def to_Fat_graph(self):
        r"""
        Convert the blossoming to FatGraph using http://igm.univ-mlv.fr/~fusy/Articles/AllGenus.pdf

        EXAMPLES::

            sage: A = Decorated_blossoming_tree((BinaryTree([[],[]]),[1,1,1],[0,0,0,0]))
            sage: A._check()
            sage: G = A.to_Fat_graph()
            sage: G
            FatGraph('(0,2,6,8)(1,4,5,3)(7,10,11,9)', '(0,3)(1,8,11,7,2,5)(4)(6,9)(10)')
            sage: G.genus()
            0

            sage: B = Decorated_blossoming_tree((BinaryTree([[],[]]),[0,0,1],[1,1,0,0]))
            sage: B._check()
            sage: G2 = B.to_Fat_graph()
            sage: G2
            FatGraph('(0,2,4,8)(1,5,6,3)(7,10,11,9)', '(0,3)(1,8,11,7,5,2,6,9,4)(10)')
            sage: G2.genus()
            1
        """
        l = [0]
        self._tree.iterative_post_order_traversal(lambda node: incr_l(l,node))
        nodes = l[0]
        m = nodes-len(self._offset)
        vp = [None for i in range(2*m)]
        opening_edges = [1]
        free_index = 1
        last = [0]
        first = [0]
        current = self._tree
        word_index = 1
        offset_index = 0
        sub_tree = []
        corner = [0]
        while current != self._tree or corner[-1] != len(self._tree):
            if not(bool(current)):
                if self._word[word_index] == 0:
                    ope_edg = opening_edges[-1-self._offset[offset_index]]
                    vp[last[-1]] = ope_edg
                    opening_edges.remove(ope_edg)
                    last[-1] = ope_edg
                    offset_index += 1
                elif self._word[word_index] == 1:
                    vp[last[-1]] = 2*free_index
                    last[-1] = 2*free_index
                    opening_edges.append(2*free_index+1)
                    free_index += 1
                else:
                    vp[last[-1]] = 2*free_index
                    last[-1] = 2*free_index
                    vp[2*free_index+1] = 2*free_index+1
                    free_index += 1
                word_index += 1
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
                
            elif corner[-1] == 0:
                first.append(last[-1])
                if bool(current[0]):
                    vp[last[-1]] = 2*free_index
                    last[-1] = 2*free_index
                    last.append(2*free_index + 1)
                    free_index += 1
                sub_tree.append(current)
                current = current[corner[-1]]
                corner.append(0)
            
            elif corner[-1] < len(current):
                if bool(current[corner[-1]]):
                    vp[last[-1]] = 2*free_index
                    last[-1] = 2*free_index
                    last.append(2*free_index + 1)
                    free_index += 1
                sub_tree.append(current)
                current = current[corner[-1]]
                corner.append(0)
            else:
                vp[last.pop()] = first.pop()
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
        vp[last.pop()] = first.pop()
        return FatGraph(vp=vp)

    def random_k_val(n,k,planar=False):
        r"""
        Sample randomly a blossoming tree with n nodes of degree 2k
        """
        l = []
        for i in range(2*k-1):
            l.append(0)
        l.append(n)
        t = classical_generation(l)[0]
        m = [0]
        t.iterative_post_order_traversal(lambda node: incr_l_leaf(m,node))
        leaves = m[0]+1
        opening_leaves = leaves//2
        w = DyckWords(opening_leaves).random_element()
        if planar:
            off = [0 for _ in range(leaves//2)]
        else:
            h = 0
            off = []
            for elt in w:
                if elt==1:
                    h += 1
                elif elt==0:
                    h -= 1
                    off.append(randint(0,h))
        res = Decorated_blossoming_tree((t,w,off))
        res._check()
        return res
        
        

    def random_4_val_reject(n,planar=False):
        r"""
        Sample randomly an eulerian tree with n nodes of degree 4.
        This work by reject.
        """
        t = BinaryTrees(n).random_element()
        b = [randint(0,1)]+[randint(0,2) for _ in range(n-1)]
        off = [0 for _ in range(n+1)]
        res = None
        work = True
        while work:
            try: 
                res = Decorated_blossoming_tree((t,b,off))
                res._check()
                work = False
            except ValueError:
                b = [randint(0,1)]+[randint(0,2) for _ in range(n-1)]
        h = 0
        if not planar:
            off = []
            for elt in res._word:
                if elt == 1:
                    h += 1
                else:
                    h -= 1
                    off.append(randint(0,h))
        res = Decorated_blossoming_tree((t,b,off))
        res._check()
        return res
        

    def random_element(n, zero_vertex=False, planar=False):
        r"""
        Sample randomly a tree, a Dyck word on the leaves of this tree (including the root) and an offset for the closing leaves.
        """

        t = OrderedTrees(n).random_element()
        l = [0]
        t.iterative_post_order_traversal(lambda node: incr_l_leaf(l,node))
        leaves = l[0]+1
        nb_zero = 0
        if zero_vertex:
            nb_zero = randint(0,leaves//2-1)
        opening_leaves = leaves//2-nb_zero
        nb_zero += leaves%2
        w = DyckWords(opening_leaves).random_element()
        word = list(w)
        for _ in range(nb_zero):
            index = randint(1,len(word))
            word = word[:index]+[None]+word[index:]
        if planar:
            off = [0 for _ in range(leaves//2)]
        else:
            h = 0
            off = []
            for elt in word:
                if elt==1:
                    h += 1
                elif elt==0:
                    h -= 1
                    off.append(randint(0,h))
        res = Decorated_blossoming_tree((t,word,off))
        res._check()
        return res







    
        