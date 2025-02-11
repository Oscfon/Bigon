r"""
    Generation of FatGraph starting from small examples. We ensure small operations to augment the size of the FatGraph.
    All of this FatGraph are 4-valency FatGraph.
"""

from surface_dynamics.misc.permutation import perm_invert
from surface_dynamics import FatGraph



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



    Note de programmation : Tirer un arbre, un mot de Dick -> blossoming tree + une fonction de matching sur chaque arête entrante !
    
    Le degré des nœuds de l'arbre est le degré des sommets de la carte !
    Il faut un arbre eulérien ie : tous les sommets ont degré pair (2k) avec k-1 feuilles
"""

class Eulerian_tree:
    """r
        Class of eulerian tree with the implemented height function 
    """
    
    def __init__(self,t,b): # t is a binary tree and b is the corner in which is the out-leaf on each node.
        current = t
        left_tree = []
        sub_tree = []
        corner = [0]
        h = 1
        out = []
        while corner[-1] != 2 or current != t:
            if current.is_empty():
                left_tree.append(['c',h])
                corner.pop()
                corner[-1] += 1
                current = sub_tree[-1]
                sub_tree.pop()
                h -= 1
                if h<0:
                    raise ValueError("The tree is not balanced")
            elif corner[-1] == 0:
                out.append(b.pop())
                if out[-1] == 0:
                    h += 1
                sub_tree.append(current)
                current = current[0]
                corner.append(0)
            elif corner[-1] == 1:
                if out[-1] == 1:
                    h += 1
                sub_tree.append(current)
                current = current[1]
                corner.append(0)
            else:
                t2 = left_tree.pop()
                t1 = left_tree.pop()
                if out[-1] == 0:
                    left_tree.append(['o',t1,t2])
                elif out[-1] == 1:
                    left_tree.append([t1,'o',t2])
                else:
                    left_tree.append([t1,t2,'o'])
                    h += 1
                out.pop()
                corner.pop()
                corner[-1] += 1
                current = sub_tree[-1]
                sub_tree.pop()
        t2 = left_tree.pop()
        t1 = left_tree.pop()
        if out[-1] == 0:
            self._tree = ['o', 'o', t1, t2]
        elif out[-1] == 1:
            self._tree = ['o', t1,'o', t2]
        else:
            self._tree = ['o',t1, t2, 'o']

    def random_4_val(n):
        t = BinaryTrees(n).random_element()
        b = [randint(0,2) for _ in range(n)]
        res = None
        while res == None:
            try: 
                res = eulerian_tree(t,b)
            except ValueError:
                b = [randint(0,2) for _ in range(n)]
        return res

def nb_leaves(t):
    if not(bool(t)):
        return 1
    else:
        res = 0
        for elt in t:
            res += nb_leaves(elt)
        return res


class Decored_blossoming_tree:
    """r
        Class of tree with a decoration to each point corresponding to the off-set 
    """

    def __init__(self, t, w, offset = None): 
        # t is a tree and w a Dyck word with length the number of leaves plus one of t. The offset is the list of decorations. If offset is None, the decorations are random
        specifize = False
        if offset != None: 
            specifize = True
            off_index = 0
            if len(offset) != len(w)//2:
                raise ValueError("Not the same lenght of offset than the number of closing leaves.")
        current = t
        left_tree = []
        sub_tree = []
        corner = [0]
        w_index = 1
        h = 1
        size = 0
        while current != t or corner[-1] != len(t):
            if not bool(current):
                if w_index == len(w):
                    raise ValueError("The Dyck word is too short.")
                elif w[w_index]==0:
                    dec = randint(0,h-1)
                    if specifize:
                        if offset[off_index]>h-1:
                            raise ValueError("Too big offset.")
                        dec = offset[off_index]
                        off_index += 1
                    left_tree.append(['c',dec])
                    h -= 1
                    w_index+=1
                    size += 1
                else:
                    left_tree.append(['o'])
                    h += 1
                    w_index+=1
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
                
            elif corner[-1] < len(current):
                sub_tree.append(current)
                current = current[corner[-1]]
                corner.append(0)
            else:
                l = []
                for i in range(len(current)):
                    l.append(left_tree.pop())
                l.reverse()
                left_tree.append(l)
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
                size += 1
        if w_index < len(w):
            raise ValueError("The Dyck word is too long")
        self._tree = [['o']]+left_tree
        self._m = size

    def to_Fat_graph(self):
        r"""
            Convert the blossoming to FatGraph using http://igm.univ-mlv.fr/~fusy/Articles/AllGenus.pdf
        """
        vp = [None for i in range(2*self._m)]
        opening_edges = []
        free_index = 0
        last = [0]
        first = [0]
        current = self._tree
        sub_tree = []
        corner = [0]
        while current != self._tree or corner[-1] != len(self._tree):
            if current[0] == 'c':
                dec = current[1]
                ope_edg = opening_edges[-dec-1]
                vp[last[-1]] = ope_edg
                opening_edges.remove(ope_edg)
                last[-1] = ope_edg
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
            
            elif current[0] == 'o':
                vp[last[-1]] = 2*free_index
                last[-1] = 2*free_index
                opening_edges.append(2*free_index+1)
                free_index += 1
                corner.pop()
                corner[-1] += 1
                current = sub_tree.pop()
                
            elif corner[-1] == 0:
                first.append(last[-1])
                if current[0][0] != 'o' and current[0][0] != 'c':
                    vp[last[-1]] = 2*free_index
                    last[-1] = 2*free_index
                    last.append(2*free_index + 1)
                    free_index += 1
                sub_tree.append(current)
                current = current[corner[-1]]
                corner.append(0)
            
            elif corner[-1] < len(current):
                if current[corner[-1]][0] != 'o' and current[corner[-1]][0] != 'c':
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












    
        