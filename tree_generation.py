r"""
Generation of trees (and forest) using the linear time algorithm from Alonso-Rémy-Schott : https://inria.hal.science/inria-00073765/document

The method takes a list of patterns and the number of occurences of all the patterns and generates a forest which respects the pattern constraint.
The pattern are coded using where 'x' denotes a vertex, 'y' denotes a classical edge, 'f' a semi-edge and 'o' a multi-edge.
The order of coding is defined as follow : if v is a vertex with children v_1, ..., v_p and edges a_1, ..., a_q then the coding of v will be :
t(v_1)t(v_2)...t(v_p)x t(a_p)...t(a_1)
"""

from sage.misc.prandom import randint
from sage.combinat.ordered_tree import OrderedTree

def mix(l):
    r"""
    Take a list of occurences and create uniformly at random a list res of integer where i appears exactly l[i] times.
    """
    l_comp = []
    for i in range(len(l)):
        for j in range(l[i]):
            l_comp.append(i)
    res = [0 for i in l_comp]
    for i in range(len(l_comp)-1,-1,-1):
        N = randint(0,i)
        res[i] = l_comp[N]
        l_comp[N] = l_comp[i]
    return res

def replace_patern(mix_list,pattern_words):
    r"""
    Take a list of index and the corresponding patterns and create a word obtained by the concatenation of this pattern.

    EXAMPLES::

        sage: mix_list = [1, 0, 1, 1]
        sage: pattern_word = ["yx", "ox"]
        sage: replace_patern(mix_list,pattern_word)
        ['(', ')', 'x', 'y', 'x', '(', ')', 'x', '(', ')', 'x']
        
    """
    res = []
    for elt in mix_list:
        for letter in pattern_words[elt]:
            if letter == 'o':
                res.append('(')
                res.append(')')
            else:
                res.append(letter)
    return res


def missing_edges(l,pattern_words,p):
    r"""
    Take the list of patterns, the number of occurences of this patterns and the number of trees in the forest and compute the number of edges that must be added to have a correct forest.
    """
    n = 0
    e = 0
    c = 0
    d = 0
    for i in range(len(l)):
        n_m = 0
        e_m = 0
        c_m = 0
        d_m = 0
        for letter in pattern_words[i]:
            if letter == 'x':
                n_m += 1
            elif letter == 'y':
                e_m += 1
            elif letter == 'f':
                c_m += 1
            elif letter == 'o':
                d_m += 1
        n += n_m*l[i]
        e += e_m*l[i]
        c += c_m*l[i]
        d += d_m*l[i]
    if n-p<e+c:
        raise ValueError("There is no forest with this constraints, not enough vertex or too much trees.")
    elif n-p==e+c:
        return 0
    elif d == 0:
        raise ValueError("There is no forest with this constraints, too much vertices and not enough edges.")
    else:
        return n-p-e-c


def adding_edges(nb_edges,w,pos_d):
    r"""
    Take a word w that encodes a tree and pos_d the list of positions of multi-edges in w and add nb_edges edge to w.
    """
    number_edges = nb_edges
    pos = 0
    d = len(pos_d)
    diff = 0
    while number_edges > 0:
        r = randint(1,d-pos+number_edges-1)
        if r > number_edges:
            pos += 1
        else:
            w = w[:pos_d[pos]+diff]+['f']+w[pos_d[pos]+diff:]
            diff += 1
            number_edges -= 1
    return w


def cyclic_permutation(w,p):
    r"""
    Take a word w that encodes the forest and p the number of tree in the forest and perform a cyclic permutation on w such that the word begins by a well formed tree.
    """
    pos = 0
    height = 0
    w_min = 0
    pos_min = 0
    for elt in w:
        if elt == 'x':
            if height <= w_min:
                pos_min = pos
                w_min = height
            height += 1
        elif elt == 'f' or elt == 'y':
            height -= 1
        pos += 1

    height_chosen = w_min + randint(0,p-1)
    pos = 0
    height = 0
    for elt in w:
        if elt == 'x':
            if height == height_chosen:
                beginning = pos
            height += 1
        elif elt == 'f' or elt == 'y':
            height -= 1
        pos+=1

    beginning -= 1
    res = []
    for i in range(len(w)):
        index = (beginning+i)%len(w)
        res.append(w[index])
    return res

def to_forest(w):
    r"""
    Convert a word into a forest seen as a list of trees.
    """
    globalstack = []
    patternstack = []
    pos = 0
    while pos < len(w):
        if w[pos] == '[':
            pos += 1
        vertex = []
        while w[pos] != 'x' and w[pos] != ']':
            value = w[pos]
            if value == 'f':
                vertex.append(globalstack.pop())
            elif value == 'y':
                vertex.append(patternstack.pop())
            pos += 1
        if w[pos] == 'x':
            pos += 1
        if w[pos] == ']':
            pos += 1
            globalstack.append(OrderedTree(vertex))
        else:
            patternstack.append(OrderedTree(vertex))
    
    res = []
    for elt in globalstack:
        res.append(OrderedTree(elt))
    return res

def forest_generation(nb_pat, patterns, nb_trees):
    r"""
    Generate a forest of nb_trees trees where the pattern patterns[i] appears exactly nb_pat[i] times in the forest. 
    """
    l = mix(nb_pat)
    w = replace_patern(l,patterns)
    me = missing_edges(nb_pat,patterns,nb_trees)
    if me>0:
        pos_d = []
        for i in range(len(w)):
            if w[i] == '(':
                pos_d.append(i+1)
        w = adding_edges(me,w,pos_d)
    w = cyclic_permutation(w,nb_trees)
    return to_forest(w)


def star(k):
    r"""
    Return the pattern corresponding to a star with k semi-edges.

    EXAMPLES::

        sage: star(3)
        ['[', 'x', 'f', 'f', 'f', ']']
    """
    l = ['[','x']
    for _ in range(k):
        l.append('f')
    l.append(']')
    return l

def classical_pattern(l):
    r"""
    Compute the list of patterns and the list of occurences corresponding to a tree with l[i] nodes of degrees i+1.

    EXAMPLES::

        sage: classical_pattern([0,12])
        ([12, 13], [['[', 'x', 'f', 'f', ']'], ['[', 'x', ']']])
        sage: classical_pattern([1,1])
        ([1, 1, 2], [['[', 'x', 'f', ']'], ['[', 'x', 'f', 'f', ']'], ['[', 'x', ']']])
        sage: 
    """
    patterns = []
    nb_pat = []
    nb_leaves = 1
    for k in range(len(l)):
        if l[k] != 0:
            patterns.append(star(k+1))
            nb_pat.append(l[k])
            nb_leaves += l[k]*(k)
    patterns.append(star(0))
    nb_pat.append(nb_leaves)
    return nb_pat, patterns
    
def classical_generation(l):
    r"""
    Generate uniformly at random a tree with l[i] nodes of degrees i+1.
    """
    (a,b)=classical_pattern(l)
    return forest_generation(a,b,1)
    





    




            