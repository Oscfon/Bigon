r"""
Generation of trees (and forest) using the linear time algorithm from Alonso-Rémy-Schott : https://inria.hal.science/inria-00073765/document
"""

from sage.misc.prandom import randint
from sage.combinat.ordered_tree import OrderedTree

def mix(l):
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


def star(k):
    l = ['[','x']
    for _ in range(k):
        l.append('f')
    l.append(']')
    return l

def classical_pattern(l):
    patterns = []
    nb_pat = []
    nb_leaves = 1
    for k in range(len(l)):
        if l[k] != 0:
            patterns.append(star(k))
            nb_pat.append(l[k])
            nb_leaves += l[k]*(k-1)
    patterns.append(star(0))
    nb_pat.append(nb_leaves)
    return nb_pat, patterns
    
def forest_generation(nb_pat, patterns, nb_trees):
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
    
def classical_generation(l):
    (a,b)=classical_pattern(l)
    return forest_generation(a,b,1)
    





    




            