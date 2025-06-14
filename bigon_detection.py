from flatsurf import *
from surface_dynamics import *
from collections import deque


def labels(G):
    r"""
    Assign a label (an integer) to each vertex and each face of G.
    Compute for each darts the vertex and the face it belongs. 
    """
    vertices = G.vertices()
    faces = G.faces()
    darts = list(G.darts())
    darts_to_vertex = [None for _ in darts]
    darts_to_face = [None for _ in darts]
    for i in range(len(vertices)):
        for e in vertices[i]:
            darts_to_vertex[e]=i
    for i in range(len(faces)):
        for e in faces[i]:
            darts_to_face[e]=i
    return darts_to_vertex, darts_to_face


def tree_co_tree(G):
    r"""
    Compute a tree/co-tree decomposition of G.
    Return a list res of length the number of darts of G.
    res[e] is 0 if the dart e is in the tree, 1 if e is in the co-tree and 2 otherwise.
    """
    darts_to_vertex, darts_to_face = labels(G)
    vp = G.vertex_permutation(copy=False)
    T_found=[False for _ in G.vertices()]
    T_found[0] = True
    res = [2 for _ in G.darts()]
    e = 0
    path = [e]

    if not T_found[darts_to_vertex[1]]:
        T_found[darts_to_vertex[1]]=True
        res[0]=0
        res[1]=0
        path.append(1)
        e=vp[1]
    else:
        e=vp[0]
    
    while len(path)>0:
        e1 = e+1 if e%2==0 else e-1
        if e == path[-1]:
            path.pop()
            e=vp[e1]
        elif not T_found[darts_to_vertex[e1]]:
            T_found[darts_to_vertex[e1]]=True
            res[e]=0
            res[e1]=0
            path.append(e1)
            e=vp[e1]
        else:
            e=vp[e]

    fp=G.face_permutation(copy=False)
    C_found=[False for _ in G.faces()]
    C_found[0] = True
    e = fp[0]
    path = [0]
    
    while len(path)>0:
        e1 = e+1 if e%2==0 else e-1
        if e == path[-1]:
            path.pop()
            e=fp[e1]
        elif res[e]==0:
            e=fp[e]
        elif not C_found[darts_to_face[e1]]:
            C_found[darts_to_face[e1]]=True
            res[e]=1
            res[e1]=1
            path.append(e1)
            e=fp[e1]
        else:
            e=fp[e]
        
    return res


def tree_contraction(G, treecotree):
    r"""
    Compute the Fatgraph F obtained from G by:
        - contracting all the edge in the tree of treecotree
        - removing all the edge in the co-tree of treecotree
    Compute a list correspondence of length the number of darts of G such that correspondence[e] is:
        - None if e is an edge of the tree of treecotree
        - d such that fp[d]=e in G after contracting the edge of the tree if e is an edge of the cotree
        - the corresponding edge of F otherwise
    Compute a list recor such that recor[f] for f a dart in F is the corresponding dart in G
    Compute a list rank such that rank[e] is:
        - None if e is an edge of the tree of treecotree
        - the index of e around the vertex if e is an edge of the cotree
        - the total number of edge of the cotree between e and the next edge in F

    INPUT:

    treecotree should be a list of lenght the number of darts in G where G[e] is:
        - 0 if e is in the tree
        - 1 if e is in the co-tree
        - 2 otherwise
    """
    fp=G.face_permutation(copy=False)
    correspondence = [None for _ in G.darts()]
    rank = [None for _ in G.darts()]
    nfp=[None for _ in range(4*G.genus())]
    recor = [None for _ in range(4*G.genus())]
    i=0
    while treecotree[i]!=2:
        i+=2
    last=0
    correspondence[i]=0
    correspondence[i+1]=1
    recor[0]=i
    recor[1]=i+1
    e=fp[i]
    free_index=2
    r=0
    
    while e!=i:
        e1=e+1 if e%2==0 else e-1
        if treecotree[e]==1:
            correspondence[e]=last
            rank[e]=r
            r+=1
            e=fp[e1]
            
        elif treecotree[e]==0:
            e=fp[e]
        else:
            if correspondence[e]==None:
                correspondence[e]=free_index
                correspondence[e1]=free_index+1
                recor[free_index]=e
                recor[free_index+1]=e1
                free_index+=2
            rank[e]=r
            r=0
            nfp[last]=correspondence[e]
            last=correspondence[e]
            e=fp[e]
    rank[i]=r
    nfp[last]=0
    F=FatGraph(fp=nfp)
    return F, correspondence, recor, rank
            
def quad_system(G, treecotree):
    r"""
    Compute the quad system associated to G and the tree/co-tree decomposition of treecotree.

    
    INPUT:

    treecotree should be a list of lenght the number of darts in G where G[e] is:
        - 0 if e is in the tree
        - 1 if e is in the co-tree
        - 2 otherwise


    OUTPUT:

    - Q the FatGraph encoding the quad system of G
    - cor2 a list such that cor2[e] for a e a dart of G is the list of darts of Q corresponding to e
    """
    
    (F,cor, _, _)=tree_contraction(G,treecotree)
    vp=F.vertex_permutation(copy=False)
    fp=F.face_permutation(copy=False)
    qvp=[None for _ in range(8*G.genus())]
    remember=[None for _ in range(4*G.genus())]
    remember2=[None for _ in range(4*G.genus())]
    e=fp[0]
    last_edge=0
    next_edge=1
    remember[0]=vp[0]
    remember2[0]=0
    while e!=0:
        qvp[2*last_edge]=2*next_edge
        remember[next_edge]=vp[e]
        remember2[e]=next_edge
        e=fp[e]
        last_edge=next_edge
        next_edge+=1
    qvp[-2]=0
    for ne in range(4*G.genus()):
        qvp[2*ne+1]=2*remember2[remember[ne]]+1
    Q=FatGraph(vp=qvp)
    qfp=Q.face_permutation(copy=False)
    cor2 = []
    for e in G.darts():
        if treecotree[e]==0:
            cor2.append([])
        elif treecotree[e]==1:
            start=fp[cor[e]]
            e1=e+1 if e%2==0 else e-1
            end=fp[cor[e1]]
            das=2*remember2[start]+1
            dae=2*remember2[end]
            cor2.append([das,dae])
        else:
            edge=cor[e]
            da=2*remember2[edge]+1
            cor2.append([da,qvp[da-1]])
    return Q,cor2


def turn(Q,e,f):
    r"""
    Compute the number of turn from dart e to dart f.
    """
    vp=Q.vertex_permutation(copy=False)
    e1=e+1 if e%2==0 else e-1
    current=e1
    if f==e1:
        return 0
    else:
        current=vp[e1]
    turn=1
    while current!=e1 and f!=current:
        turn+=1
        current=vp[current]
    if current==e1:
        raise ValueError("The edge {} does not follow the edge {}.".format(f,e))
    return turn

def turn_add(t,turn,num):
    if num<=0:
        return
    if len(t)==0:
        t.append((turn,num))
    else:
        if t[-1][0]==turn:
            (a,b)=t.pop()
            t.append((a,b+num))
        else:
            t.append((turn,num))

def turn_add_left(t,turn,num):
    if num<=0:
        return
    if len(t)==0:
        t.append((turn,num))
    else:
        if t[0][0]==turn:
            (a,b)=t.popleft()
            t.appendleft((a,b+num))
        else:
            t.appendleft((turn,num))

def turn_modif(t,x,d):
    if len(t)==0:
        return
    (t1,n1)=t.pop()
    r=(t1+x)%d
    if n1>1:
        t.append((t1,n1-1))
        t.append((r,1))
    else:
        turn_add(t,r,1)

def turn_modif_left(t,x,d):
    if len(t)==0:
        return
    (t1,n1)=t.popleft()
    r=(t1+x)%d
    if n1>1:
        t.appendleft((t1,n1-1))
        t.appendleft((r,1))
    else:
        turn_add_left(t,r,1)

def simplification(Q,d,c,t,e):
    r"""
    Update the geodesic representant c and the turn sequence t after adding the edge e to the path.

    INPUT:

    - Q -- underlying quad system 
    - d -- the degree of Q (we assume all vertices of Q have the same degree d)
    - c -- a geodesic path in Q (as a deque)
    - t -- the turn sequence corresponding to c (as a deque)
    - e -- an edge following c

    OUTPUT

    Return a boolean saying if we flatten a bracket or remove a spur.
    This function update the value of c and t.
    """
    fp=Q.face_permutation(copy=False)
    if len(c)==0:
        c.append(e)
        return False
    elif len(c)==1:
        newturn = turn(Q,c[-1],e)
        if newturn==0:
            c.pop()
            return True
        else:
            t.append((newturn,1))
            c.append(e)
            return False
    if len(t)==0:
        raise ValueError("The length of t for c = {} shouldn't be 0.".format(c))
    last=c[-1]
    (last_t, n_last)=t.pop()
    newturn=turn(Q,last,e)
    if newturn==0:
        c.pop()
        if n_last!=1:
            t.append((last_t,n_last-1))
        return True
    elif newturn==1 and last_t==1:
        bracket_removal(Q, c, t, True, 0, d)
        return True
    elif newturn==d-1 and last_t==d-1:
        bracket_removal(Q, c, t, False, 0, d)
        return True
    elif len(t)>=1 and newturn==1 and last_t==2 and t[-1][0]==1:
        bracket_removal(Q, c, t, True, n_last, d)
        return True
    elif len(t)>=1 and newturn==d-1 and last_t==d-2 and t[-1][0]==d-1:
        bracket_removal(Q, c, t, False, n_last, d)
        return True
    elif newturn==last_t:
        c.append(e)
        t.append((newturn,n_last+1))
        return False
    else:
        c.append(e)
        t.append((last_t,n_last))
        t.append((newturn,1))
        return False


def spur_removal(c,t):
    c.pop()
    c=c.popleft()
    (t1,n1)=t.pop()
    if n1>1:
        t.append(t1,n1-1)
    (t2,n2)=t[0]
    if n2>1:
        t[0]=(t2,n2)
    else:
        t=t.popleft()

def bracket_removal(Q, c, t, positive, length, d):
    fp=Q.face_permutation(copy=False)
    if positive:
        l=deque([])
        for i in range(length+1,0,-1):
            edge=c[-i]
            edge1=edge+1 if edge%2==0 else edge-1
            l.append(fp[fp[edge1]])
        for i in range(length+2):
            c.pop()
        c=c.extend(l)
        if length!=0:
            (t2,n2)=t.pop()
            if n2!=1:
                raise ValueError("The path wasn't geodesic")
        turn_modif(t,-1,d)
        turn_add(t,d-2,length)
    else:
        l=deque([])
        for i in range(length+1,0,-1):
            edge=fp[fp[c[-i]]]
            edge1=edge+1 if edge%2==0 else edge-1
            l.append(edge1)
        for i in range(length+2):
            c.pop()
        c=c.extend(l)
        if length!=0:
            (t2,n2)=t.pop()
            if n2!=1:
                raise ValueError("The path wasn't geodesic")
        turn_modif(t,1,d)
        turn_add(t,2,length)

def origin_simplification(Q, d, c, t):
    r"""
    Perform the geodesic simplification at the based point of c
    """
    if len(c)==0:
        return
    elif len(c)==1:
        raise ValueError("The path wasn't geodesic")
    first_turn = turn(Q, c[-1], c[0])
    if first_turn == 0: #case 1
        spur_removal(c,t)
        origin_simplification(Q, d, c, t)
    elif first_turn==1 and len(t)==1 and t[0][0]==2: #case 2.1
        c.pop()
        c=c.popleft()
        l=deque([])
        for e in c:
            e1=e+1 if e%2==0 else e-1
            l.append(fp[fp[e1]])
        c.clear()
        c.extend(l)
        t=deque([(d-2, t[0][1]-2)])
    elif first_turn==d-1 and len(t)==1 and t[0][0]==d-2: #case 2.2
        c.pop()
        c=c.popleft()
        l=deque([])
        for e in c:
            e1=e+1 if e%2==0 else e-1
            l.append(fp[fp[e1]])
        c.clear()
        c.extend(l)
        t=deque([(2, t[0][1]-2)])
    elif first_turn==1 or first_turn==d-1: #case 4
        a = c.popleft()
        (t1,n1) = t.popleft()
        if n1 != 1:
            t.appendleft((t1,n1-1))
        bremoved=simplification(Q,d,c,t,a)
        if bremoved:
            first_turn=turn(Q, c[-1], c[0])
            while first_turn == 0:
                spur_removal(c,t)
                first_turn=turn(Q, c[-1], c[0])
        else:
            new_turn=turn(Q, c[-1], c[0])
            if new_turn+first_turn==d:
                return
            else:
                origin_simplification(Q,d,c,t)
    elif first_turn==2 and len(t)>=3: #case 3
        while first_turn==2:
            c.append(c.popleft())
            (t1,n1)=t.popleft()
            if n1!=1:
                t.appendleft((t1,n1-1))
            (t2,n2)=t.pop()
            if t2==2:
                t.append((2,n2+1))
            else:
                t.append((t2,n2))
                t.append((2,1))
            first_turn=t1
        if first_turn==1:
            a=c.popleft()
            (t1,n1) = t.popleft()
            if n1 != 1:
                t.appendleft((t1,n1-1))
            simplification(Q,d,c,t,a)
    elif first_turn == d-2 and len(t)>=3: #case 3
        while first_turn==d-2:
            c.append(c.popleft())
            (t1,n1)=t.popleft()
            if n1!=1:
                t.appendleft((t1,n1-1))
            (t2,n2)=t.pop()
            if t2==2:
                t.append((2,n2+1))
            else:
                t.append((t2,n2))
                t.append((2,1))
            first_turn=t1
        if first_turn==d-1:
            a=c.popleft()
            (t1,n1) = t.popleft()
            if n1 != 1:
                t.appendleft((t1,n1-1))
            simplification(Q,d,c,t,a)

def rightpush(Q,d,c,t):
    r"""
    Perform the right push to the cyclic geodesic path c and its turn sequence t to make a canonical reprensentative c.
    """
    if len(c)==0:
        return
    elif len(c)==1:
        raise ValueError("The path wasn't geodesic")
    first_turn=turn(Q, c[-1], c[0])
    fp=Q.face_permutation(copy=False)
    vp=Q.vertex_permutation(copy=False)
    if first_turn==2 and len(t)==1 and t[0][0]==2:
        l=deque([])
        for e in c:
            e1=e+1 if e%2==0 else e-1
            l.append(fp[fp[e1]])
        c.clear()
        c.extend(l)
        m=t[0][1]
        t.clear()
        t.append((d-2,m))
        return
    length=len(c)
    i = 0
    while i<=length:
        if first_turn==1:
            if t[0][0]==2 and t[-1][0]==2:
                (t1,n1)=t.pop()
                (t2,n2)=t.popleft()
            elif t[0][0]==2:
                n1=0
                (t2,n2)=t.popleft()
            elif t[-1][0]==2:
                n2=0
                (t1,n1)=t.pop()
            else:
                n1=0
                n2=0
            if len(t)==0:
                turn_add(t,t2,n2)
                turn_add(t,t1,n1)
            elif len(t)==1 and t[0][1]==1:
                (t3,n3)=t.pop()
                if t3==3:
                    turn_add(t,d-2,n2)
                    turn_add(t,d-1,1)
                    turn_add(t,d-2,n1)
                    l=deque([])
                    e=c[0]
                    e1=e+1 if e%2==0 else e-1
                    for j in range(n2):
                        e=c[1+j]
                        e1=e+1 if e%2==0 else e-1
                        l.append(fp[fp[e1]])
                    x=fp[vp[e1]]
                    l.append(x)
                    l.append(fp[x])
                    for j in range(n1+1,1,-1):
                        e=c[-j]
                        e1=e+1 if e%2==0 else e-1
                        l.append(fp[fp[e1]])
                    c.clear()
                    c.extend(l)
                    return
                else:
                    turn_add(t,d-2,n2-1)
                    if n2 != 0:
                        turn_add(t,d-1,1)
                    turn_add(t,t3-2,1)
                    if n1 != 0:
                        turn_add(t,d-1,1)
                    turn_add(t,d-2,n1-1)
                    l=deque([])
                    e=c[0]
                    e1=e+1 if e%2==0 else e-1
                    for j in range(n2):
                        e=c[1+j]
                        e1=e+1 if e%2==0 else e-1
                        l.append(fp[fp[e1]])
                    x=vp[e1]
                    l.append(x+1 if x%2==0 else x-1)
                    y=c[n2+1]
                    l.append(fp[y+1 if y%2==0 else y-1])
                    for j in range(n1+1,1,-1):
                        e=c[-j]
                        e1=e+1 if e%2==0 else e-1
                        l.append(fp[fp[e1]])
                    c.clear()
                    c.extend(l)
                    return
            else:
                turn_modif(t,-1,d)
                if n1!=0:
                    turn_add(t,d-1,1)
                turn_add(t,d-2,n1-1)
                if n1!=0 and n2!=0:
                    first_turn=d-3
                elif n1!=0 or n2!=0:
                    first_turn=d-2
                else:
                    first_turn=d-1
                turn_modif_left(t,1,d)
                if n2!=0:
                    turn_add_left(t,d-1,1)
                turn_add_left(t,d-2,n2-1)
                l=deque()
                for _ in range(n2):
                    e=c.popleft()
                    e1=e+1 if e%2==0 else e-1
                    l.appendleft(fp[vp[e1]])
                e=c.popleft()
                e1=e+1 if e%2==0 else e-1
                e2=vp[e1]
                e3=e2+1 if e2%2==0 else e2-1
                c.extendleft(l)
                c.appendleft(e3)
                l=deque()
                e=c.pop()
                e1=e+1 if e%2==0 else e-1
                for _ in range(n1):
                    e=c.pop()
                    e1=e+1 if e%2==0 else e-1
                    l.append(fp[fp[e1]])
                c.extend(l)
                c.append(fp[e1])
        i+=1
        c.rotate()
        turn_add_left(t,first_turn,1)
        (t1,n1)=t.pop()
        first_turn=t1
        if n1!=1:
            t.append((t1,n1-1))

def reprensentative(w, Q, cor):
    r"""
    Compute the geodesic representative of w in the quad system Q with correspondance function cor.
    """
    c = deque([])
    t = deque([])
    d = 4*Q.genus()
    for e in w:
        for f in cor[e]:
            simplification(Q,d,c,t,f)
    origin_simplification(Q,d,c,t)
    rightpush(Q,d,c,t)
    return c, t

def is_homotopic(G,w1,w2,Q=None, cor=None):
    r"""
    Test whether two loops w1 and w2 are homotopic in G.
    If Q and cor are provided the representant will be computed in Q, else the function compute a new quad system.

    TEST::

        sage: G=FatGraph(vp=[2, 13, 4, 7, 6, 8, 0, 3, 10, 15, 12, 25, 14, 27, 5, 16, 18, 26, 20, 23, 9, 22, 24, 1, 21, 17, 11, 19])
        sage: is_homotopic(G, [], [])
        True
        sage: is_homotopic(G, [], [15, 8])
        False
        sage: is_homotopic(G, [], [8, 20, 24, 11])
        True
        sage: is_homotopic(G, [16, 26, 19], [9, 14])
        False
        sage: is_homotopic(G, [14, 9], [12, 1, 6, 3, 0, 23, 21, 18, 27, 17, 9])
        True
        sage: is_homotopic(G, [8, 20, 22, 13], [5, 0, 27, 11])
        False
        sage: is_homotopic(G, [15, 5, 6, 3, 4, 12, 19], [20, 24, 17, 15, 12, 1, 6, 3, 0, 19, 16, 25, 21])
        True
    """
    if Q is None or cor is None:
        treecotree=tree_co_tree(G)
        Q, cor = quad_system(G, treecotree)
    c1, t1=reprensentative(w1, Q, cor)
    c2, t2=reprensentative(w2, Q, cor)
    if len(c1)!=len(c2):
        return False
    else:
        l=c2.copy()
        l.extend(l)
        return test_KMP(c1,l)

def test_KMP(u,v):
    r"""
    Test if u is a subword of v in O(|u|+|v|).
    """
    if len(u)==0:
        return True
    elif len(v)==0:
        return False
    cnd=0
    T=[-1]
    for i in range(1,len(u)):
        if u[i] == u[cnd]:
            T.append(u[cnd])
        else:
            T.append(cnd)
            while cnd>=0 and u[i]==u[cnd]:
                cnd=T[cnd]
            cnd+=1
    j=0
    k=0
    res=False
    while j<len(v) and not res:
        if u[k]==v[j]:
            j+=1
            k+=1
            if k==len(u):
                res=True
        else:
            k=T[k]
            if k==-1:
                k+=1
                j+=1
    return res


def forms(D):
    r"""
    Compute two one forms on the edges of D such that omega wedge eta = 2
    """
    omega = [0 for _ in D.darts()]
    eta = [0 for _ in D.darts()]
    omega[0]=1
    omega[1]=-1
    fp=D.face_permutation(copy=False)
    found1=False
    e=0
    goodedge= [0 for _ in D.darts()]
    for _ in D.darts():
        e=fp[e]
        e1=e+1 if e%2==0 else e-1
        if e==1:
            found1=True
        elif found1:
            goodedge[e]-=1
            goodedge[e1]+=1
        else:
            goodedge[e]+=1
            goodedge[e1]-=1
    e=0
    while goodedge[e]!=2:
        e+=1
    eta[e]=1
    e1=e+1 if e%2==0 else e-1
    eta[e1]=-1
    return omega,eta


def area_precomuputation(G, treecotree):
    r"""
    Precompute four lists that store the information to compute the area of a contractible path.
    """
    
    F, cor, recor, rank = tree_contraction(G, treecotree)
    fpF=F.face_permutation(copy=False)
    omegaF, etaF = forms(F)
    local_area=[0 for _ in G.darts()]
    omega=[0 for _ in G.darts()]
    eta=[0 for _ in G.darts()]
    integral=[0 for _ in G.darts()]
    for e in G.darts():
        if treecotree[e]==2:
            omega[e]=omegaF[cor[e]]
            eta[e]=etaF[cor[e]]
            integral[e]=omegaF[cor[e]]*etaF[cor[e]]/2
        elif treecotree[e]==1:
            start=cor[e]
            e1=e+1 if e%2==0 else e-1
            end=cor[e1]
            current = start
            om=0
            et=0
            it=0
            alpha=0
            if current!=end or rank[e1]>rank[e]:
                la=rank[e1]-rank[e]-1
            else:
                current=fpF[current]
                om+=omegaF[current]
                et+=etaF[current]
                it+= (alpha+omegaF[current]/2)*etaF[current]
                alpha+=omegaF[current]
                la=rank[recor[current]]-rank[e]+rank[e1]-1
            while current != end:
                current=fpF[current]
                om+=omegaF[current]
                et+=etaF[current]
                it+= (alpha+omegaF[current]/2)*etaF[current]
                alpha+=omegaF[current]
                la+=rank[recor[current]]
            omega[e]=om
            eta[e]=et
            integral[e]=it
            local_area[e]=la//2+1
    return local_area, omega, eta, integral
                
def area(G, w, local_area, omega, eta, integral):
    r"""
    Compute the area of a contractible path
    """
    dom=0
    area=0
    alpha=0
    for e in w:
        dom += integral[e]+ eta[e]*alpha
        area += local_area[e]
        alpha += omega[e]
    return dom*G.num_faces()-area


def bigon_test_slow(M, e, l):
    r"""
    Test if there is a bigon in the 4-valency map M starting in the corner next to e of length l.
    If there is such a bigon, return the path formed by the boundary of this bigon starting from e.

    Should be deprecated soon to a more efficient version!
    """
    vp=M.vertex_permutation(copy=False)
    d=[e]
    f=vp[e]
    u=[f]
    for _ in range(l-1):
        e=d[-1]
        e1=e+1 if e%2==0 else e-1
        d.append(vp[vp[e1]])
        f=u[-1]
        f1=f+1 if f%2==0 else f-1
        u.append(vp[vp[f1]])
    endD=d[-1]+1 if d[-1]%2==0 else d[-1]-1
    endU=u[-1]+1 if u[-1]%2==0 else u[-1]-1
    if vp[endU] != endD:
        return False, []
    for i in range(len(u)-1,-1,-1):
        d.append(u[i]+1 if u[i]%2==0 else u[i]-1)
    return is_homotopic(M, d, []), d


def turn_update_right(t, last, e, Q, cor, d):
    r"""
    Update the turn sequence t by adding the edge e at the right.
    """
    if last[1]==None:
        last[1]=e
        last[0]=e
        return
    nt=turn(Q, last[1], e)
    if len(t)==0 and nt==0:
        last=[None, None]
    elif len(t)==0:
        turn_add(t, nt, 1)
        last[1]=e
    elif nt==0:
        (t1,n1)=t.pop()
        if n1!=1:
            t.append((t1,n1-1))
        cur=last[1]
        t2=t[-1][0]
        for _ in range(d-t2):
            cur=vp[cur]
        cur=cur+1 if cur%2==0 else cur-1
        last[1]=cur
    elif nt==1 and t[-1][0]==1:
        (t1,n1)=t.pop()
        if n1!=1:
            t.append((t1,n1-1))
        cur=last[1]
        cur=cur+1 if cur%2==0 else cur-1
        last[1]=fp[fp[cur]]
        if len(t)==0:
            last[0]=last[1]
    elif nt==d-1 and t[-1][0]==d-1:
        (t1,n1)=t.pop()
        if n1!=1:
            t.append((t1,n1-1))
        cur=last[1]
        cur=fp[fp[cur]]
        cur=cur+1 if cur%2==0 else cur-1
        last[1]=cur
        if len(t)==0:
            last[0]=last[1]
    elif len(t)>=2 and nt==1 and t[-1][0]==2 and t1[-2][0]==1:
        (t1, n1)=t.pop()
        (t2, n2)=t.pop()
        if n2!=1:
            t.append((t2,n2-1))
        t.append((d-2,n1))
        cur=last[1]
        cur=cur+1 if cur%2==0 else cur-1
        last[1]=fp[fp[cur]]
        if len(t)==1:
            cur=last[0]
            cur=cur+1 if cur%2==0 else cur-1
            last[0]=fp[fp[cur]]
    elif len(t)>=2 and nt==d-1 and t[-1][0]==d-2 and t[-2][0]==d-1:
        (t1, n1)=t.pop()
        (t2, n2)=t.pop()
        if n2!=1:
            t.append((t2,n2-1))
        t.append((2,n1))
        cur=last[1]
        cur=fp[fp[cur]]
        cur=cur+1 if cur%2==0 else cur-1
        last[1]=cur
        if len(t)==1:
            cur=last[0]
            cur=fp[fp[cur]]
            last[0]=cur+1 if cur%2==0 else cur-1
    else:
        last[1]=e
        turn_add(t, nt, 1)


def turn_update_left((t, last, e, Q, cor, d):
    r"""
    Update the turn sequence t by adding the edge e at the left.
    """
    if last[1]==None:
        last[1]=e
        last[0]=e
        return
    nt=turn(Q, e, last[0])
    if len(t)==0 and nt==0:
        last=[None, None]
    elif len(t)==0:
        turn_add_left(t, nt, 1)
        last[0]=e
    elif nt==0:
        (t1,n1)=t.popleft()
        if n1!=1:
            t.appendleft((t1,n1-1))
        cur=last[0]
        cur=cur+1 if cur%2==0 else cur-1
        t2=t[-1][0]
        for _ in range(t2):
            cur=vp[cur]
        last[0]=cur
    elif nt==1 and t[0][0]==1:
        (t1,n1)=t.popleft()
        if n1!=1:
            t.appendleft((t1,n1-1))
        cur=last[0]
        cur=cur+1 if cur%2==0 else cur-1
        last[0]=fp[fp[cur]]
        if len(t)==0:
            last[1]=last[0]
    elif nt==d-1 and t[0][0]==d-1:
        (t1,n1)=t.popleft()
        if n1!=1:
            t.appendleft((t1,n1-1))
        cur=last[0]
        cur=fp[fp[cur]]
        cur=cur+1 if cur%2==0 else cur-1
        last[0]=cur
        if len(t)==0:
            last[1]=last[0]
    elif len(t)>=2 and nt==1 and t[0][0]==2 and t1[1][0]==1:
        (t1, n1)=t.popleft()
        (t2, n2)=t.popleft()
        if n2!=1:
            t.appendleft((t2,n2-1))
        t.appendleft((d-2,n1))
        cur=last[0]
        cur=cur+1 if cur%2==0 else cur-1
        last[0]=fp[fp[cur]]
        if len(t)==1:
            cur=last[1]
            cur=cur+1 if cur%2==0 else cur-1
            last[1]=fp[fp[cur]]
    elif len(t)>=2 and nt==d-1 and t[0][0]==d-2 and t[1][0]==d-1:
        (t1, n1)=t.popleft()
        (t2, n2)=t.popleft()
        if n2!=1:
            t.appendleft((t2,n2-1))
        t.appendleft((2,n1))
        cur=last[0]
        cur=fp[fp[cur]]
        cur=cur+1 if cur%2==0 else cur-1
        last[0]=cur
        if len(t)==1:
            cur=last[1]
            cur=fp[fp[cur]]
            last[1]=cur+1 if cur%2==0 else cur-1
    else:
        last[0]=e
        turn_add_left(t, nt, 1)


def bigon_test(M, e, l, Q=None, cor=None):
    r"""
    Test if there is a bigon in the 4-valency map M in the corner next to e of length at most l.
    If there is such a bigon, it detect the one with minimal length.
    Return the path formed by the boundary of this bigon starting from e.

    INPUT:
    Q and cor is a quad system corresponding to M. If no such quad system is given, one is computed.
    """
    vp=M.vertex_permutation(copy=False)
    deg=4*M.genus()
    if Q is None or cor is None:
        treecotree=tree_co_tree(G)
        Q, cor = quad_system(G, treecotree)
    d=[e]
    f=vp[e]
    f1=f+1 if f%2==0 else f-1
    u=[f]
    turn_seq=deque([])
    last=[None, None]
    for rep in cor[f1]:
        turn_update_right(turn_seq, last, rep, Q, cor, deg)
    for rep in cor[e]:
        turn_update_right(turn_seq, last, rep, Q, cor, deg)
    endD=d[-1]+1 if d[-1]%2==0 else d[-1]-1
    endU=u[-1]+1 if u[-1]%2==0 else u[-1]-1
    i=0
    bigon=(vp[endU]==endD and len(turn_seq)==0)
    while i<l-1 and not bigon:
        e=d[-1]
        e1=e+1 if e%2==0 else e-1
        next_d=vp[vp[e1]]
        f=u[-1]
        f1=f+1 if f%2==0 else f-1
        next_u=vp[vp[f1]]
        next_u1=next_u+1 if next_u%2==0 else next_u-1
        for rep in cor[next_d]:
            turn_update_right(turn_seq, last, rep, Q, cor, deg)
        for j in range(len(cor[next_u])):
            rep=cor[next_u][-j-1]
            turn_update_left(turn_seq, last, rep, Q, cor, deg)
        d.append(next_d)
        u.append(next_u)
        i+=1
        endD=d[-1]+1 if d[-1]%2==0 else d[-1]-1
        endU=u[-1]+1 if u[-1]%2==0 else u[-1]-1
        bigon=(vp[endU]==endD and len(turn_seq)==0)
    if bigon:
        for i in range(len(u)-1,-1,-1):
            d.append(u[i]+1 if u[i]%2==0 else u[i]-1)
        return True, d
    else:
        return False, []
    


def minimal_bigon(M):
    r"""
    Find the minimal bigon in M. If no such bigon exist, this function return None.
    """
    n=M.num_vertices()
    treecotree=tree_co_tree(M)
    la, om, et, it=area_precomuputation(M,treecotree)
    bigons=[]
    area_bigons=[]
    for e in M.darts():
        l=1
        found=False
        while not found and l<=12*n:
            found,path=bigon_test_slow(M,e,l)
            l+=1
        if found:
            bigons.append(path)
            area_bigons.append(area(H, path, la, om, et, it))
    if len(bigons)==0:
        return None
    else:
        min_ind=0
        value=area_bigons[0]
        for j in range(len(bigons)):
            if value>area_bigons[j]:
                value=area_bigons[j]
                min_ind=j
        return bigons[min_ind], area_bigons[min_ind]







