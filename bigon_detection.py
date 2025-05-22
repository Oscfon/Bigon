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
    if num==0:
        return
    if len(t)==0:
        t.append((turn,num))
    else:
        if t[-1][0]==turn:
            (a,b)=t.pop()
            t.append((a,b+1))
        else:
            t.append((turn,num))

def turn_add_left(t,turn,num):
    if num==0:
        return
    if len(t)==0:
        t.append((turn,num))
    else:
        if t[0][0]==turn:
            (a,b)=t.popleft()
            t.appendleft((a,b+1))
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

    Return a boolean saying if we flatten a bracket or did a spur removal.
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
    last=c[-1]
    (last_t, n_last)=t.pop()
    newturn=turn(Q,last,e)
    if newturn==0:
        c.pop()
        if n_last!=1:
            t.append((last_t,n_last-1))
        return True
    elif newturn==1 and last_t==1:
        bracket_removal(c, t, True, 0, d)
    elif newturn==d-1 and last_t==d-1:
        bracket_removal(c, t, False, 0, d)
        return True
    elif len(t)>=2 and newturn==1 and last_t==2 and t[-1][0]==1:
        bracket_removal(c, t, True, n_last, d)
        return True
    elif len(t)>=2 and newturn==d-1 and last_t==d-2 and t[-1][0]==d-1:
        bracket_removal(c, t, False, n_last, d)
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

def bracket_removal(c, t, positive, length, d):
    if positive:
        l=deque([])
        for i in range(length+1,0,-1):
            edge=c[-i]
            edge1=edge+1 if edge%2==0 else edge-1
            l.append(fp[fp[edge1]])
        for i in range(length+2):
            c.pop()
        c=c.extend(l)
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
        c=l
        t=deque([(d-2, t[0][1]-2)])
    elif first_turn==d-1 and len(t)==1 and t[0][0]==d-2: #case 2.2
        c.pop()
        c=c.popleft()
        l=deque([])
        for e in c:
            e1=e+1 if e%2==0 else e-1
            l.append(fp[fp[e1]])
        c=l
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
            c.rotate()
            if n1==1:
                t.appendleft((t1,n1))
            else:
                t.popleft()
                t.appendleft((t1,n1))
            (t2,n2)=t.pop()
            if n2 != 1:
                t.append((t2,n2-1))
            c.reverse()
            t.reverse()
            b = c.popleft()
            (t3,n3) = t.popleft()
            if n3 != 1:
                t.appendleft((t3,n3-1))
            bremoved2=simplification(Q,d,c,t,b)
            if bremoved2:
                c.reverse()
                t.reverse()
                first_turn=turn(Q, c[-1], c[0])
                while first_turn == 0:
                    spur_removal(c,t)
                    first_turn=turn(Q, c[-1], c[0])
            else:
                c.rotate()
                if n3==1:
                    t.appendleft((t3,n3))
                else:
                    t.popleft()
                    t.appendleft((t3,n3))
                (t4,n4)=t.pop()
                if n4 != 1:
                    t.append((t4,n4-1))
                c.reverse()
                t.reverse()
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
    elif fisrt_turn == d-2 and len(t)>=3: #case 3
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

def rightpush(Q,c,t):
    r"""
    Perform the right push to the cyclic geodesic path c and its turn sequence t to make a canonical reprensentative c.
    """
    if len(c)==0:
        return
    elif len(c)==1:
        raise ValueError("The path wasn't geodesic")
    first_turn=turn(Q, c[-1], c[0])
    fp=Q.face_permutation(copy=False)
    if first_turn==2 and len(t)==1 and t[0][0]==2:
        l=deque([])
        for e in c:
            e1=e+1 if e%2==0 else e-1
            l.append(fp[fp[e1]])
        c=l
        t=deque((d-2,t[0][1]))
        return
    length=len(c)
    i = 0
    while i<length:
        if first_turn==1:
            (t1,n1)=t.pop()
            (t2,n2)=t.popleft()
            if t1==2 and t2==2 and len(t)==1:
                (t3,n3)=t.pop()
                if n3==1 and t3==3:
                    t.append((d-2,n2))
                    t.append((d-1,1))
                    t.append((d-2,n1))
                    l=deque([])
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
                    c=l
                    return
                elif n3==1:
                    turn_add(t,d-2,n2-1)
                    turn_add(t,d-1,1)
                    turn_add(t,n3-2,1)
                    turn_add(t,d-1,1)
                    turn_add(t,d-2,n1-1)
                    l=deque([])
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
                    c=l
                    return
                
        else:
            i+=1
            c.rotate()
            turn_add_left(t,first_turn,1)
            (t1,n1)=t.pop()
            first_turn=t1
            if n1!=1:
                t.append((t1,n1-1))
            
                
        
    

def reprensentative(w, Q, cor):
    r"""
    Compute the geodesic representant of w in the quad system Q with correspondance function cor.
    """
    c = []
    t = []
    d = 4*Q.genus()
    for e in w:
        for f in cor[e]:
            simplification(Q,d,c,t,f)
    origin_simplification(Q,d,c,t)
    rightpush(Q,c,t)
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
    """
    if Q is None or cor is None:
        treecotree=tree_co_tree(G)
        Q, cor = quad_system(G, treecotree)
    c1, t1=reprensentative(w1, Q, cor)
    c2, t2=reprensentative(w2, Q, cor)
    return c1==c2


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
            la=rank[e1]-rank[e]-1
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








