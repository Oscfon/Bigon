from flatsurf import *
from surface_dynamics import *


def labels(G):
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
    fp=G.face_permutation(copy=False)
    correspondence = [None for _ in G.darts()]
    nfp=[None for _ in range(4*G.genus())]
    i=0
    while treecotree[i]!=2:
        i+=2
    last=0
    correspondence[i]=0
    correspondence[i+1]=1
    e=fp[i]
    free_index=2
    
    while e!=i:
        e1=e+1 if e%2==0 else e-1
        if treecotree[e]==1:
            correspondence[e]=last
            e=fp[e1]
        elif treecotree[e]==0:
            e=fp[e]
        else:
            if correspondence[e]==None:
                correspondence[e]=free_index
                correspondence[e1]=free_index+1
                free_index+=2
            nfp[last]=correspondence[e]
            last=correspondence[e]
            e=fp[e]
    nfp[last]=0
    F=FatGraph(fp=nfp)
    return F,correspondence
            
def quad_system(G, treecotree):
    (F,cor)=tree_contraction(G,treecotree)
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
            start=qfp[cor[e]]//2
            e1=e+1 if e%2==0 else e-1
            end=qfp[cor[e1]]//2
            das=2*remember2[start]+1
            dae=2*remember2[end]
            cor2.append([das,dae])
        else:
            edge=cor[e]
            da=2*remember2[edge]+1
            cor2.append([da,qvp[da-1]])
    return Q,cor2


def turn(Q,e,f):
    vp=Q.vertex_permutation(copy=False)
    i=0
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
    if len(t)==0:
        t.append((turn,num))
    else:
        if t[-1][0]==turn:
            t[-1][1]+=num
        else:
            t.append((turn,num))

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
    fp=Q.face_permutation(copy=False)
    if len(c)==0:
        c.append(e)
        return
    last=c[-1]
    (last_t, n_last)=t.pop()
    newturn=turn(Q,last,e)
    if newturn==0:
        c.pop()
        if n_last!=1:
            t.append((last_t,n_last-1))
    elif newturn==1 and last_t==1:
        c.pop()
        edge=c.pop()
        edge1=edge+1 if edge%2==0 else edge-1
        c.pop()
        c.append(fp[fp[edge1]])
        if n_last!=1:
            raise ValueError("The path wasn't geodesic")
        turn_modif(t,-1,d)
    elif newturn==d-1 and last_t==d-1:
        c.pop()
        edge=fp[fp[c.pop()]]
        edge1=edge+1 if edge%2==0 else edge-1
        c.pop()
        c.append(edge1)
        if n_last!=1:
            raise ValueError("The path wasn't geodesic")
        turn_modif(t,1,d)
    elif len(t)>=2 and newturn==1 and last_t==2 and t[-1][0]==1:
        c.pop()
        l=[]
        for _ in range(n_last+1):
            edge=c.pop()
            edge1=edge+1 if edge%2==0 else edge-1
            l.append(fp[fp[edge1]])
        c.pop()
        c=c+l
        (t2,n2)=t.pop()
        if n2!=1:
            raise ValueError("The path wasn't geodesic")
        turn_modif(t,-1,d)
        turn_add(t,d-2,n_last)
    elif len(t)>=2 and newturn==d-1 and last_t==d-2 and t[-1][0]==d-1:
        c.pop()
        l=[]
        for _ in range(n_last+1):
            edge=fp[fp[c.pop()]]
            edge1=edge+1 if edge%2==0 else edge-1
            l.append(edge1)
        c.pop()
        c=c+l
        (t2,n2)=t.pop()
        if n2!=1:
            raise ValueError("The path wasn't geodesic")
        turn_modif(t,1,d)
        turn_add(t,2,n_last)
    elif newturn==last_t:
        c.append(e)
        t.append(newturn,n_last+1)
    else:
        c.append(e)
        t.append((last_t,n_last))
        t.append((newturn,1))


def representant(w, Q, cor):
    c = []
    t = []
    d = 4*Q.genus()
    for e in w:
        for f in cor[e]:
            simplification(Q,d,c,t,e)
    return c, t

def is_homotopic(G,w1,w2,Q=None, cor=None):
    if Q is None:
        treecotree=tree_co_tree(G)
        Q, cor = quad_system(G, treecotree)
    c1, t1=representant(w1, Q, cor)
    c2, t2=representant(w2, Q, cor)
    return c1==c2


