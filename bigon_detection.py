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