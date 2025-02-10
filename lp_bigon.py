"""

LP programming to find the smallest bigon of a 4-valency map.

"""

from surface_dynamics import *
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.numerical.mip import MIPSolverException

def valency4(G):
    res = True
    for elt in G.vertex_degrees():
        res &= (elt==4)
    return res

def rank(l,e):
    mino = 0
    app = False
    for elt in l:
        if elt == e:
            app = True
        elif elt<e:
            mino += 1
    if not app:
        raise ValueError('element not in list')
    return mino
        
def lp_k_gon(G,k, boundary_faces = None): #G is the fat graph, k an integer (0,1,2) number of corners, boundary faces the list of faces where homotopy is not allowed

    r"""
    Example for boundary_faces :
        boundary_faces=[1,3]
        x = [0,2,4,...]
    """
    
    vertices = G.vertices()
    vp = G.vertex_permutation(copy = False)
    faces = G.faces()
    fp = G.face_permutation(copy = False)
    n = len(vertices)
    f = len(faces)
    p = MixedIntegerLinearProgram(maximization = False)
    indices_f = set(range(f)) 
    if boundary_faces:
        indices_f.difference_update(boundary_faces)
    x = p.new_variable(real = True, nonnegative = True, indices = indices_f)
    alpha = p.new_variable(real = True, nonnegative = True, indices = range(4*n))
    beta = p.new_variable(real = True, nonnegative = True, indices = range(4*n))
    gamma = p.new_variable(real = True, nonnegative = True, indices = range(n))

    r"""list of constraints"""
    p.add_constraint(p.sum(alpha[i] for i in range(4*n))==k, name = "number of corners")
    p.add_constraint(2 <= 2*p.sum(x[i] for i in indices_f)-sum(beta[j] for j in range(4*n))-2*sum(gamma[j] for j in range(n)), name = "Euler characteristic")

    for i in indices_f:
        for de in faces[i]:
            for j in range(n):
                if de in vertices[j]:
                    r = rank(vertices[j],de)
                    if de%2==0:
                        d = de+1
                    else:
                        d = de-1
                    de2 = fp[d]
                    r2 = rank(vertices[j],de2)
                    p.add_constraint(x[i] == gamma[j]+alpha[4*j+r]+beta[4*j+r]+beta[4*j+r2], name = "face constraint")

    for j in range(n):
        for de in vertices[j]:
            for j2 in range(j+1,n):
                for de2 in vertices[j2]:
                    if de//2 == de2//2:
                        r = rank(vertices[j],de)
                        r2 = rank(vertices[j2], fp[de])
                        r3 = rank(vertices[j2], vp[vp[de2]])
                        p.add_constraint(alpha[4*j+r]+beta[4*j+r]==alpha[4*j2+r2]+beta[4*j2+r3], name = "edge constraint")
                        r4 = rank(vertices[j],fp[de2])
                        r5 = rank(vertices[j], vp[vp[de]])
                        r6 = rank(vertices[j2], de2)
                        p.add_constraint(alpha[4*j+r4]+beta[4*j+r5]==alpha[4*j2+r6]+beta[4*j2+r6], name = "edge constraint") 
                        p.add_constraint(beta[4*j+r4]+gamma[j]==beta[4*j2+r2]+gamma[j2], name = "edge constraint")

                        if fp[fp[fp[de]]]==de:
                            de3 = fp[de]
                            de4 = fp[de3]
                            for z in range(n):
                                if de4 in vertices[z]:
                                    j3 = z
                            r7 = rank(vertices[j3], de4)
                            p.add_constraint(beta[4*j3+r7]<=alpha[4*j+r]+gamma[j2], name = "new")
                            p.add_constraint(gamma[j3]<=gamma[j]+beta[4*j2+r3], name = "new")

    r"""optimization"""
    p.set_objective(p.sum(x[i] for i in indices_f))
    try: 
        p.solve()
    except MIPSolverException:
        raise ValueError('The following system of curves has no {}-gon'.format(k))
    
    return p, x, alpha, beta, gamma












        
    