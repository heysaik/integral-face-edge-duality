"""Check every local identity for both orientations with exact polynomials."""

from itertools import combinations

import sympy as sp
from sympy.combinatorics import Permutation


u,v,p,w,t = sp.symbols('u v p w t')
count = 0
for e in (1,-1):
    x = [u,u+v,v,v-e*t]
    y = [p-e*t,p,p+w,w]
    phi = []
    for k in range(4):
        vertices = [j for j in range(4) if j!=k]
        phi.append(dict(zip(vertices,[0,x[k],x[k]-y[k]])))
    f = [0,-x[2]-x[3],x[0]-x[2]-x[3]+y[3],x[0]-2*x[3]+y[3]]
    s = x[0]+y[3]+sp.Rational(e+1,2)*t
    expected = {(0,1):t,(2,3):t,(0,3):e*(x[0]+y[3]),
                (1,2):e*(x[0]+y[3]),(0,2):-t-e*(x[0]+y[3]),
                (1,3):-t-e*(x[0]+y[3])}
    for j,m in combinations(range(4),2):
        k,l = [a for a in range(4) if a not in (j,m)]
        hdiff = phi[k][m]-phi[k][j]-phi[l][m]+phi[l][j]
        sign = Permutation([j,m,k,l]).signature()
        q = e*sign*hdiff
        assert sp.expand(q-expected[j,m])==0
        assert sp.expand(phi[k][j]+phi[k][m]-phi[l][j]-phi[l][m]-f[l]+f[k])==0
        for vertex,other in [(j,m),(m,j)]:
            arc_sign = e*Permutation([vertex,other,k,l]).signature()
            arc = f[l]+2*phi[l][vertex]-f[k]-2*phi[k][vertex]
            assert sp.expand(arc-arc_sign*q)==0
            count += 1
        count += 2
print(f'PASS: {count} exact local edge, potential, and cusp-arc identities.')
