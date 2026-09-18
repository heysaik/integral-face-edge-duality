"""Exhaust all globally ordered face-pairing tables through three tetrahedra."""

import json
from pathlib import Path

import regina


def matchings(items):
    if not items:
        yield []
        return
    a=items[0]
    for k in range(1,len(items)):
        b=items[k]
        for rest in matchings(items[1:k]+items[k+1:]):
            yield [(a,b)]+rest


def check(n):
    counts=dict(n=n,face_pairings=0,oriented_torus_homology_cases=0,
        no_low_degree_cases=0,zero_edge_matrix_cases=0,strict_angle_cases=0)
    for matching in matchings(list(range(4*n))):
        counts['face_pairings']+=1
        tri=regina.Triangulation3()
        for _ in range(n):
            tri.newTetrahedron()
        for a,b in matching:
            i,k=divmod(a,4)
            j,l=divmod(b,4)
            permutation=[None]*4
            permutation[k]=l
            for v,w in zip([v for v in range(4) if v!=k],
                           [w for w in range(4) if w!=l]):
                permutation[v]=w
            tri.tetrahedron(i).join(k,tri.tetrahedron(j),regina.Perm4(*permutation))
        if not (tri.isConnected() and tri.isValid() and tri.isOrientable()
                and tri.isIdeal() and tri.countVertices()==1
                and tri.countEdges()==n and tri.homology().rank()==1):
            continue
        assert tri.isOrdered()
        counts['oriented_torus_homology_cases']+=1
        # Vanishing is independent of orientation-dependent column signs.
        zero=True
        for tet in tri.tetrahedra():
            column=[0]*n
            for a,b,sign in ((0,2,1),(1,3,1),(0,3,-1),(1,2,-1)):
                column[tet.edge(regina.Edge3.edgeNumber[a][b]).index()]+=sign
            if any(column):
                zero=False
        no_low_degree=all(e.degree()>=3 for e in tri.edges())
        counts['no_low_degree_cases']+=int(no_low_degree)
        counts['zero_edge_matrix_cases']+=int(zero)
        assert not (no_low_degree and zero)
        angles=tri.hasStrictAngleStructure()
        counts['strict_angle_cases']+=int(angles)
        assert not angles or (no_low_degree and not zero)
    return counts


if __name__=='__main__':
    results=[check(n) for n in (1,2,3)]
    (Path(__file__).resolve().parent/'rank-pairing-results.json').write_text(
        json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
