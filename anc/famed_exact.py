"""Exact exploratory FAMED matrices. No numerical rank decisions.

Definitions: Ben Aribi--Guilloux--Wong, arXiv:2512.17437, section 2.
SnapPy column convention follows upstream FAMEDExploration/NZ_data.py.
Geometricity is deliberately not inferred from a strict angle structure.
"""

from itertools import product
from math import gcd

import regina
import snappy
import sympy as sp


EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def ordered_versions(tri, quotient_reversal=True):
    """Enumerate global edge orientations, retain transitive tournaments."""
    n = tri.size()
    incidence = [[None] * 6 for _ in range(n)]
    for edge in tri.edges():
        for emb in edge.embeddings():
            u, v = emb.vertices()[0], emb.vertices()[1]
            local = EDGES.index(tuple(sorted((u, v))))
            incidence[emb.simplex().index()][local] = (edge.index(), 1 if u < v else -1)
    offset = 1 if quotient_reversal else 0
    for tail in product((1, -1), repeat=tri.countEdges() - offset):
        signs = (1,) + tail if quotient_reversal else tail
        iso = regina.Isomorphism3(n)
        for i, inc in enumerate(incidence):
            indegree = [0] * 4
            for (u, v), (edge, relative) in zip(EDGES, inc):
                indegree[v if signs[edge] * relative == 1 else u] += 1
            if sorted(indegree) != [0, 1, 2, 3]:
                break
            iso.setTetImage(i, i)
            iso.setFacePerm(i, regina.Perm4(*indegree))
        else:
            ordered = iso(tri)
            assert ordered.isOrdered()
            yield ordered


def face_matrices(tri):
    assert tri.isOrdered() and tri.isOrientable() and tri.isIdeal()
    n = tri.size()
    xs = [sp.zeros(n, tri.countTriangles()) for _ in range(4)]
    for tet in tri.tetrahedra():
        for k in range(4):
            xs[k][tet.index(), tet.triangle(k).index()] = 1
    acal = (xs[0] - xs[1] + xs[2]).col_join(xs[2] - xs[3])
    eps = sp.diag(*[tet.orientation() for tet in tri.tetrahedra()])
    bcal = sp.zeros(n).col_join(eps)
    return xs, acal, bcal, eps


def longitude(manifold):
    """Primitive rationally null-homologous slope for a one-cusped b1=1 M.

    For knot exteriors H1=Z, this is the preferred longitude up to sign.
    Uses the same abelianization principle as SnapPy's Sage-only routine.
    """
    group = manifold.fundamental_group()
    gens = group.generators()
    def exponents(word):
        return sp.Matrix([[word.count(g) - word.count(g.upper()) for g in gens]])
    relations = sp.zeros(0, len(gens))
    for word in group.relators():
        relations = relations.col_join(exponents(word))
    basis = relations.nullspace()
    assert len(basis) == 1, "Requires first Betti number one"
    mer, lon = group.peripheral_curves()[0]
    m = (exponents(mer) * basis[0])[0]
    l = (exponents(lon) * basis[0])[0]
    denom = sp.ilcm(m.q, l.q)
    m, l = int(m * denom), int(l * denom)
    d = gcd(m, l)
    assert d
    slope = (l // d, -m // d)
    if slope[0] < 0 or (slope[0] == 0 and slope[1] < 0):
        slope = tuple(-v for v in slope)
    return slope


def nz_matrices(tri):
    manifold = snappy.Manifold(tri)
    n = tri.size()
    raw = manifold.gluing_equations()
    nz = sp.Matrix([[int(x) for x in row] for row in raw])
    p, q = longitude(manifold)
    eq = nz[:n-1, :].col_join(p * nz[-2, :] + q * nz[-1, :])
    g, gp, gpp = [eq[:, list(range(k, 3*n, 3))] for k in (0, 2, 1)]
    a, b = g - gp, gpp - gp
    assert a * b.T == b * a.T
    assert a.row_join(b).rank() == n
    return a, b, manifold, (p, q)


def graph_reduction(tri):
    """Collapse face equalities f2=f3 by union-find."""
    n = tri.size()
    parents = list(range(2*n))
    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    cycles = 0
    for tet in tri.tetrahedra():
        a = find(tet.triangle(2).index())
        b = find(tet.triangle(3).index())
        if a == b:
            cycles += 1
        else:
            parents[a] = b
    roots = sorted({find(i) for i in range(2*n)})
    component = [roots.index(find(i)) for i in range(2*n)]
    c = sp.zeros(n, len(roots))
    for tet in tri.tetrahedra():
        for k, sign in ((0, 1), (1, -1), (2, 1)):
            c[tet.index(), component[tet.triangle(k).index()]] += sign
    return c, cycles, component


def analyze(tri, with_nz=True):
    n = tri.size()
    xs, acal, bcal, eps = face_matrices(tri)
    c, cycles, _ = graph_reduction(tri)
    null_a = acal.cols - acal.rank()
    assert null_a == c.cols - c.rank()
    result = dict(n=n, isosig=tri.isoSig(), null_a=null_a, cycles=cycles,
                  reduced_row_defect=n-c.rank(), det_a=int(acal.det()),
                  strict_angles=tri.hasStrictAngleStructure())
    if with_nz:
        a, b, manifold, slope = nz_matrices(tri)
        null_b = b.cols - b.rank()
        result.update(null_b=null_b, det_b=int(b.det()), slope=slope,
                      rank_conjecture=(null_b == 2*null_a),
                      numerical_solution=manifold.solution_type())
        if null_a == null_b == 0:
            m = xs[0] * acal.inv() * bcal
            s = m + m.T + (eps + sp.eye(n))/2
            result['famed_identity'] = s == b.inv()*a
            if not result['famed_identity']:
                result['residual'] = str(s - b.inv()*a)
    return result


def generalized_inverse(mat):
    """A rational reflexive generalized inverse from a full-rank minor."""
    cols = list(mat.rref()[1])
    rows = list(mat.T.rref()[1])
    out = sp.zeros(mat.cols, mat.rows)
    if cols:
        minor_inverse = mat.extract(rows, cols).inv()
        for i, c in enumerate(cols):
            for j, r in enumerate(rows):
                out[c, r] = minor_inverse[i, j]
    assert mat * out * mat == mat
    return out


def generalized_data(tri):
    xs, c, d, eps = face_matrices(tri)
    a, b, manifold, slope = nz_matrices(tri)
    n = tri.size()
    k = sp.Matrix.hstack(*c.nullspace()) if c.nullspace() else sp.zeros(2*n, 0)
    left = sp.Matrix.hstack(*c.T.nullspace()) if c.T.nullspace() else sp.zeros(2*n, 0)
    constraint = (xs[0]*k).T.col_join(left.T*d)
    g = xs[0] * generalized_inverse(c) * d
    s = g + g.T + (eps+sp.eye(n))/2
    left_b = sp.Matrix.hstack(*b.T.nullspace()) if b.T.nullspace() else sp.zeros(n, 0)
    nz_constraint = left_b.T*a
    null_a = c.cols-c.rank()
    null_b = b.cols-b.rank()
    w = sp.Matrix.hstack(*constraint.nullspace()) if constraint.nullspace() else sp.zeros(n, 0)
    nz_s = generalized_inverse(b)*a
    rows_equal = constraint.rank() == nz_constraint.rank() == constraint.col_join(nz_constraint).rank()
    forms_equal = w.T*(s-nz_s)*w == sp.zeros(w.cols)
    return dict(null_a=null_a, null_b=null_b, constraint_rank=constraint.rank(),
                rows_equal=rows_equal, forms_equal=forms_equal,
                rank_conjecture=(null_b==2*null_a)), (constraint,s,a,b,w)
