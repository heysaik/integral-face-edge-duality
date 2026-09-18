"""Independent exact checks for the face-cochain manuscript.

Run with Python 3.13 after installing requirements.txt. This is an exact
finite verification, not a proof of the universal mathematical assertions.
"""

import json
from pathlib import Path

import regina
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from famed_exact import face_matrices, generalized_data, nz_matrices


def from_faces(table):
    tri=regina.Triangulation3()
    for _ in table:
        tri.newTetrahedron()
    occurrences={}
    for i,row in enumerate(table):
        for k,face in enumerate(row):
            occurrences.setdefault(face,[]).append((i,k))
    assert all(len(v)==2 for v in occurrences.values())
    for (i,k),(j,l) in occurrences.values():
        source=[a for a in range(4) if a!=k]
        target=[a for a in range(4) if a!=l]
        permutation=[None]*4
        permutation[k]=l
        for a,b in zip(source,target):
            permutation[a]=b
        tri.tetrahedron(i).join(k,tri.tetrahedron(j),regina.Perm4(*permutation))
    assert tri.isValid() and tri.isOrdered() and tri.isOrientable()
    assert tri.isIdeal() and tri.countVertices()==1
    assert str(tri.homology())=='Z'
    return tri


def invariants(matrix):
    normal=smith_normal_form(matrix,domain=sp.ZZ)
    return [abs(int(normal[i,i])) for i in range(min(normal.shape))]


def integer_kernel(matrix):
    """A saturated integer basis, using unimodular Smith transformations."""
    normal,left,right=smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(sp.ZZ))
    normal,left,right=(m.to_Matrix() for m in (normal,left,right))
    assert normal==left*matrix*right
    assert abs(left.det())==abs(right.det())==1
    rank=sum(normal[i,i]!=0 for i in range(min(normal.shape)))
    basis=right[:,rank:]
    assert matrix*basis==sp.zeros(matrix.rows,basis.cols)
    return basis


def modular_rank(matrix,prime):
    rows=[[int(v)%prime for v in row] for row in matrix.tolist()]
    rank=0
    for column in range(matrix.cols):
        pivot=next((i for i in range(rank,matrix.rows) if rows[i][column]),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inverse=pow(rows[rank][column],-1,prime)
        rows[rank]=[(v*inverse)%prime for v in rows[rank]]
        for i in range(rank+1,matrix.rows):
            factor=rows[i][column]
            rows[i]=[(a-factor*b)%prime for a,b in zip(rows[i],rows[rank])]
        rank+=1
        if rank==matrix.rows:
            break
    return rank


def valuation_list(smith,prime):
    valuations=[]
    for factor in smith:
        if not factor:
            continue
        exponent=0
        while factor%prime==0:
            factor//=prime
            exponent+=1
        if exponent:
            valuations.append(exponent)
    return sorted(valuations)


def check(specimen):
    tri=from_faces(specimen['face_table'])
    assert tri.isoSig()==specimen['isosig']
    n=tri.size()
    xs,c,d,eps=face_matrices(tri)
    cb=(xs[3]-xs[2]+xs[1]).col_join(xs[1]-xs[0])
    delta=(eps+sp.eye(n))/2
    joint=c.row_join(sp.zeros(2*n)).row_join(-d).col_join(
        sp.zeros(2*n).row_join(cb).row_join(-d))
    lifts=c.row_join(sp.zeros(2*n)).col_join(sp.zeros(2*n).row_join(cb)).col_join(
        xs[0].row_join(xs[3]))
    projection=sp.zeros(n,4*n).row_join(sp.eye(n)).col_join(
        xs[0].row_join(xs[3]).row_join(delta))
    domain=sp.Matrix.hstack(*joint.nullspace())
    out=projection*domain
    a,b,manifold,slope=nz_matrices(tri)
    assert domain.cols==out.rank()==n
    assert a.row_join(-b)*out==sp.zeros(n,n)
    assert invariants(c)==invariants(cb)
    assert invariants(joint)==[1]*(4*n)
    assert invariants(lifts)==[1]*(4*n)
    assert invariants(a.row_join(b))==[1]*(n-1)+[2]
    assert 2*n-c.rank()==specimen['expected_null_c']
    assert abs(int(c.det()))==specimen['expected_absdet_c']
    assert n-b.rank()==2*(2*n-c.rank())
    assert abs(b.det())==2*c.det()**2
    kernel_c=integer_kernel(c)
    kernel_cb=integer_kernel(cb)
    kernel_b=integer_kernel(b)
    face_kernel=(xs[0]*kernel_c).row_join(xs[3]*kernel_cb)
    assert face_kernel.cols==face_kernel.rank()==kernel_b.cols
    assert hermite_normal_form(face_kernel)==hermite_normal_form(kernel_b)
    smith_c,smith_b=invariants(c),invariants(b)
    assert sum(factor>1 for factor in smith_c)<=1
    modular_checks={}
    for prime in (2,3,5,7,11):
        null_c=2*n-modular_rank(c,prime)
        null_b=n-modular_rank(b,prime)
        if prime==2:
            assert 2*null_c<=null_b<=2*null_c+1
        else:
            assert null_b==2*null_c
            assert valuation_list(smith_b,prime)==sorted(2*valuation_list(smith_c,prime))
        modular_checks[str(prime)]=dict(null_c=null_c,null_b=null_b)
    if c.det():
        m=xs[0]*c.inv()*d
        assert xs[3]*cb.inv()*d==m.T
        assert b.inv()*a==m+m.T+delta
    general,_=generalized_data(tri)
    assert general['rows_equal'] and general['forms_equal']
    assert general['constraint_rank']==2*general['null_a']
    return dict(name=specimen['name'],n=n,null_c=general['null_a'],
        absdet_c=abs(int(c.det())),absdet_b=abs(int(b.det())),
        exact_checks='PASS',strict_angles=tri.hasStrictAngleStructure(),
        longitude_slope=slope,smith_c=smith_c,smith_b=smith_b,
        integer_kernel_decomposition='PASS',modular_checks=modular_checks)


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    specimens=json.loads((here/'specimens.json').read_text())
    results=[check(s) for s in specimens]
    (here/'verification-results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
    print(f'PASS: {len(results)} exact fixed-specimen verifications.')
