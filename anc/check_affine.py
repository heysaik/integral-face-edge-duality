"""Exact convention check using the figure-eight data printed in RIMS 2026.

This finite calculation checks the affine constants and face map. It does not
prove the general face-gluing theorem or existence for other triangulations.
"""

import sympy as sp


def check():
    face_table=[[1,3,2,0],[2,0,1,3]]
    xs=[sp.zeros(2,4) for _ in range(4)]
    for i,row in enumerate(face_table):
        for k,face in enumerate(row):
            xs[k][i,face]=1
    c=(xs[0]-xs[1]+xs[2]).col_join(xs[2]-xs[3])
    cb=(xs[3]-xs[2]+xs[1]).col_join(xs[1]-xs[0])
    e=sp.diag(1,-1)
    delta=(e+sp.eye(2))/2
    d=sp.zeros(2).col_join(e)
    g=sp.Matrix([[2,0],[0,0]])
    gp=sp.Matrix([[1,2],[0,-4]])
    gpp=sp.Matrix([[0,1],[0,-2]])
    a,b=g-gp,gpp-gp
    kappa=sp.Matrix([2,-2])
    eta=kappa-gp*sp.ones(2,1)
    p0=sp.ones(2,1)/3
    r0=-p0
    q=sp.Matrix([sp.Rational(1,2),sp.Rational(1,4)])
    r=sp.zeros(2,1)
    assert eta==sp.Matrix([-1,2])
    assert a*p0-b*r0==eta
    assert a*q-b*r==sp.Matrix([0,1])
    assert c.det()==-1 and b.det()==-2
    face_relations=(xs[0]-xs[1]+xs[2]).row_join(sp.zeros(2,4)).col_join(
        sp.zeros(2,4).row_join(xs[3]-xs[2]+xs[1])).col_join(
        (xs[2]-xs[3]).row_join(xs[0]-xs[1]))
    t=e*(xs[2]-xs[3])
    projection=t.row_join(sp.zeros(2,4)).col_join(
        (xs[0]+delta*t).row_join(xs[3]))
    face_basis=sp.Matrix.hstack(*face_relations.nullspace())
    image=projection*face_basis
    assert face_relations.rank()==6
    assert face_basis.cols==image.rank()==2
    assert a.row_join(-b)*image==sp.zeros(2)
    m=xs[0]*c.inv()*d
    assert xs[3]*cb.inv()*d==m.T
    assert b.inv()*a==m+m.T+delta
    z_log=sp.I*sp.pi*p0
    s_log=sp.I*sp.pi*r0
    zp_log=sp.I*sp.pi*sp.ones(2,1)-z_log+s_log
    zpp_log=-s_log
    assert g*z_log+gp*zp_log+gpp*zpp_log==sp.I*sp.pi*kappa
    assert z_log==zp_log==zpp_log
    z=(1+sp.I*sp.sqrt(3))/2
    assert sp.simplify(z*(1-z)-1)==0
    assert sp.simplify(1/(1-z)-z)==0
    assert sp.simplify(1-1/z-z)==0
    meridian_b=sp.Matrix([[-1,-1],[1,1]])
    assert meridian_b.rank()==1
    print('PASS: exact affine constants, face isomorphism, regular shapes, and peripheral-choice counterexample.')


if __name__=='__main__':
    check()
