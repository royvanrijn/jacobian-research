#!/usr/bin/env sage-python
"""Bounded normalizer discovery and exact replay for the literal det1092 T.

This proves explicit integral Clifford/normalizer/discriminant computations.
Identifying a period curve additionally uses the maximal-order normalizer and
ternary-spin correspondence. It supplies no K3 equation or non-CM point.
"""
import argparse
import itertools
import json
from hashlib import sha256
from math import isqrt
from pathlib import Path
from sage.all import (QQ, ZZ, matrix, vector, QuadraticForm, CliffordAlgebra,
                      QuaternionAlgebra, identity_matrix, gcd)

def rows(m):
    return [list(map(str, r)) for r in m.rows()]

def build(witness=None):
    T = matrix(ZZ, [[-2, 1, 0], [1, 2, 2], [0, 2, 220]])
    assert T.det() == -1092
    C = CliffordAlgebra(QuadraticForm(QQ, T))
    e = C.gens()
    basis = [C.one(), e[0]*e[1], e[0]*e[2], e[1]*e[2]]
    keys = list(C.basis().keys())
    def coords(x):
        v = vector(QQ, [x.monomial_coefficients().get(keys[i], 0)
                        for i in (0, 4, 5, 6)])
        assert sum((v[i]*basis[i] for i in range(4)), C.zero()) == x
        return v
    def elt(v):
        return sum((QQ(v[i])*basis[i] for i in range(4)), C.zero())
    def norm(x):
        v = coords(x*x.clifford_conjugate())
        assert not any(v[1:])
        return v[0]
    table = [rows(matrix(QQ, [coords(x*y) for y in basis])) for x in basis]
    assert all(a in ZZ for x in basis for y in basis for a in coords(x*y))
    trace = matrix(ZZ, 4, 4, lambda i,j:
        coords(basis[i]*basis[j]+(basis[i]*basis[j]).clifford_conjugate())[0])
    assert trace.det() == -546**2
    quaternion = QuaternionAlgebra(QQ, 5, QQ(2184)/5)
    assert quaternion.discriminant() == 546
    # Equality of the reduced order and algebra discriminants proves maximality.
    diag = [int(norm(x)) for x in basis]
    cross = {(i,j): int(norm(basis[i]+basis[j])-diag[i]-diag[j])
             for i in range(4) for j in range(i+1,4)}
    assert diag[0] == 1
    volume = e[0]*e[1]*e[2]-QQ(1)/2*e[2]-e[0]
    images = [volume*x for x in e]
    V = matrix(QQ, [coords(x) for x in images]).transpose()
    left = (V.transpose()*V).inverse()*V.transpose()
    assert left*V == identity_matrix(QQ,3)
    smith,U,W = T.smith_form()
    assert U*T*W == smith and list(smith.diagonal()) == [1,1,1092]
    dual = T.inverse()*U.inverse()*vector(ZZ,[0,0,1])
    q = dual*T*dual/2
    units = [u for u in range(1092) if gcd(u,1092)==1 and (u*u-1)*q in ZZ]

    def check(p,v):
        x = elt(v)
        assert norm(x) == p
        inv = x.clifford_conjugate()/p
        A = matrix(QQ,[coords(x*b*inv) for b in basis]).transpose()
        assert all(a in ZZ for a in A.list()) and abs(A.det())==1
        M = matrix(QQ,[left*coords(x*b*inv) for b in images]).transpose()
        assert V*M == matrix(QQ,[coords(x*b*inv) for b in images]).transpose()
        assert all(a in ZZ for a in M.list()) and M.det()==1
        assert M.transpose()*T*M==T
        S = U*M.inverse().transpose()*U.inverse()
        assert all(S[2,j]%1092==0 for j in (0,1))
        u = int(S[2,2]%1092)
        assert M*dual-u*dual in ZZ**3
        return dict(prime=p,coordinates=list(map(int,v)),norm=p,
                    order_action=rows(A),T_action=rows(M),discriminant_unit=u)

    records = {}
    visited = 0
    if witness is None:
        # Fixed finite cube; a is solved exactly, never numerically rounded.
        for b,c,d in itertools.product(range(-24,25),repeat=3):
            visited += 1
            tail=[b,c,d]
            linear=sum(cross[0,j]*tail[j-1] for j in range(1,4))
            constant=sum(diag[j]*tail[j-1]**2 for j in range(1,4))
            constant+=sum(cross[i,j]*tail[i-1]*tail[j-1]
                          for i in range(1,4) for j in range(i+1,4))
            for p in (2,3,7,13):
                if p in records: continue
                disc=linear*linear-4*(constant-p)
                if disc<0: continue
                rt=isqrt(disc)
                if rt*rt!=disc or (rt-linear)%2: continue
                a=(rt-linear)//2
                records[p]=check(p,[a,b,c,d])
                print('normalizer',p,records[p]['coordinates'],flush=True)
            if len(records)==4: break
        assert len(records)==4, ('bounded normalizer miss',sorted(records))
    else:
        for r in witness['prime_normalizers']:
            p=int(r['prime']); records[p]=check(p,r['coordinates'])
        assert sorted(records)==[2,3,7,13]
    actions={1:1}
    for p in (2,3,7,13):
        actions.update({m*p:u*records[p]['discriminant_unit']%1092
                        for m,u in list(actions.items())})
    assert sorted(actions.values())==units and len(units)==16
    stable=[m for m,u in sorted(actions.items()) if u in (1,1091)]
    assert stable==[1,546]
    result=dict(schema='det1092.literal-marking.v1',T=rows(T),
        clifford_multiplication_table=table,trace_pairing=rows(trace),
        reduced_order_discriminant=546,quaternion_discriminant=546,
        even_order_is_maximal=True,discriminant_generator=list(map(str,dual)),
        discriminant_quadratic_value=str(q),discriminant_orthogonal_units=units,
        prime_normalizers=[records[p] for p in (2,3,7,13)],
        atkin_lehner_actions={str(m):u for m,u in sorted(actions.items())},
        projectively_stable_labels=stable,
        period_curve_with_named_correspondence='X(546)/<w546>',
        boundary='Exact integral normalizer arithmetic; period interpretation uses external theorems. No individual non-CM point, rational marked K3 equation, or new elliptic fibre is proved.')
    if witness is not None:
        assert result==witness
    return result,visited

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    args=parser.parse_args()
    assert bool(args.output)!=bool(args.verify)
    witness=json.loads(args.verify.read_text()) if args.verify else None
    result,visited=build(witness)
    if args.output:
        with args.output.open('x') as out: json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(dict(status='PASS',visited_triples=visited,
                         stable_labels=result['projectively_stable_labels'],
                         checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())))
