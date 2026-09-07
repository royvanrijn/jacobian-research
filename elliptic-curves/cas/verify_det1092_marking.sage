#!/usr/bin/env sage-python
"""Portable replay using word reduction, not Sage's Clifford implementation."""
import argparse
import json
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from sage.all import QQ, ZZ, matrix, vector, identity_matrix, QuaternionAlgebra, gcd

T = matrix(ZZ, [[-2,1,0],[1,2,2],[0,2,220]])
words = [(),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]
index = {w:i for i,w in enumerate(words)}

@lru_cache(None)
def reduce_word(w):
    for i in range(len(w)-1):
        a,b=w[i:i+2]
        if a==b:
            return QQ(T[a,a])/2*reduce_word(w[:i]+w[i+2:])
        if a>b:
            return (T[a,b]*reduce_word(w[:i]+w[i+2:])
                    -reduce_word(w[:i]+(b,a)+w[i+2:]))
    v=vector(QQ,8);v[index[w]]=1
    return v

def mul(a,b):
    return sum((a[i]*b[j]*reduce_word(words[i]+words[j])
                for i in range(8) for j in range(8) if a[i] and b[j]),vector(QQ,8))

def conjugate(a):
    return sum((a[i]*(-1)**len(words[i])*reduce_word(tuple(reversed(words[i])))
                for i in range(8) if a[i]),vector(QQ,8))

def build(data):
    assert matrix(ZZ,data['T'])==T and T.det()==-1092
    basis=[reduce_word(w) for w in [(),(0,1),(0,2),(1,2)]]
    def even(v):
        assert all(v[i]==0 for i in (1,2,3,7))
        return vector(QQ,[v[i] for i in (0,4,5,6)])
    def elt(v):
        return sum((QQ(a)*b for a,b in zip(v,basis)),vector(QQ,8))
    # Check every product by the defining quadratic relations, not its table.
    for i in range(4):
        for j in range(4):
            assert even(mul(basis[i],basis[j]))==vector(QQ,data['clifford_multiplication_table'][i][j])
    trace=matrix(QQ,4,4,lambda i,j:
                 (mul(basis[i],basis[j])+conjugate(mul(basis[i],basis[j])))[0])
    assert trace==matrix(QQ,data['trace_pairing']) and trace.det()==-546**2
    # Construct a Hilbert presentation in the independently reduced algebra.
    i=2*basis[1]-basis[0]
    j=basis[2]-QQ(2)/5*i
    assert mul(i,i)==5*basis[0]
    assert mul(j,j)==QQ(546)/5*basis[0]
    assert mul(i,j)+mul(j,i)==vector(QQ,8)
    assert matrix(QQ,[basis[0],i,j,mul(i,j)]).rank()==4
    assert QuaternionAlgebra(QQ,5,QQ(546)/5).discriminant()==546
    assert data['reduced_order_discriminant']==data['quaternion_discriminant']==546
    assert data['even_order_is_maximal'] is True

    volume=reduce_word((0,1,2))-QQ(1)/2*reduce_word((2,))-reduce_word((0,))
    V=matrix(QQ,[even(mul(volume,reduce_word((k,)))) for k in range(3)]).transpose()
    dual=vector(QQ,data['discriminant_generator'])
    assert T*dual in ZZ**3
    assert max(a.denominator() for a in dual)==1092
    q=dual*T*dual/2
    assert q==QQ(data['discriminant_quadratic_value'])
    units=[u for u in range(1092) if gcd(u,1092)==1 and (u*u-1)*q in ZZ]
    assert units==data['discriminant_orthogonal_units'] and len(units)==16
    actions={1:1}
    for r in data['prime_normalizers']:
        p=r['prime'];x=elt(r['coordinates']);n=mul(x,conjugate(x))
        assert n==p*basis[0] and r['norm']==p
        columns=[even(mul(mul(x,b),conjugate(x)/p)) for b in basis]
        A=matrix(QQ,columns).transpose()
        assert A==matrix(ZZ,r['order_action']) and abs(A.det())==1
        M=matrix(ZZ,r['T_action'])
        assert A*V==V*M and M.det()==1 and M.transpose()*T*M==T
        u=r['discriminant_unit']
        assert u in units and M*dual-u*dual in ZZ**3
        actions.update({p*m:u*v%1092 for m,v in list(actions.items())})
    assert sorted(r['prime'] for r in data['prime_normalizers'])==[2,3,7,13]
    assert {str(m):u for m,u in actions.items()}==data['atkin_lehner_actions']
    assert sorted(actions.values())==units
    stable=[m for m,u in sorted(actions.items()) if u in (1,1091)]
    assert stable==data['projectively_stable_labels']==[1,546]
    return dict(status='PASS_INDEPENDENT_WORD_REDUCTION',stable_labels=stable,
                algebra_discriminant=546,order_reduced_discriminant=546,
                verified_prime_normalizers=4,verified_discriminant_actions=16,
                boundary='Period correspondence is external; no individual non-CM point or K3 equation is certified.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=build(json.loads(args.input.read_text()))
    result.update(input_sha256=sha256(args.input.read_bytes()).hexdigest(),
                  checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    with args.output.open('x') as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))
