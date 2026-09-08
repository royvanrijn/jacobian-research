#!/usr/bin/env sage-python
"""Bounded torsion-module obstruction from the ten immutable RR sextics."""
import hashlib,json,signal
from itertools import product
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix,vector,zero_matrix
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
SRC=ART/'det1092_marked_kummer_transport_v1.json'
OUT=ART/'det1092_rr_torsion_transport_gate_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    source=json.loads(SRC.read_text());assert len(source['cases'])==10
    F=GF(2);permutation=[1,2,3,4,0,5];columns=[]
    for i in range(4):
        v=vector(F,6);v[permutation[i]]+=1;v[permutation[5]]+=1
        # In even subsets modulo complementation, choose representative v4=0.
        v+=vector(F,[v[4]]*6);assert v[4]==0 and sum(v)==0
        columns.append(v[:4])
    J=matrix(F,columns).transpose();R=PolynomialRing(F,'z');z=R.gen()
    cp=J.charpoly('z');assert cp==z**4+z**3+z**2+z+1 and cp.is_irreducible()
    assert J**5==matrix.identity(F,4)
    intertwining=[]
    for entries in product([0,1],repeat=4):
        E=matrix(F,2,2,entries)
        if not E.is_invertible():continue
        constraints=[];reverse=[]
        for i in range(8):
            T=zero_matrix(F,2,4);T[i//4,i%4]=1
            constraints.append(vector(F,(T*J-E*T).list()))
            U=zero_matrix(F,4,2);U[i//2,i%2]=1
            reverse.append(vector(F,(J*U-U*E).list()))
        M=matrix(F,constraints).transpose();N=matrix(F,reverse).transpose()
        assert M.rank()==N.rank()==8
        intertwining.append({'elliptic_action':[[int(a) for a in row] for row in E.rows()],
                             'J_to_E_constraints':[[int(a) for a in row] for row in M.rows()],
                             'E_to_J_constraints':[[int(a) for a in row] for row in N.rows()],
                             'J_to_E_rank':8,'E_to_J_rank':8,'Hom_dimensions':[0,0]})
    assert len(intertwining)==6
    cases=[]
    for row in source['cases']:
        p=row['theta_primes']['one_plus_five'];R=PolynomialRing(GF(p),'T')
        q=R([QQ(a) for a in row['q']]);assert q.degree()==6 and q.gcd(q.derivative())==1
        fac=q.factor();assert sorted(int(g.degree()) for g,e in fac)==[1,5]
        assert all(e==1 for g,e in fac)
        cases.append({'label':row['label'],'q':row['q'],'p':p,
                      'mod_p_factors':[{'coefficients':[int(a) for a in g], 'multiplicity':int(e)} for g,e in fac],
                      'J2_GQ_irreducible':True,'Q_simple':True,
                      'Hom_GQ_J2_E2':'ZERO for every elliptic curve over Q',
                      'Hom_Q_J_E':'ZERO for every elliptic curve over Q'})
    return {'classification':'verified application and new deduction',
            'status':'PASS_TEN_RR_JACOBIANS_Q_SIMPLE_AND_NO_ELLIPTIC_TORSION_TRANSPORT',
            'sage_version':sage_version,'permutation':permutation,
            'torsion_model':'Even subsets of the six branch points modulo complementation',
            'basis':'e0+e5,e1+e5,e2+e5,e3+e5, with e4+e5 equal to their sum',
            'J2_action':[[int(a) for a in row] for row in J.rows()],
            'characteristic_polynomial':[int(a) for a in cp],
            'elliptic_intertwining_tests':intertwining,'cases':cases,
            'deductions':{'Q_simple':'An elliptic subvariety over Q would give a GQ-stable two-dimensional subspace of J[2], impossible by irreducibility.',
                          'no_elliptic_quotient':'A nonzero J-to-elliptic homomorphism has an elliptic connected kernel over Q, contradicting Q-simplicity.',
                          'no_curve_to_elliptic_map':'With the inherited rational basepoint, a nonconstant RR-curve-to-elliptic map would factor through a nonzero Jacobian homomorphism.',
                          'no_coefficient_induced_Selmer_transport':'The zero Hom_GQ(J[2],E[2]) group excludes nonzero maps induced on H1 and Selmer groups by a torsion-module homomorphism.'},
            'scope':'Over Q only. Absolute simplicity, maps after number-field extension, arbitrary maps between finite Selmer groups, and nonlinear incidence tests are not decided. Selmer groups themselves are not computed.',
            'limits':{'wall_seconds':25,'sextics':10,'retained_prime_witnesses':10,'new_prime_searches':0,
                      'GL2_F2_actions':6,'point_searches':0,'Selmer_runs':0,'class_group_runs':0,'pilot_changes':0},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [SRC,Path(__file__)]}}
if __name__=='__main__':
    signal.alarm(25);d=build();payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(d['status'],flush=True)
