#!/usr/bin/env python3
"""Independent matrix modules and DFS/Sage cochain kernels for tower controls."""
import argparse
from pathlib import Path
from sage.all import GF,VectorSpace,matrix,vector
import retrospective as r

INPUT=r.OUT/'rank_jump_generic_two_power_tower_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_two_power_tower_verification_v1.json'


def verify_module(saved):
    F=GF(2);W=VectorSpace(F,4)
    S=matrix(F,[[0,1],[1,0]]);T=matrix(F,[[0,1],[1,1]])
    def mat(v):return matrix(F,2,2,list(v))
    def vec(A):return vector(F,A.list())
    def tr(v):return vec(T*mat(v)*T.inverse())
    def e1(v):return tr(v)+tr(tr(v))
    D=W.subspace([e1(v) for v in W.basis()]);assert D.dimension()==2
    assert all(mat(v).trace()==0 and e1(vec(mat(v)+mat(v)**2))==v for v in D)
    stable=[];total=0
    for d in range(5):
        for H in W.subspaces(d):
            total+=1
            if not all(tr(v) in H and vec(S*mat(v)*S) in H for v in H.basis()):continue
            dim=W.subspace([e1(v) for v in H.basis()]).dimension()
            assert dim in (0,2)
            if dim:assert D.is_subspace(H)
            stable.append({'subspace':sorted(r.pack(v) for v in H),'standard_dimension':dim})
    assert total==saved['all_subspaces_checked']==67
    assert sorted(stable,key=lambda x:x['subspace'])==saved['invariant_subspaces']
    assert sorted(r.pack(v) for v in D)==saved['standard_matrices']
    return {'subspaces':total,'stable_subspaces':len(stable),'standard_dimension':2}


def verify_control(saved):
    n=saved['modulus'];gens=[tuple(x) for x in saved['generators']];identity=(1,0,0,1)
    def product(A,B):
        a,b,c,d=A;e,f,g,h=B
        return ((a*e+b*g)%n,(a*f+b*h)%n,(c*e+d*g)%n,(c*f+d*h)%n)
    expressions={identity:(0,0)};pending=[identity];equations=[]
    # DFS gives different word representatives from the worker's BFS.
    while pending:
        A=pending.pop()
        for j,B in enumerate(gens):
            AB=product(A,B)
            values=[]
            for i in range(2):
                value=expressions[A][i]
                if A[2*i]%2:value^=1<<(2*j)
                if A[2*i+1]%2:value^=1<<(2*j+1)
                values.append(value)
            if AB not in expressions:
                expressions[AB]=tuple(values);pending.append(AB)
            else:equations.extend(x^y for x,y in zip(values,expressions[AB]))
    assert len(expressions)==saved['group_order']
    M=matrix(GF(2),[[int(v>>i&1) for i in range(2*len(gens))] for v in set(equations)])
    ker=M.right_kernel();solutions=sorted(r.pack(v) for v in ker)
    assert solutions==saved['cochain_assignments']
    assert ker.dimension()==saved['Z1_dimension']
    cob=[]
    for v in [(0,0),(1,0),(0,1),(1,1)]:
        cob.append(r.pack([(A[2*i]*v[0]+A[2*i+1]*v[1]+v[i])%2 for A in gens for i in range(2)]))
    assert len(set(cob))==4 and all(v in solutions for v in cob)
    assert len(solutions)<=8
    return {'modulus':n,'group_order':len(expressions),'H1_dimension':int(ker.dimension())-2,'status':'PASS'}


def verify():
    x=r.read(INPUT)
    for path,d in x['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==d
    mod=verify_module(x['module']);controls=[verify_control(c) for c in x['proper_image_controls']]
    assert len(x['rows'])==16 and all(c['new_Selmer_dimension_in_entire_generic_two_power_tower']==0 for c in x['rows'])
    return {'schema':'rank-jump.generic-two-power-tower-verification.v1','status':'PASS',
        'module':mod,'proper_image_controls':controls,'panel_rows':16,
        'proof_boundary':'Independent finite checks support the written arbitrary-depth filtration proof; finite controls alone do not prove it.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),INPUT]}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','check']);args=p.parse_args();result=verify()
    if args.mode=='capture':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],result['module'],[c['H1_dimension'] for c in result['proper_image_controls']])
