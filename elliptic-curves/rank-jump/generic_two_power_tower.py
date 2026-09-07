#!/usr/bin/env python3
"""Uniform congruence-module checks and proper-image cochain controls."""
import argparse
from collections import deque
from pathlib import Path
import retrospective as r

PROTOCOL=Path(__file__).with_name('GENERIC_TWO_POWER_TOWER_PROTOCOL.json')
PRIOR=r.OUT/'rank_jump_generic_fourth_division_v1.json'
OUT=r.OUT/'rank_jump_generic_two_power_tower_v1.json'
I=(1,0,0,1)
S=(0,1,1,0)
T=(0,-1,1,-1)


def mul(A,B,n):return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))%n for i in range(2) for j in range(2))
def add(A,B,n):return tuple((a+b)%n for a,b in zip(A,B))
def bits(A):return r.pack([x%2 for x in A])
def mat(v):return tuple((v>>i)&1 for i in range(4))
def tau(A):return mul(mul(T,A,2),mul(T,T,2),2)
def project(A):return add(tau(A),tau(tau(A)),2)


def module_checks():
    D=sorted({bits(project(mat(v))) for v in range(16)})
    assert len(D)==4
    assert all((mat(v)[0]+mat(v)[3])%2==0 for v in D)
    subspaces={frozenset([0])};pending=list(subspaces)
    while pending:
        H=pending.pop()
        for v in range(16):
            J=H|frozenset(x^v for x in H)
            if J not in subspaces:subspaces.add(J);pending.append(J)
    assert len(subspaces)==67
    stable=[]
    for H in subspaces:
        if not all(bits(tau(mat(v))) in H and bits(mul(mul(S,mat(v),2),S,2)) in H for v in H):continue
        E={bits(project(mat(v))) for v in H}
        assert E in ({0},set(D))
        if E!={0}:
            assert set(D)<=H
            squared={bits(project(add(mat(v),mul(mat(v),mat(v),2),2))) for v in D}
            assert squared==set(D)
        stable.append({'subspace':sorted(H),'standard_dimension':0 if E=={0} else 2})
    # Check the exact lowest-depth identity for all lifts A mod4, not just D.
    for v in range(256):
        A=tuple((v>>(2*i))&3 for i in range(4))
        X=tuple(I[i]+2*A[i] for i in range(4))
        sq=mul(X,X,8)
        residue=tuple(((sq[i]-I[i])//4)%2 for i in range(4))
        assert residue==add(A,mul(A,A,2),2)
    return {'standard_matrices':D,'invariant_subspaces':sorted(stable,key=lambda x:x['subspace']),
            'all_subspaces_checked':67,'lowest_depth_lifts_checked':256,
            'standard_component_of_squaring_is_identity':True}


def proper_control(n,depth,D):
    generators=[tuple(x%n for x in S),tuple(x%n for x in T)]
    generators += [tuple((I[i]+(2**depth)*mat(v)[i])%n for i in range(4))
                   for v in sorted(r.basis(D).values())]
    expr={I:(0,0)};pending=deque([I]);constraints=[]
    while pending:
        A=pending.popleft()
        for j,B in enumerate(generators):
            C=mul(A,B,n)
            e=tuple(expr[A][i]^sum((A[2*i+k]%2)<<(2*j+k) for k in range(2)) for i in range(2))
            if C not in expr:expr[C]=e;pending.append(C)
            else:constraints.extend(x^y for x,y in zip(expr[C],e))
    zdim=2*len(generators)-r.rank(constraints)
    assert zdim-2<=1
    solutions=[v for v in range(1<<(2*len(generators))) if all((v&w).bit_count()%2==0 for w in r.basis(constraints).values())]
    assert len(solutions)==2**zdim
    return {'modulus':n,'first_standard_depth':depth,'generators':list(map(list,generators)),
            'group_order':len(expr),'Z1_dimension':zdim,'B1_dimension':2,'H1_dimension':zdim-2,
            'cochain_assignments':solutions,'constraint_basis':list(r.basis(constraints).values())}


def calculate():
    old=r.read(PRIOR);assert old['status']=='PASS'
    for path,d in old['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==d
    mod=module_checks()
    controls=[proper_control(n,d,mod['standard_matrices']) for n,d in r.read(PROTOCOL)['controls']]
    rows=[]
    for row in old['rows']:
        m=row['generic_dimension'];assert row['split_Kummer_dimension']==m+1
        rows.append({'token':row['token'],'generic_dimension':m,
            'all_n_at_least_2_torsion_split_Kummer_dimension':1,
            'all_n_at_least_2_generic_division_split_Kummer_dimension':m+1,
            'all_n_at_least_2_generic_division_split_Selmer_dimension':m,
            'all_n_at_least_2_full_translation_rank_over_Z_mod_2n':2*m,
            'new_Selmer_dimension_in_entire_generic_two_power_tower':0,
            'actual_higher_linear_images':'NOT_COMPUTED_OR_REQUIRED_FOR_THIS_BOUND',
            'proof_basis':'uniform congruence filtration theorem plus nonzero inadmissible level-four derivative'})
    return {'schema':'rank-jump.generic-two-power-tower.v1','status':'PASS',
            'module':mod,'proper_image_controls':controls,'rows':rows,
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,PRIOR]},
            'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','check']);args=p.parse_args()
    result=calculate()
    if args.mode=='capture':r.write_new(OUT,result)
    else:assert result==r.read(OUT)
    print(result['status'],[(x['modulus'],x['group_order'],x['H1_dimension']) for x in result['proper_image_controls']])
