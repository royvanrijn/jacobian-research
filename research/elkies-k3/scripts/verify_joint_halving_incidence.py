#!/usr/bin/env python3
"""Finite monodromy and genus premises for simultaneous branch incidence.
No elliptic-curve search or low-degree-point enumeration is performed.
"""
import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import resource
import time


def compose(a,b):
    return tuple(a[b[i]] for i in range(len(a)))


def generated(gens):
    identity=tuple(range(len(gens[0])))
    group={identity}; todo=[identity]
    while todo:
        g=todo.pop()
        for h in gens:
            z=compose(g,h)
            if z not in group:group.add(z);todo.append(z)
    return group


def cycle_lengths(g):
    remaining=set(range(len(g)));lengths=[]
    while remaining:
        i=min(remaining);n=0
        while i in remaining:remaining.remove(i);n+=1;i=g[i]
        lengths.append(n)
    return sorted(lengths)


def verify():
    resource.setrlimit(resource.RLIMIT_CPU,(10,10))
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    start=time.monotonic()
    linear=[(0,*p) for p in permutations((1,2,3))]
    assert all(g[x^y]==g[x]^g[y] for g in linear for x,y in product(range(4),repeat=2))
    single={tuple(g[x]^a for x in range(4)) for g in linear for a in range(4)}
    assert len(single)==24 and single==set(permutations(range(4)))
    # Every invariant subspace of V+V; no semisimplicity assumption is hidden.
    subspaces={frozenset([0])};todo=list(subspaces)
    while todo:
        W=todo.pop()
        for v in range(16):
            if v in W:continue
            U=W|frozenset(x^v for x in W)
            if U not in subspaces:subspaces.add(U);todo.append(U)
    assert len(subspaces)==67
    def action(g,x):return g[x&3]+4*g[x>>2]
    invariant=[W for W in subspaces if all(action(g,x) in W for g in linear for x in W)]
    assert sorted(map(len,invariant))==[1,4,4,4,16]
    for W in invariant:
        if len(W)==16:continue
        assert any(all(((x&3) if a else 0)^((x>>2) if b else 0)==0 for x in W)
                   for a,b in [(1,0),(0,1),(1,1)])
    # Optional direct cross-check of H^1(S3,V)=0: all cocycles are coboundaries.
    ident=tuple(range(4));linear.sort();identity_index=linear.index(ident)
    table=[[linear.index(compose(g,h)) for h in linear] for g in linear]
    cocycles=[]
    for c in product(range(4),repeat=6):
        if c[identity_index]:continue
        if all(c[table[i][j]]==c[i]^linear[i][c[j]] for i,j in product(range(6),repeat=2)):
            cocycles.append(c)
    coboundaries={tuple(v^g[v] for g in linear) for v in range(4)}
    assert set(cocycles)==coboundaries and len(cocycles)==4
    joint={tuple(action(g,x)^a for x in range(16)) for g in linear for a in range(16)}
    assert len(joint)==96
    H=[g for g in joint if g[0]==0]
    H1=[g for g in joint if (g[0]&3)==0]
    assert len(H)==6 and len(H1)==24
    # Stabilizer of one half has a primitive four-point action on the other.
    for g in H1:
        if g not in H:assert generated(H+[g])==set(H1)
    assert {g[0] for g in joint}==set(range(16))
    mixed_points=list(product(range(4),range(1,4)))
    mixed={tuple(mixed_points.index((g[x]^a,g[e])) for x,e in mixed_points)
           for g in linear for a in range(4)}
    assert len(mixed)==24 and {g[0] for g in mixed}==set(range(12))
    inertia=(0,2,1,3)
    representations={
        'single_half':inertia,
        'nonzero_two_torsion':tuple([1,2,3].index(inertia[x]) for x in [1,2,3]),
        'two_independent_halves':tuple(action(inertia,x) for x in range(16)),
        'half_and_nonzero_two_torsion':tuple(mixed_points.index((inertia[x],inertia[e])) for x,e in mixed_points),
    }
    expected={'single_half':9,'nonzero_two_torsion':10,'two_independent_halves':57,'half_and_nonzero_two_torsion':49}
    data={}
    for name,g in representations.items():
        cycles=cycle_lengths(g);ramification=len(g)-len(cycles)
        genus=1-len(g)+12*ramification
        assert genus==expected[name]
        data[name]={'degree':len(g),'inertia_cycle_lengths':cycles,'ramification_per_I1':ramification,'genus':genus}
    # Castelnuovo-Severi bounds against a map of degree d to genus0/1.
    bounds={}
    for name,n in [('two_independent_halves',4),('half_and_nonzero_two_torsion',3)]:
        bounds[name]=[n*9+d+(n-1)*(d-1) for d in range(1,5)]
        assert max(bounds[name])<data[name]['genus']
    # Kadets-Vogt Theorem1.3 alternative-(2) genus bounds for d<=4.
    kv={}
    for d in range(1,5):
        m=(d+1)//2-1;eps=3*d-1-6*m
        kv[d]=max(d*(d-1)//2+1,3*m*(m-1)+m*eps)
    assert kv=={1:1,2:2,3:4,4:7}
    return {'schema':'joint-halving-incidence-v1','status':'FINITE_PREMISES_VERIFIED',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'linear_group_order':6,'single_affine_group_order':24,'joint_affine_group_order':96,
        'invariant_subspaces':[sorted(W) for W in sorted(invariant,key=lambda W:(len(W),sorted(W)))],
        'cocycles':cocycles,'joint_stabilizer_order':6,'single_half_stabilizer_order':24,
        'single_half_stabilizer_maximal':True,'curves':data,'CS_bounds_d1_to_d4_target_genus1':bounds,
        'Kadets_Vogt_1_3_genus_bounds':kv,
        'limits':{'cpu_seconds':10,'memory_bytes':512*1024**2},'elapsed_seconds':time.monotonic()-start,
        'proof_boundary':'Finite group and integer premises only. Connectedness, local geometry, curve inequalities, published low-degree-point theorems and finiteness deduction require the written proof. No exceptional point enumeration or new rank gain.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record',type=Path,required=True);args=p.parse_args()
    result=verify();args.record.parent.mkdir(parents=True,exist_ok=True)
    with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'genera':{k:v['genus'] for k,v in result['curves'].items()},'elapsed_seconds':result['elapsed_seconds']}))
