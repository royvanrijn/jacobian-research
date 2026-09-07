#!/usr/bin/env python3
"""Retrospective CT blocks adapted to the marked generic strict subgroup."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import scalar_cup
import strict_class_blocks as strict
import production_minus_twist_completion as completed
import verify_production_minus_twist as verified

OUTPUT=r.OUT/'rank_jump_production_twist_blocks_v1.json'


def pairing(M,x,y):
    return sum((row&y).bit_count() for i,row in enumerate(M) if x>>i&1)%2


def decompose(M,first_pair=None):
    """Symplectic Gram--Schmidt, optionally preserving one designated pair."""
    remaining=[1<<i for i in range(len(M))];pairs=[]
    while True:
        pair=first_pair if not pairs and first_pair else next(
            ((x,y) for x in remaining for y in remaining if pairing(M,x,y)),None)
        if pair is None:break
        x,y=pair;assert pairing(M,x,y)==1;pairs.append([x,y])
        projected=[v ^ (x if pairing(M,v,y) else 0) ^ (y if pairing(M,v,x) else 0) for v in remaining]
        remaining=list(r.basis(projected).values())
    flat=[x for pair in pairs for x in pair]+remaining
    assert r.rank(flat)==len(M)
    for i,x in enumerate(flat):
        for j,y in enumerate(flat):
            assert pairing(M,x,y)==int(i//2==j//2 and i!=j and i<2*len(pairs) and j<2*len(pairs))
    return pairs,remaining


def build(check=False):
    rows=[]
    for cup in r.read(scalar_cup.OUTPUT)['production_cases']:
        index=cup['case_index'];old=r.read(strict.OUTPUT)['rows'][index]
        space=old['witness_strict_kernel_masks'];generic=old['generic_strict_kernel_masks']
        G=[lc.coordinates(x,space) for x in generic];M=list(map(r.pack,cup['scalar_cup_matrix']))
        first_pair=(G[0],next(1<<j for j in range(len(M)) if pairing(M,G[0],1<<j)))
        if len(G)==2:assert pairing(M,*G)==1;first_pair=tuple(G)
        pairs,radical=decompose(M,first_pair)
        m=old['generic_dimension']
        def encode(v):
            witness=lc.lift(v,space)
            return {'strict_basis_mask':v,'witness_basis_mask':witness,'exceptional_quotient_mask':witness>>m}
        rank_generic=r.rank([r.pack(pairing(M,x,y) for y in G) for x in G])
        rank_cross=r.rank([r.pack(pairing(M,x,1<<j) for j in range(len(M))) for x in G])
        tail=[v for pair in pairs[1:] for v in pair]+radical
        rows.append({'case_index':index,'id':cup['id'],'generic_rank':m,
            'known_independent_rank':old['witness_dimension'],
            'observed_quotient_rank':old['witness_dimension']-m,
            'generic_strict_dimension':len(G),'retained_strict_dimension':len(space),
            'strict_quotient_dimension':len(space)-len(G),
            'generic_restricted_CT_rank':rank_generic,'generic_cross_CT_rank':rank_cross,
            'full_retained_CT_rank':2*len(pairs),
            'hyperbolic_pairs':[[encode(x),encode(y)] for x,y in pairs],
            'radical_basis':[encode(v) for v in radical],
            'after_first_pair_dimension':len(tail),'after_first_pair_CT_rank':2*(len(pairs)-1),
            'after_first_pair_exceptional_quotient_dimension':r.rank([lc.lift(v,space)>>m for v in tail]),
            'classification':'solubility: CT form on the transported retained strict space of the -1 twist',
            'boundary':'Basis decomposition is not a canonical partition of original rational points or a proof that radical covers are rational. All oracle coordinates are retrospective.'})
    bindings={str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),scalar_cup.OUTPUT,strict.OUTPUT,completed.OUTPUT,verified.OUTPUT)}
    result={'schema':'rank-jump.production-twist-blocks.v1','bindings':bindings,'rows':rows}
    if check:assert r.read(OUTPUT)==result;print('PASS generic-adapted production CT blocks')
    else:r.write_new(OUTPUT,result)
    for row in rows:print(row['id'],row['strict_quotient_dimension'],row['full_retained_CT_rank'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();build(args.mode=='check')
