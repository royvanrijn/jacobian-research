#!/usr/bin/env sage-python
"""One bounded exact conic solve for a smaller, equation-derived M18 seed.

This constructive lane does not alter or refill the ten-million-address funnel.
No quartic point search, public point, rank label or catalogue is read.
"""
import argparse
from fractions import Fraction as F
import json
from math import gcd,isqrt
from pathlib import Path
import sys

from sage.all import ZZ,matrix,pari
import det1092_funnel as f
from det1092_funnel_worker import prepare,conic_points,first_seed,seal_seed,amplifier_context
from v3_warm_support import atomic,bindings,read,require,sha


def main(parent_run,output):
    parent=read(parent_run/'protocol.json');bindings(f.ROOT,parent['inputs']);bindings(f.ROOT,parent['sources'])
    arithmetic=f.Arithmetic(parent)
    common=gcd(*arithmetic.conic);root=isqrt(common)
    require(root*root==common,'conic content is not a square')
    q0,q1,q2=[v//common for v in arithmetic.conic]
    G=matrix(ZZ,[[2*q2,q1,0],[q1,2*q0,0],[0,0,-2]])
    p=dict(parent)
    p.update(schema='det1092-small-conic-seed.v1',profile='constructive_conic',
             population=0,maximum_draws=0,domain='det1092-small-conic-seed-v1',
             sources=f.source_bindings(),
             scope='One equation-only conic solve, exact rank18 confirmation, and optional unchanged V3. Separate from broad-population yield.',
             constructive=dict(method='PARI qfsolve on the primitive ternary form',
                               primitive_conic=[str(v) for v in [q0,q1,q2]],square_content=str(common),
                               gram=[list(map(str,r)) for r in G.rows()],
                               maximum_solutions=1,constructor_wall_seconds=90))
    # f.protocol's Python metadata refers to the regular-Python controller;
    # checked_protocol validates this Sage worker against worker_software.
    atomic(output/'protocol.json',p,immutable=True)
    proof_path=output/'conic-solve.json'
    if proof_path.exists():
        solve=read(proof_path);v=matrix(ZZ,3,1,list(map(ZZ,solve['vector'])))
    else:
        raw=pari.qfsolve(G)
        require(raw.type()=='t_COL','qfsolve did not return a rational point')
        v=matrix(ZZ,3,1,[ZZ(x) for x in raw])
        require(v[1,0]!=0 and (v.transpose()*G*v)[0,0]==0,'bad isotropic vector')
        atomic(proof_path,dict(status='EXACT_RATIONAL_CONIC_POINT',vector=list(map(str,v.column(0))),
                              gram=[list(map(str,r)) for r in G.rows()],
                              protocol_sha256=sha(output/'protocol.json')),immutable=True)
    require(v[1,0]!=0 and (v.transpose()*G*v)[0,0]==0,'retained conic identity failed')
    s=F(int(v[0,0]),int(v[1,0]))
    row=arithmetic.candidate(arithmetic.record(0,s.numerator,s.denominator),'separate_constructive_conic')
    row['id']='conic-small-01'
    require(row['conic_splitting']=='SPLIT','primitive and original conic splitting differ')
    atomic(output/'selection.json',dict(status='SEALED_CONSTRUCTIVE_SINGLETON',seed_inputs=[row],
                                       arithmetic_candidates=[row],parameters=1,
                                       protocol_digest=f.digest(f.packed(p))),immutable=True)
    folder=output/'seeds'/row['id'];prepared=prepare(folder,row)
    if prepared is None:
        print('UNRESOLVED_INHERITED_RANK',flush=True);return
    model,base,proof=prepared;points=conic_points(row)
    atomic(folder/'conic-points.json',[list(map(str,P)) for P in points],immutable=True)
    found=first_seed(model,base,proof,points)
    if found is None:
        atomic(folder/'result.json',dict(status='CONIC_POINT_INDEPENDENCE_UNRESOLVED',parameter=row['parameter']),immutable=True)
        print('UNRESOLVED_EXTRA_DIRECTION',flush=True);return
    seal_seed(folder,row,model,base,*found,evidence=dict(kind='primitive_conic_solve',sha256=sha(proof_path)))
    # Exact V3 map only; any actual cascade has a separate supervision record.
    from v3_warm_engine import preflight
    ctx=amplifier_context(output,row['id']);preflight(ctx)
    result=dict(status='CERTIFIED_SMALL_CONIC_M18',parameter=row['parameter'],
                original_parameter=row['original_parameter'],
                j_numerator_bits=row['j_numerator_bits'],j_denominator_bits=row['j_denominator_bits'],
                rank_lower_bound=18,seed_sha256=sha(folder/'m18.json'),
                counted_as_broad_population_yield=False,
                boundary='Explicit rank18 subgroup on a smaller constructive fibre. No higher rank, conductor, novelty or amplification outcome claimed.')
    atomic(output/'constructive-result.json',result,immutable=True)
    print(json.dumps(result,sort_keys=True),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parent-run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();main(a.parent_run.resolve(),a.output.resolve())
