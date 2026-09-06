#!/usr/bin/env python3
"""Coordinate witness for a derivative reciprocity constraint on the frozen panel."""
import argparse
from pathlib import Path
import sys
import retrospective as r
import fresh_strict_boundary as direct
import fresh_governing_panel as base
import matched103b2_class_boundary as prior

DIRECT_FILE=Path(direct.__file__)
PROTOCOL=Path(__file__).with_name('FRESH_STRICT_BOUNDARY_COORDINATE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_strict_boundary_coordinates_v1.json'
REPORT=r.OUT/'rank_jump_fresh_strict_boundary_coordinate_comparison_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-strict-boundary-coordinates-v1'


def bindings():
    paths=(Path(__file__),PROTOCOL,DIRECT_FILE,base.INPUT,base.OUTPUT,
           Path(base.__file__),Path(prior.__file__),Path(r.__file__),base.LOCAL)
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def worker(token):
    from sage.all import pari,AA,EllipticCurve,GF,matrix
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old,f,pts,primes,nf=prior.setup(token);loc=old['local']
    assert loc['status'] in ('PASS','PARTIAL') and 'S_finite' in loc
    theta=pari.Mod('z',pari(f));delta=f.discriminant();beta=-pari(delta)*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==delta**4
    E=EllipticCurve([0,0,0,f[1],f[0]])
    S=[2]+[p for p in primes if p!=2 and E.local_data(p,proof=True).conductor_valuation()>0]
    assert S==loc['S_finite'];omitted=[]
    for p in sorted(set(primes)-set(S)):
        assert p!=2
        vals=[int(pari.idealval(nf,beta,P)) for P in pari.idealprimedec(nf,p)]
        assert all(v%2==0 for v in vals);omitted.append({'prime':p,'valuations':vals})
    roots=f.roots(AA,multiplicities=False);signs=[int(-delta*f.derivative()(a)<0) for a in roots]
    real_dim=int(len(roots)==3);real_basis=[[0,1,1]] if real_dim else []
    assert (delta>0)==bool(real_dim)
    witness={'place':'infinity','type':'real Hilbert dot product'} if real_dim else None
    if real_dim:assert sum(a*b for a,b in zip(signs,real_basis[0]))%2==1
    joint=[0]*len(pts);width=0;ell=real_dim;rows=[]
    for p in S:
        saved=next(x for x in loc['local'] if x['place']==p)
        dim=len(pari.idealprimedec(nf,p))-1+int(p==2)
        assert dim==saved['point_kummer_dimension']
        sigs=list(map(r.pack,saved['signatures']));rank=r.rank(sigs);assert rank<=dim
        for i,s in enumerate(sigs):joint[i]|=s<<width
        width+=len(saved['signatures'][0]);ell+=dim
        record={'place':p,'point_dimension':dim,'generic_local_rank':rank,'tested_derivative':False}
        if witness is None and rank==dim:
            chars=LocalSquareclasses(nf,p)
            # Match the character convention to the retained generic signatures.
            assert [list(chars.signature(pari(x)-theta)) for x,y in pts]==saved['signatures']
            ds=r.pack(chars.signature(beta));outside=r.rank(sigs+[ds])>rank
            record.update(tested_derivative=True,derivative_signature=ds,outside_full_point_image=outside)
            if outside:witness={'place':p,'type':'outside full generic local point span'}
        rows.append(record)
    real_sigs=[[int(x<a) for a in roots] for x,y in pts]
    assert all(s in ([0,0,0],[0,1,1]) if real_dim else s==[0] for s in real_sigs)
    for i,s in enumerate(real_sigs):joint[i]|=r.pack(s)<<width
    g=int(matrix(GF(2),[[v>>j&1 for j in range(width+len(roots))] for v in joint]).rank())
    m=len(pts);k=m-g;assert k==loc['strict_generic_dimension'];b=int(witness is not None)
    upper=ell-b;assert upper>=g
    return {'status':'PASS','S_finite':S,'polynomial_discriminant':str(delta),
       'derivative_coefficients':[str(pari.lift(beta).polcoef(i)) for i in range(3)],
       'derivative_norm':str(delta**4),'omitted_good_prime_valuations':omitted,'local':rows,
       'real_derivative_signs':signs,'complete_real_point_basis':real_basis,'witness':witness,
       'generic_dimension':m,'generic_local_dimension':g,'generic_strict_dimension':k,
       'local_point_product_dimension':ell,'reciprocity_constraint_rank_lower_bound':b,
       'Selmer_boundary_dimension_interval':[g,upper],'additional_boundary_capacity_upper_bound':upper-g,
       'additional_strict_dimension':'UNKNOWN','additional_quotient_CT':'UNKNOWN',
       'boundary':'A witnessed equation-defined constraint, or a conservative unreduced boundary. No new class representative or point supplied.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','report']);p.add_argument('--token');args=p.parse_args()
    # Reuse orchestration and label joining unchanged, with this frozen worker and paths.
    direct.PROTOCOL=PROTOCOL;direct.OUTPUT=OUTPUT;direct.REPORT=REPORT;direct.WORK=WORK
    direct.bindings=bindings;direct.__file__=__file__
    if args.mode=='capture':direct.capture()
    elif args.mode=='report':direct.report()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        r.write_new(WORK/f'{args.token}.json',{'bindings':bindings(),**worker(args.token)})
