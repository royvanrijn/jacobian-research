#!/usr/bin/env sage-python
"""Exact coherent seed intake and prospective geometry for twelve fixed fibres."""
import argparse
import sys
from pathlib import Path
from dataclasses import asdict
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,EllipticCurve,matrix,pari
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import mestre_parent_calibration as control
import certify_compact_r17_candidates as cert
from mestre_parent_adapter import specialize
from research_runtime.store import checkpoint,digest
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from memory_rank_certificate import checked_rank

def intake(index):
    p=control.intake_protocol();row=p['rows'][index];out=control.BATCH/row['id']/'seed.json'
    if out.exists():raise FileExistsError('preserve own11 seed intake')
    checkpoint(out,{'status':'RUNNING','intake_protocol_hash':digest(p),'row':row})
    model,images=specialize(row['outer_u'],row['fibre_T']);E=EllipticCurve(QQ,[QQ(str(c)) for c in model])
    P=[E([QQ(str(z)) for z in point]) for point in images];cloud=[images[0]]
    for Q in P[1:]:
        halves=(Q-P[0]).division_points(2)
        if len(halves)!=1 or halves[0].is_zero() or 2*halves[0]!=Q-P[0]:raise ArithmeticError('unique rational divisor half required')
        cloud.append(tuple(cert.F(str(z)) for z in halves[0].xy()))
    cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];torsion=None
    for prime in _primes_up_to(997):
        if prime==2:continue
        if torsion is None and cert.short_curve_has_no_rational_2_torsion_modular_certificate(model,prime):torsion=prime
        try:sig=signature(cache,model,cloud,prime)
        except ValueError:continue
        before=len(pivots)
        for bits in sig.rows:insert(pivots,bits)
        if len(pivots)>before:signatures.append(asdict(sig))
    if len(pivots)!=11 or torsion is None:raise ArithmeticError('fixed coherent11 seed gate failed; no refill')
    chosen=sorted(pivots);basis=[cloud[i] for i in chosen]
    proof=checked_rank(model,basis,[s['prime'] for s in signatures],torsion)
    checkpoint(out,{'status':'PASS','intake_protocol_hash':digest(p),'family':row['family'],
      'outer_u':row['outer_u'],'parameter':row['fibre_T'],'fibre_T':row['fibre_T'],
      'curve':list(map(str,model)),'points':[list(map(str,P)) for P in basis],
      'generic_points':[list(map(str,P)) for P in basis],
      'covariant_images':[list(map(str,P)) for P in images],
      'divisor_cloud':[list(map(str,P)) for P in cloud],
      'independent_column_indices':chosen,'rank_certificate':proof,'rank_lower_bound':11,
      'model_coefficient_bits':max(max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for c in model),
      'scope':'Exactly checked coherent covariant images, rational halves of their differences and finite-quotient independence. The11 seed points specialize rational divisor/covariant sections. No exact-rank claim.'})
    print(row['id'],'PASS11','bits',max(abs(c.numerator).bit_length() for c in model),flush=True)

def maps(index):
    control.configure(index);p=control.protocol();out=control.D/'maps.json';rank=11
    if out.exists():raise FileExistsError('preserve prospective maps')
    geometry=SourceFileLoader('mestre_geometry',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
    mapper=SourceFileLoader('mestre_mapper',str(CAS/'factor_free_pari_mapping.sage')).load_module()
    seed=cert.read(control.SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    gram,asym=geometry.canonical_height_gram(model,points);g=matrix(ZZ,geometry.rounded_gram(gram,1000000))
    u=matrix(ZZ,pari(g).qflllgram()).transpose();inverse=u.inverse()
    if abs(u.det())!=1:raise ArithmeticError('unimodular metric change required')
    reduced=u*g*u.transpose();oracle=geometry.CosetOracle(reduced.rows());sample=[]
    data={'status':'RUNNING_SAMPLE','protocol_hash':digest(p),'metric_gram':[[str(v) for v in r] for r in gram],
      'maximum_gram_asymmetry':str(asym),'rounded_gram':[list(map(int,r)) for r in g.rows()],
      'change_of_basis':[list(map(int,r)) for r in u.rows()],'reduced_gram':[list(map(int,r)) for r in reduced.rows()],
      'sample':sample,'rows':[]};checkpoint(out,data)
    for i,mask in enumerate(control.masks(p)):
        residue=matrix(ZZ,1,rank,[(mask>>j)&1 for j in range(rank)]);target=[int(v)%2 for v in (residue*inverse).row(0)]
        norm,rep,error=oracle.solve(target);word=list(map(int,(matrix(ZZ,1,rank,rep)*u).row(0)))
        if any((word[j]-(mask>>j))%2 for j in range(rank)) or sum(word[j]*g[j,k]*word[k] for j in range(rank) for k in range(rank))!=norm:raise ArithmeticError('exact parity/norm transport differs')
        sample.append({'parity':mask,'representative':word,'metric_norm':norm,'cvp_error':error,'reduced_representative':list(rep)})
        if (i+1)%511==0:checkpoint(out,data);print('PARITIES',i+1,'of2047',flush=True)
    data['centres']=sorted(sample,key=lambda c:(-c['metric_norm'],c['parity']))[:49];data['status']='RUNNING_MAPS';checkpoint(out,data)
    pari.allocatemem(256000000,silent=True)
    for c in data['centres']:
        data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data)
    print(control.ROW['id'],'FROZEN49 MAPS FROM2047 PARITIES',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['intake','maps']);parser.add_argument('--index',type=int,required=True);a=parser.parse_args();globals()[a.stage](a.index)
