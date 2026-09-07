#!/usr/bin/env sage-python
"""Exact specialization and frozen generic-centre PARI maps for fresh finalists."""
import argparse,sys
from pathlib import Path
from math import lcm,isqrt
from importlib.machinery import SourceFileLoader
from sage.all import QQ,EllipticCurve,pari
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from half_lattice_pointed_sieve import linear_combination
from search_observability import transform,multiply
import higher24_r17_pari_batch as batch
geometry=SourceFileLoader('fresh_r17_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()

def mapping(model,points,c):
    x,y=linear_combination(model,points,c['representative']);raw=(-3*x*x-4*model[3],-8*y,-6*x,cert.F(0),cert.F(1));den=lcm(*(v.denominator for v in raw));polynomial='+'.join(f'({int(v*den*den)})*x^{i}' for i,v in enumerate(raw))
    result=pari('my(m1,m2,C1,C2);C1=hyperellminimalmodel(['+polynomial+',0],&m1);C2=hyperellred(C1,&m2);[C2,m1,m2]')
    matrices=[tuple(cert.F(str(m[i,j])) for i in range(2) for j in range(2)) for m in (result[1][1],result[2][1])];M=multiply(*matrices)
    P=[cert.F(str(result[0][0].polcoef(i))) for i in range(5)];Q=[cert.F(str(result[0][1].polcoef(i))) for i in range(3)];disc=[4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)];transformed=transform(raw,M);k=next(i for i,v in enumerate(disc) if v);ratio=transformed[k]/disc[k]
    if ratio<=0 or any(a!=ratio*b for a,b in zip(transformed,disc)) or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator:raise ArithmeticError('PARI quartic map identity failed')
    return {'centre':c,'raw_coefficients':list(map(str,raw)),'denominator_clearing':str(den),'first_matrix':list(map(str,matrices[0])),'second_matrix':list(map(str,matrices[1])),'matrix':list(map(str,M)),'reduced_P':list(map(str,P)),'reduced_Q':list(map(str,Q)),'discriminant_quartic':list(map(str,disc)),'square_ratio':str(ratio),'coordinate_policy':{'kind':'raw','matrix':list(map(str,M))}}

def prepare(index):
    p=batch.protocol();row=p['rows'][index];folder=batch.D/row['id'];out=folder/'maps.json'
    if out.exists():raise FileExistsError('preserve exact map attempt')
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']==row['family']);original,generic=spec.specialize(f,row['parameter']);curve=EllipticCurve(QQ,[QQ(str(v)) for v in original]);minimal=curve.local_data(2).minimal_model();iso=curve.isomorphism_to(minimal);u=cert.F(str(iso.u)) if hasattr(iso,'u') else cert.F(str(iso.tuple()[0]));a1,a2,a3,a4,a6=minimal.a_invariants();c2,c4,c6=a2+a1*a1/4,a4+a1*a3/2,a6+a3*a3/4;model=tuple(cert.F(str(v)) for v in (0,0,0,c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27));points=[]
    for x,y in generic:
        X,Y=iso(curve([QQ(str(x)),QQ(str(y))])).xy();points.append((cert.F(str(X+c2/3)),cert.F(str(Y+(a1*X+a3)/2))))
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or tuple(points)!=tuple((x/u**2,y/u**3) for x,y in generic):raise ArithmeticError('exact family scale failed')
    gram,asymmetry=geometry.canonical_height_gram(model,points);rounded=geometry.rounded_gram(gram,1000000);oracle=geometry.CosetOracle(rounded);centres=[]
    for mask in p['generic_masks'][row['family']]:
        norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(17)))
        if len(rep)!=17 or any((rep[j]-(mask>>j))%2 for j in range(17)):raise ArithmeticError('specialized parity mismatch')
        centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
    centres.sort(key=lambda c:(-c['metric_norm'],c['mask']));data={'status':'RUNNING','protocol_hash':digest(p),'family':row['family'],'parameter':row['parameter'],'curve':list(map(str,model)),'generic_points':[list(map(str,P)) for P in points],'family_to_curve_scale_u':str(u),'metric_gram':[[str(v) for v in r] for r in gram],'maximum_gram_asymmetry':str(asymmetry),'centres':centres,'rows':[]};checkpoint(out,data);pari.allocatemem(256000000,silent=True)
    for c in centres:data['rows'].append(mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('HIGHER24 R17 MAPS',row['id'],len(centres),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);prepare(p.parse_args().index)
