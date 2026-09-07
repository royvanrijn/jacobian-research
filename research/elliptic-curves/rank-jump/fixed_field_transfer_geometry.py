#!/usr/bin/env python3
"""Equation-only obstruction to fixing the descent field by low-genus base change."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_symbolic_discriminant as prior
import fresh_governing_panel as panel

PROTOCOL=Path(__file__).with_name('FIXED_FIELD_TRANSFER_GEOMETRY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_fixed_field_transfer_geometry_inputs_v1.json'
PROVENANCE=r.OUT/'rank_jump_fixed_field_transfer_geometry_provenance_v1.json'
OUTPUT=r.OUT/'rank_jump_fixed_field_transfer_geometry_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fixed-field-transfer-geometry-v1'


def export():
    specs=r.read(PROTOCOL)['families'];rows=[];sources={}
    parameters={'074d9':'2818/1535','103b2':'3726/881','11952':'-2448/11',
                'a1-fibration-01':'-1867/270','a1-fibration-05':'3/17'}
    # Frozen parameter provenance only; no exceptional coordinates are projected.
    manifest=r.read(panel.MANIFEST)
    for fam,t in parameters.items():
        if fam!='a1-fibration-05':assert any(x['family']==fam and x['parameter']==t for x in manifest['rows'])
    ref=r.read(r.OUT/'rank_jump_bounded_gain_reference_provenance_v1.json');assert ref['parameter']=='3/17'
    for p in prior.ATLASES:
        sources[str(p.relative_to(r.ROOT))]=r.digest(p.read_bytes())
        for f in r.read(p)['families']:
            family=f.get('family',f.get('fibration_id'))
            if family in specs:rows.append({'family':family,'A':f['A_coefficients_low_to_high'],
                'B':f['B_coefficients_low_to_high'],'irreducibility_witness_parameter':parameters[family]})
    assert sorted(x['family'] for x in rows)==sorted(specs)
    r.write_new(INPUT,{'schema':'rank-jump.fixed-field-transfer-geometry-inputs.v1','families':rows})
    for p in (Path(__file__),PROTOCOL,INPUT,panel.MANIFEST,r.OUT/'rank_jump_bounded_gain_reference_provenance_v1.json'):
        sources[str(p.relative_to(r.ROOT))]=r.digest(p.read_bytes())
    r.write_new(PROVENANCE,{'schema':'rank-jump.fixed-field-transfer-geometry-provenance.v1','bindings':sources,
        'whitelist':'Compact A,B coefficient arrays and already-retained rational parameters only. No points, classes or rank labels.'})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,Path(r.__file__))}


def worker(family):
    from sage.all import QQ,GF,PolynomialRing,prime_range
    data=next(x for x in r.read(INPUT)['families'] if x['family']==family)
    R=PolynomialRing(QQ,'t');A=R(data['A']);B=R(data['B']);D=-16*(4*A**3+27*B**2)
    assert A.degree()<=8 and B.degree()<=12 and D.degree()==24
    assert D.gcd(A)==1 # All finite singular fibres are multiplicative.
    fac=D.squarefree_decomposition();product=R(fac.unit());factors=[];odd=0
    for q,e in fac:
        assert q.gcd(q.derivative())==1
        for v in factors:assert q.gcd(R(v['coefficients_ascending']))==1
        product*=q**e
        factors.append({'coefficients_ascending':list(map(str,q.list())),'multiplicity':int(e),'degree':int(q.degree())})
        if e%2:odd+=int(q.degree())
    assert product==D and odd>0 and odd%2==0
    # Good infinity in the degree(8,12) compact model.
    infA=A[8];infB=B[12];infD=-16*(4*infA**3+27*infB**2)
    assert infD==D[24] and infD
    t0=QQ(data['irreducibility_witness_parameter']);a,b=A(t0),B(t0);assert D(t0)
    witness=None
    for p in prime_range(3,r.read(PROTOCOL)['limits']['modular_prime_bound']+1):
        if a.denominator()%p==0 or b.denominator()%p==0:continue
        ap,bp=GF(p)(a),GF(p)(b)
        if 4*ap**3+27*bp**2==0:continue
        roots=[int(x) for x in GF(p) if x**3+ap*x+bp==0]
        if not roots:
            witness={'parameter':str(t0),'prime':int(p),'A_mod_p':int(ap),'B_mod_p':int(bp),'roots_mod_p':roots};break
    assert witness is not None
    return {'status':'PASS','family':family,'discriminant_ascending':list(map(str,D.list())),
        'discriminant_unit':str(fac.unit()),'squarefree_factors':factors,
        'c4_discriminant_gcd':'1','infinity_coefficients':list(map(str,(infA,infB,infD))),
        'cubic_irreducibility_witness':witness,'geometric_monodromy':'S3',
        'geometric_transposition_branch_points':odd,'root_cover_degree':3,'root_cover_genus':(odd-4)//2,
        'quadratic_resolvent_degree':2,'quadratic_resolvent_genus':(odd-2)//2,
        'constant_cubic_isomorphism_cover_degree':6,'constant_cubic_isomorphism_cover_genus':(3*odd-10)//2,
        'constant_field_basechange_degree_lower_bound':6,'constant_field_basechange_genus_lower_bound':(3*odd-10)//2,
        'rational_or_genus_one_constant_field_basechange_exists':False,
        'boundary':'Geometric base-change obstruction. No arithmetic rank, Selmer dimension, local solubility or rational-point count computed.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for family in r.read(PROTOCOL)['families']:
        path=WORK/f'{family}.json'
        if not path.exists():
            error=None
            with (WORK/f'{family}.log').open('x') as log:
                try:
                    p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--family',family],stdout=log,stderr=log,timeout=30)
                    if p.returncode:error='Worker failure'
                except subprocess.TimeoutExpired:error='Bounded worker timeout'
            if error:r.write_new(path,{'family':family,'status':'UNKNOWN','reason':error,'bindings':bindings()})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(family,row['status'],row.get('geometric_transposition_branch_points'),row.get('constant_cubic_isomorphism_cover_genus'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fixed-field-transfer-geometry.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':
        from sage.all import pari
        pari.allocatemem(64000000,268435456,silent=True)
        r.write_new(WORK/f'{args.family}.json',{'bindings':bindings(),**worker(args.family)})
    else:globals()[args.mode]()
