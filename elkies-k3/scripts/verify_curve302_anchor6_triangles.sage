#!/usr/bin/env sage-python
"""Replay the complete P6 height-six rank gate and all high-rank inverses.

One worker, 600 seconds. No parameter-height sweep. The two supplied QQ
coordinate witnesses are verified algebraically; all88 finite-field j-map
screens and all880 lattice cases are reconstructed. One independent
specialized-fibre normalization is checked for each of the88 MW14 pencils.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
import sys
from sage.all import QQ,ZZ,GF,PolynomialRing

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results'
PREPARE=ROOT/'elkies-k3/scripts/certify_curve302_anchor6_mw14.sage'
GATE=ROOT/'elkies-k3/scripts/certify_curve302_anchor6_triangle_gate.sage'
PROTOCOL=ART/'elkies-k3-curve302-anchor6-mw14-protocol-v1.json'
SCREEN=ART/'elkies-k3-curve302-anchor6-mw14-v1.json'
GATE_OUT=ART/'elkies-k3-curve302-anchor6-triangle-gate-v1.json'
MW15=ART/'elkies-k3-curve302-six-mw15-triangles-complete-v1.json'
CONVERT=ROOT/'elkies-k3/scripts/certify_curve302_triangle_generic_conversion.sage'
BRANCH=ROOT/'elkies-k3/scripts/certify_curve302_triangle_branch_separation.sage'
FIBRE=ROOT/'elkies-k3/scripts/probe_curve302_triangle_mw14.sage'
OUT=ART/'elkies-k3-curve302-anchor6-triangles-complete-v1.json'

def digest(p): return sha256(p.read_bytes()).hexdigest()

def build(check):
    prep=runpy.run_path(str(PREPARE))
    protocol=prep['prepare'](); assert protocol==json.loads(PROTOCOL.read_text())
    screen=prep['probe'](); assert screen==json.loads(SCREEN.read_text())
    assert screen['excluded_count']==86
    assert [r['index'] for r in screen['records'] if r['status']=='UNKNOWN']==[4,33]
    gate=runpy.run_path(str(GATE))['build'](); assert gate==json.loads(GATE_OUT.read_text())
    prior=json.loads(MW15.read_text())
    assert prior['status']=='SIX_INEQUIVALENT_MW15_TRIANGLES_EXCLUDE302'
    for p,h in prior['input_sha256'].items(): assert digest(ROOT/p)==h,p
    old=json.loads(OUT.read_text()) if check else None
    convert=runpy.run_path(str(CONVERT)); rational=[]
    for index in [4,33]:
        witness=next(r['generic_conversion'] for r in old['rational_inverses'] if r['index']==index) if check else json.loads((ROOT/f'artifacts/local/elkies-k3/curve302-anchor6-row{index}-qq.json').read_text())
        data=convert['convert'](QQ,protocol['models'][index],provided=witness['forward_plane_v_coefficients_in_old_x'])
        packet=convert['packet'](data); assert packet==witness
        inverse=convert['inverse'](data['j']); assert inverse['prime']==61
        f=list(map(ZZ,inverse['primitive_comparison'])); p=inverse['prime']
        assert len(f)==25 and f[-1]%p
        assert all(sum(c*pow(a,i,p) for i,c in enumerate(f))%p for a in range(p))
        rational.append(dict(index=index,generic_conversion=packet,inverse=inverse))
        print('EXACT_QQ_INVERSE',index,p,flush=True)
    branch=runpy.run_path(str(BRANCH)); fibre=runpy.run_path(str(FIBRE))['fibre_j']
    compiler=runpy.run_path(str(prep['COMPILER'])); field=GF(1013); R=PolynomialRing(field,'s')
    branches=[]
    for model,row in zip(protocol['models'],screen['records']):
        d=row['results'][0]; assert d['prime']==1013 and 'j_numerator' in d
        n,den=R(d['j_numerator']),R(d['j_denominator'])
        h,samples=branch['branch_polynomial'](n,den); assert h.degree()==43
        K,f,c,z,meet=compiler['make_pencil'](field,model); independent=None
        for a in range(8):
            value=fibre(field,K,f,a)
            if value is not None:
                assert den(a) and n(a)/den(a)==field(value)
                independent=[a,value]; break
        assert independent is not None
        branches.append(dict(index=model['index'],normalized_branch_polynomial=list(map(int,h.list())),
                             binary_discriminant_samples=samples,independent_fibre_value=independent))
        print('BRANCH_AND_FIBRE',model['index'],flush=True)
    assert len({tuple(r['normalized_branch_polynomial']) for r in branches})==88
    return dict(schema='curve302.anchor6-triangles.complete.v1',status='NO_MW14_OR_HIGHER_302_PARENT_IN_COMPLETE_P6_TRIANGLE_CLASS',
                input_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),PREPARE,GATE,PROTOCOL,SCREEN,GATE_OUT,MW15,CONVERT,BRANCH,FIBRE]},
                rational_inverses=rational,branch_records=branches,
                triangle_count=880,rank_at_most13_count=786,excluded_MW14_count=88,excluded_MW15_count=6,
                boundary='Complete only for triangles O+P6+Q with height(Q)=6 and pairing(P6,Q)=3 on the pinned 11952 model. The low-rank786 are not excluded as fibres of302. Other anchors, nontriangle degree-three fibrations and other surfaces remain open. No302 parent or full Weierstrass-coordinate MW basis is constructed.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--build',action='store_true'); args=parser.parse_args()
    sys.set_int_max_str_digits(20000)
    signal.alarm(600); result=build(not args.build)
    if args.build:
        assert not OUT.exists(); OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    assert result==json.loads(OUT.read_text())
    print(result['status'],flush=True)
