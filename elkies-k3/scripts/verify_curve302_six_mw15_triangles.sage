#!/usr/bin/env sage-python
"""Replay six MW15 pencils, six inequivalence witnesses and final302 inverse.

Build accepts a supplied rational coordinate-change witness and verifies it
exactly. Replay reads that witness from the final certificate. 120 seconds,
one worker; independent specialized-fibre checks use three values per model.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import QQ, GF, PolynomialRing, ZZ

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
PREPARE = ROOT/'elkies-k3/scripts/certify_curve302_six_mw15_triangles.sage'
PROTOCOL = ART/'elkies-k3-curve302-six-mw15-triangles-protocol-v1.json'
SCREEN = ART/'elkies-k3-curve302-six-mw15-triangles-v1.json'
CONVERT = ROOT/'elkies-k3/scripts/certify_curve302_triangle_generic_conversion.sage'
BRANCH = ROOT/'elkies-k3/scripts/certify_curve302_triangle_branch_separation.sage'
FIBRE = ROOT/'elkies-k3/scripts/probe_curve302_triangle_mw14.sage'
OUT = ART/'elkies-k3-curve302-six-mw15-triangles-complete-v1.json'

def digest(p): return sha256(p.read_bytes()).hexdigest()

def build(check):
    prep = runpy.run_path(str(PREPARE))
    protocol = prep['prepare'](); assert protocol == json.loads(PROTOCOL.read_text())
    screen = prep['probe'](); assert screen == json.loads(SCREEN.read_text())
    assert screen['excluded_count']==5 and screen['records'][4]['status']=='UNKNOWN'
    old = json.loads(OUT.read_text()) if check else None
    witness = old['row4_generic_conversion'] if check else json.loads((ROOT/'artifacts/local/elkies-k3/curve302-mw15-row4-qq.json').read_text())
    conversion = runpy.run_path(str(CONVERT)); model = protocol['models'][4]
    data = conversion['convert'](QQ,model,provided=witness['forward_plane_v_coefficients_in_old_x'])
    packet = conversion['packet'](data); assert packet==witness
    inverse = conversion['inverse'](data['j']); assert inverse['prime']==97
    # Separately replay the exact integer comparison with basic modular arithmetic.
    f = list(map(ZZ,inverse['primitive_comparison'])); p=inverse['prime']
    assert len(f)==25 and f[-1]%p and all(sum(c*pow(a,i,p) for i,c in enumerate(f))%p for a in range(p))
    inverse['maximum_coefficient_bits']=max(abs(c).nbits() for c in f)
    branch = runpy.run_path(str(BRANCH)); fibre = runpy.run_path(str(FIBRE))['fibre_j']
    compiler = runpy.run_path(str(prep['COMPILER']))
    field=GF(1013); R=PolynomialRing(field,'s'); records=[]
    for model,row in zip(protocol['models'],screen['records']):
        d=row['results'][0]; assert d['prime']==1013
        n,den=R(d['j_numerator']),R(d['j_denominator'])
        h,samples=branch['branch_polynomial'](n,den)
        assert h.degree()==44
        K,f,c,z,meet=compiler['make_pencil'](field,model)
        independent=[]
        for a in range(8):
            value=fibre(field,K,f,a)
            if value is not None:
                assert den(a) and n(a)/den(a)==field(value)
                independent.append([a,value])
            if len(independent)==3: break
        assert len(independent)==3
        records.append(dict(index=model['index'],normalized_branch_polynomial=list(map(int,h.list())),
                            binary_discriminant_samples=samples,independent_fibre_values=independent))
        print('BRANCH_AND_FIBRE_REPLAY',model['index'],flush=True)
    assert len({tuple(r['normalized_branch_polynomial']) for r in records})==6
    return dict(schema='curve302.six-mw15-triangles.complete.v1',
                status='SIX_INEQUIVALENT_MW15_TRIANGLES_EXCLUDE302',
                input_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),PREPARE,PROTOCOL,SCREEN,CONVERT,BRANCH,FIBRE]},
                row4_generic_conversion=packet,row4_exact_inverse=inverse,branch_records=records,
                boundary='Six inequivalent elliptic fibrations on one X948 surface, all excluded at every rational parameter. Abstract full generic MW lattices are proved; Weierstrass-coordinate MW bases are not constructed. This finite roster does not classify degree-three fibrations or exclude all302 parents.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--build',action='store_true'); args=parser.parse_args()
    signal.alarm(120); result=build(not args.build)
    if args.build:
        assert not OUT.exists(),'Preserve existing evidence'
        OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else: assert result==json.loads(OUT.read_text())
    print(result['status'],flush=True)
