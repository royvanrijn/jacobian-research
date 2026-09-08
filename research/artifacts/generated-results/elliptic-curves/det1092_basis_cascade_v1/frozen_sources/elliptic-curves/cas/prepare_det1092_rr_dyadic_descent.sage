#!/usr/bin/env sage-python
"""Bounded real/dyadic preflight: exact input, retained approximations.

A p-adic approximation is not promoted to a local-image certificate here.
"""
import argparse,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,pari,prod,floor
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_preflight_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def prepare(i):
    paths=[ART/'det1092_rr_good_local_images_v1'/('case-%02d.json'%j) for j in range(10)]
    protocol={'classification':'exact local arithmetic preflight, not a full Selmer group',
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
              'limits':{'cases':10,'seconds_per_case':20,'p':2,'maximum_absolute_precision':512,
                        'factorpadic_calls_per_case':1,'point_searches':0,'class_groups':0,
                        'parameter_searches':0,'remote_arithmetic':0,'pilot_changes':0},
              'precision_policy':'N=max(64,2*v2(discriminant(integral_monic_model))+16), capped at512; exceeding cap is unresolved.',
              'script_sha256':sha(Path(__file__))}
    OUT.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    source=json.loads(paths[i].read_text());R=PolynomialRing(QQ,'Z');Z=R.gen();q=R(source['q'])
    f=q.monic();s=min([ZZ(0)]+[floor(a.valuation(2)/(6-j)) for j,a in enumerate(f.list()[:-1]) if a])
    f=R(f(2**s*Z)/2**(6*s));assert f.is_monic() and all(a.valuation(2)>=0 for a in f if a)
    vd=f.discriminant().valuation(2);N=max(64,2*int(vd)+16)
    prefix={'classification':'verified exact normalization; p-adic factor output requires independent certification',
            'case_index':i,'source':str(paths[i].relative_to(ROOT)),'source_sha256':sha(paths[i]),
            'q':source['q'],'T_equals_two_power_times_Z':int(s),'integral_monic_model':list(map(str,f.list())),
            'model_discriminant_v2':int(vd),'requested_absolute_precision':N,
            'real_root_count_by_sturm':int(pari(q).polsturm()),
            'limits':protocol['limits'],'software':{'sage':sage_version,'pari':str(pari.version())}}
    retain(OUT/('case-%02d-input.json'%i),prefix)
    if N>512:
        result={**prefix,'status':'UNRESOLVED_PRECISION_CAP'}
    else:
        raw=pari(f).factorpadic(2,N);factors=[];multiplicities=[]
        for j in range(raw.nrows()):
            h=raw[j,0];factors.append(R([QQ(h.polcoef(k).lift()) for k in range(int(h.poldegree())+1)]))
            multiplicities.append(int(raw[j,1]))
        assert all(e==1 for e in multiplicities) and all(h.is_monic() for h in factors)
        error=f-prod(factors);ve=min([a.valuation(2) for a in error if a]+[ZZ(N)])
        assert ve>=N
        result={**prefix,'status':'PASS_DYADIC_FACTORIZATION_APPROXIMATION_NOT_YET_CERTIFIED',
                'factor_degrees':[int(h.degree()) for h in factors],
                'factors':[list(map(str,h.list())) for h in factors],
                'factor_coefficient_valuations':[[int(a.valuation(2)) if a else None for a in h] for h in factors],
                'factor_product_error_v2_at_least':int(ve),
                'raw_factorpadic':str(raw),
                'full_local_Kummer_image':'NOT_COMPUTED','full_global_Selmer_group':'NOT_COMPUTED'}
    retain(OUT/('case-%02d.json'%i),result)
    print('case',i,result['status'],'degrees',result.get('factor_degrees'),'real',prefix['real_root_count_by_sturm'],
          'vdisc',vd,'precision',N,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(20);prepare(args.case)
