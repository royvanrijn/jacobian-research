#!/usr/bin/env python3
"""Family discriminant factor gate on five frozen incomplete fibres."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base
import fresh_retained_factors as retained

PROTOCOL=Path(__file__).with_name('FRESH_SYMBOLIC_DISCRIMINANT_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_fresh_symbolic_discriminant_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_fresh_symbolic_discriminant_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-symbolic-discriminant-v1'
ATLASES=[r.OUT/'compact_six_r17_atlas_v1.json',r.OUT/'compact_five_mw16_atlas_v1.json']


def export():
    families=[];spec=r.read(PROTOCOL);sources={}
    for path in ATLASES:
        sources[str(path.relative_to(r.ROOT))]=r.digest(path.read_bytes())
        for f in r.read(path)['families']:
            label=f.get('family',f.get('fibration_id'))
            if label in spec['families']:
                families.append({'family':label,'A':f['A_coefficients_low_to_high'],'B':f['B_coefficients_low_to_high']})
    assert len(families)==3
    labels={x['token']:x for x in r.read(base.MANIFEST)['rows']};models={x['token']:x for x in r.read(base.INPUT)['cases']}
    remainders={x['token']:x['factor']['unresolved_cofactor'] for x in r.read(retained.OUTPUT)['rows']}
    cases=[{'token':t,'family':labels[t]['family'],'parameter':labels[t]['parameter'],
            'model':models[t]['model'],'unresolved_cofactor':remainders[t]} for t in spec['cases']]
    files=(Path(__file__),PROTOCOL,base.INPUT,base.MANIFEST,retained.OUTPUT)
    sources.update({str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files})
    r.write_new(INPUT,{'schema':'rank-jump.fresh-symbolic-discriminant-inputs.v1','families':families,'cases':cases,'source_hashes':sources})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,Path(r.__file__))}


def worker(family):
    from sage.all import QQ,ZZ,GF,PolynomialRing,prime_range
    inp=r.read(INPUT);spec=r.read(PROTOCOL)['limits'];f=next(x for x in inp['families'] if x['family']==family)
    R=PolynomialRing(QQ,'t');A=R(list(map(QQ,f['A'])));B=R(list(map(QQ,f['B'])));D=-16*(4*A**3+27*B**2)
    fac=D.factor();product=R(fac.unit());factors=[]
    for q,e in fac:
        product*=q**e;possible=set(range(1,int(q.degree())));modular=[]
        primitive=q*q.denominator()
        content=ZZ(0)
        for c in primitive:content=content.gcd(ZZ(c))
        primitive/=content
        for p in prime_range(spec['modular_prime_bound']+1):
            if not possible:break
            fp=primitive.change_ring(GF(p))
            if fp.degree()!=q.degree() or fp.gcd(fp.derivative())!=1:continue
            factors_p=fp.factor();degrees=[int(g.degree()) for g,k in factors_p for _ in range(k)]
            sums={0}
            for d in degrees:sums|={s+d for s in list(sums)}
            possible&=sums-{0,int(q.degree())}
            modular.append({'prime':int(p),'factors_ascending':[[int(c) for c in g.monic()] for g,k in factors_p],
                            'factor_degrees':degrees,'remaining_possible_proper_factor_degrees':sorted(possible)})
        factors.append({'coefficients_ascending':list(map(str,q.list())),'exponent':int(e),
                        'primitive_integer_coefficients_ascending':list(map(str,primitive.list())),
                        'degree':int(q.degree()),'modular_irreducibility_status':'PASS' if not possible else 'UNKNOWN',
                        'remaining_possible_proper_factor_degrees':sorted(possible),'modular':modular})
    assert product==D
    rows=[]
    for c in inp['cases']:
        if c['family']!=family:continue
        t=QQ(c['parameter']);den=t.denominator();rawA=A(t)*den**8;rawB=B(t)*den**12
        model,_=r.short(c['model'],[]);targetA=QQ(str(model[3]));targetB=QQ(str(model[4]))
        u2=rawB*targetA/(targetB*rawA);assert u2>0 and u2.is_square();u=u2.sqrt()
        assert rawA==u**4*targetA and rawB==u**6*targetB
        assert D(t)*den**24==u**12*(-16*(4*targetA**3+27*targetB**2))
        N=ZZ(c['unresolved_cofactor']);gcds=[];parts=[N]
        for factor in factors:
            q=R(list(map(QQ,factor['coefficients_ascending'])));value=q(t)*den**q.degree()
            value=ZZ(value.numerator());g=N.gcd(value)
            gcds.append({'factor_degree':factor['degree'],'homogeneous_numerator':str(value),'gcd':str(g)})
            refined=[]
            for n in parts:
                g=n.gcd(value)
                refined.extend([g,n//g] if 1<g<n else [n])
            parts=refined
        assert __import__('math').prod(parts)==N
        rows.append({'token':c['token'],'parameter':str(t),'scale_to_frozen_short_model':str(u),
                     'factor_value_gcds':gcds,'cofactor_parts':list(map(str,parts)),
                     'proper_split_count':len(parts)-1})
    return {'status':'PASS','family':family,'discriminant_coefficients_ascending':list(map(str,D.list())),
            'discriminant_degree':int(D.degree()),'unit':str(fac.unit()),'factors':factors,'specializations':rows,
            'boundary':'Polynomial decomposition and integer gcd only. No class incidence or rational-solubility assertion.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for family in r.read(PROTOCOL)['families']:
        path=WORK/f'{family}.json'
        if not path.exists():
            with (WORK/f'{family}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--family',family],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['worker_seconds_per_family'])
                    error=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:error='bounded timeout'
            if error:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','family':family,'reason':error})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(family,row['status'],[(x['degree'],x['modular_irreducibility_status']) for x in row.get('factors',[])],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-symbolic-discriminant.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        r.write_new(WORK/f'{args.family}.json',{'bindings':bindings(),**worker(args.family)})
