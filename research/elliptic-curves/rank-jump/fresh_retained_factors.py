#!/usr/bin/env python3
"""Recover equation-only factor hints for incomplete frozen-panel rows."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base
import fresh_governing_completion as completion
import fresh_strict_boundary_coordinates as coords

PROTOCOL=Path(__file__).with_name('FRESH_RETAINED_FACTOR_PROTOCOL.json')
HINTS=r.OUT/'rank_jump_fresh_retained_factor_hints_v1.json'
OUTPUT=r.OUT/'rank_jump_fresh_retained_factor_supplement_v1.json'
REPORT=r.OUT/'rank_jump_fresh_retained_factor_comparison_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-retained-factor-v1'
PRIMES=r.OUT/'record_prime_factor_proofs_20260904.json'
HISTORIC=r.OUT/'exceptional_selmer_feasibility_v1.json'
INVENTORY=r.OUT/'inventory187_conductor_bounds_v2.json'


def export():
    labels={x['token']:x for x in r.read(base.MANIFEST)['rows']}
    inventory={x['id']:x for x in r.read(INVENTORY)['rows']}
    partial={x['token']:x for x in r.read(completion.OUTPUT)['rows']};rows=[]
    for token in r.read(PROTOCOL)['cases']:
        label=labels[token];hints=[2,3]
        if token=='case-13':
            values=[p for p,e in r.read(PRIMES)['records']['356']['factorizations']['DISCRIMINANT_FACTORIZATION']]
            path=PRIMES;pointer='/records/356/factorizations/DISCRIMINANT_FACTORIZATION'
        elif token=='case-14':
            data=r.read(HISTORIC)['historical_probes'];index=next(i for i,x in enumerate(data) if x['curve_id']==385)
            values=data[index]['record']['input']['factor_hint_primes']
            path=HISTORIC;pointer=f'/historical_probes/{index}/record/input/factor_hint_primes'
        elif label['id'] in inventory:
            row=inventory[label['id']];values=[x['prime'] for x in row['local_data']]+[row['remaining_cofactor']]
            path=INVENTORY;pointer=f'rows[id={label["id"]}]:local_data.prime,remaining_cofactor'
        else:
            values=partial[token]['tested_bad_primes'];path=completion.OUTPUT;pointer=f'rows[token={token}].tested_bad_primes'
        hints=sorted(set(hints+list(map(int,values))))
        rows.append({'token':token,'hints':hints,'source':str(path.relative_to(r.ROOT)),
                     'source_sha256':r.digest(path.read_bytes()),'whitelist_pointer':pointer})
    paths=(Path(__file__),PROTOCOL,base.MANIFEST)
    r.write_new(HINTS,{'schema':'rank-jump.fresh-retained-factor-hints.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'scope':'Only equation factor hints projected. Hints are not assumed prime or complete. Rank labels and exceptional points are excluded from the worker payload.'})


def bindings():
    paths=(Path(__file__),PROTOCOL,HINTS,base.INPUT,Path(base.__file__),base.PROTOCOL,base.LOCAL,
           Path(coords.__file__),Path(coords.prior.__file__),Path(r.__file__))
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def factor(token):
    from sage.all import ZZ
    f,pts,scale=base.model_data(token);D=ZZ(16*f.discriminant());remaining=abs(D);factors=[];trials=[]
    hints=next(x['hints'] for x in r.read(HINTS)['rows'] if x['token']==token)
    for p in hints:
        if p<=1 or remaining%p:continue
        prime=bool(ZZ(p).is_prime(proof=True));trials.append({'hint':p,'proved_prime':prime})
        if not prime:continue
        e=0
        while remaining%p==0:remaining//=p;e+=1
        factors.append([p,e])
    if remaining>1:
        prime=bool(remaining.is_prime(proof=True));trials.append({'remaining':str(remaining),'proved_prime':prime})
        if prime:factors.append([int(remaining),1]);remaining=ZZ(1)
    factors.sort();assert remaining*__import__('math').prod(p**e for p,e in factors)==abs(D)
    return {'status':'PASS' if remaining==1 else 'UNKNOWN','integral_cubic_ascending':list(map(str,f.list())),
        'scale':str(scale),'elliptic_discriminant':str(D),'factors':factors,'unresolved_cofactor':str(remaining),
        'primality_checks':trials,'reason':None if remaining==1 else 'Retained hints leave a composite cofactor; no general factorization attempted'}


def setup(token):
    from sage.all import ZZ,pari
    fact=r.read(WORK/f'{token}-factor.json');local=r.read(WORK/f'{token}-local.json')
    assert fact['status']=='PASS' and local['status'] in ('PASS','PARTIAL')
    f,pts,scale=base.model_data(token);product=ZZ(1)
    for p,e in fact['factors']:assert ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert product==abs(16*f.discriminant());primes=[p for p,e in fact['factors']]
    nf=pari.nfinit([pari(f),primes])
    return {'factor':fact,'local':local},f,pts,primes,nf


def boundary(token):
    coords.prior.setup=setup
    return coords.worker(token)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(PROTOCOL)['cases']:
        row={'token':token}
        for stage in ('factor','local','boundary'):
            path=WORK/f'{token}-{stage}.json'
            if not path.exists():
                if stage!='factor' and row['factor']['status']!='PASS':
                    r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':'Factor coverage incomplete'})
                elif stage=='boundary' and row['local']['status'] not in ('PASS','PARTIAL'):
                    r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':'Generic local coverage incomplete'})
                else:
                    with (WORK/f'{token}-{stage}.log').open('x') as log:
                        try:
                            proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token,'--stage',stage],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits'][stage+'_seconds'])
                            error=None if proc.returncode==0 else 'worker failure'
                        except subprocess.TimeoutExpired:error='bounded worker timeout'
                    if error:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':error})
            value=r.read(path);assert value['bindings']==bindings();row[stage]=value
            print(token,stage,value['status'],value.get('additional_boundary_capacity_upper_bound',value.get('strict_generic_dimension')),flush=True)
        rows.append(row)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-retained-factor-supplement.v1','bindings':bindings(),'rows':rows})


def report():
    labels={x['token']:x for x in r.read(base.MANIFEST)['rows']};rows=[]
    for x in r.read(OUTPUT)['rows']:
        y=labels[x['token']];b=x['boundary'];local=x['local']
        row={k:y[k] for k in ('token','id','generic_rank','retained_rank_lower_bound','observed_quotient_rank')}
        if b['status']!='PASS':row.update(status='UNKNOWN',reason=b['reason'],unresolved_cofactor=x['factor'].get('unresolved_cofactor'))
        else:
            k=b['generic_strict_dimension'];h=b['Selmer_boundary_dimension_interval'][1];a=b['additional_boundary_capacity_upper_bound'];R=y['retained_rank_lower_bound']
            row.update(status='PASS',generic_strict_dimension=k,local_product_dimension=b['local_point_product_dimension'],
                boundary_upper_bound=h,additional_boundary_capacity_upper_bound=a,
                necessary_additional_strict_rational_dimension=max(0,R-h-k),
                rank_label_dependent_S_class_dimension_lower_bound=max(k,R-h),
                inherited_CT_switch_rank=local['minus_twist_CT_rank'],
                additional_Selmer_dimension_formula=f'(c_S - {k}) + e, 0 <= e <= {a}',
                rank_upper_bound_formula=f'c_S + {h}',additional_quotient_CT='UNKNOWN')
        rows.append(row)
    paths=(Path(__file__),OUTPUT,base.MANIFEST)
    r.write_new(REPORT,{'schema':'rank-jump.fresh-retained-factor-comparison.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'boundary':'Rank-derived necessities only. No exceptional class representative or full additional CT supplied.'})
    for row in rows:print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','report']);p.add_argument('--token');p.add_argument('--stage');args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='capture':capture()
    elif args.mode=='report':report()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        base.WORK=WORK
        result={'factor':factor,'local':base.local_worker,'boundary':boundary}[args.stage](args.token)
        r.write_new(WORK/f'{args.token}-{args.stage}.json',{'bindings':bindings(),**result})
