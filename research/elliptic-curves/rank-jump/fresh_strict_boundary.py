#!/usr/bin/env python3
"""Masked local-boundary capacity and separate retrospective strict-block accounting."""
import argparse
from pathlib import Path
from math import prod
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base
import matched103b2_class_boundary as prior

PROTOCOL=Path(__file__).with_name('FRESH_STRICT_BOUNDARY_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_strict_boundary_v1.json'
REPORT=r.OUT/'rank_jump_fresh_strict_boundary_comparison_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-strict-boundary-v1'


def bindings():
    paths=(Path(__file__),PROTOCOL,base.INPUT,base.OUTPUT,Path(base.__file__),Path(prior.__file__),Path(r.__file__))
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def worker(token):
    from sage.all import pari,AA,EllipticCurve,GF,matrix
    old,f,pts,primes,nf=prior.setup(token);loc=old['local']
    assert loc['status'] in ('PASS','PARTIAL') and 'S_finite' in loc
    theta=pari.Mod('z',pari(f));delta=f.discriminant();beta=-pari(delta)*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==delta**4
    E=EllipticCurve([0,0,0,f[1],f[0]])
    S=[2]+[p for p in primes if p!=2 and E.local_data(p,proof=True).conductor_valuation()>0]
    assert S==loc['S_finite']
    omitted=[]
    for p in sorted(set(primes)-set(S)):
        assert p!=2
        vals=[int(pari.idealval(nf,beta,P)) for P in pari.idealprimedec(nf,p)]
        assert all(v%2==0 for v in vals);omitted.append({'prime':p,'valuations':vals})
    generic=[pari(x)-theta for x,y in pts];rows=[];joint=[0]*len(pts);width=0;ell=0;total=[0]*len(pts)
    for p in S:
        place=next(x for x in loc['local'] if x['place']==p)
        decomposition=list(pari.idealprimedec(nf,p));dim=len(decomposition)-1+int(p==2)
        assert dim==place['point_kummer_dimension']
        hilbert=[[int(pari.nfhilbert(nf,beta,g,P)) for P in decomposition] for g in generic]
        assert all(v in (-1,1) for h in hilbert for v in h)
        bits=[int(prod(h)==-1) for h in hilbert]
        for i,bit in enumerate(bits):total[i]^=bit
        sigs=list(map(r.pack,place['signatures']));w=len(place['signatures'][0])
        for i,s in enumerate(sigs):joint[i]|=s<<width
        width+=w;ell+=dim
        rows.append({'place':p,'point_dimension':dim,'hilbert_symbols':hilbert,'functional_on_generic_local_images':bits})
    roots=f.roots(AA,multiplicities=False);signs=[int(-delta*f.derivative()(a)<0) for a in roots]
    real_sigs=[[int(x<a) for a in roots] for x,y in pts]
    real_bits=[sum(a*b for a,b in zip(signs,s))%2 for s in real_sigs]
    for i,bit in enumerate(real_bits):total[i]^=bit;joint[i]|=r.pack(real_sigs[i])<<width
    assert total==[0]*len(pts) # complete reciprocity check on each generic global class
    real_dim=int(len(roots)==3);ell+=real_dim
    assert (delta>0)==bool(real_dim)
    real_basis=[[0,1,1]] if real_dim else []
    real_witness=any(sum(a*b for a,b in zip(signs,s))%2 for s in real_basis)
    rows.append({'place':'infinity','point_dimension':real_dim,'derivative_signs':signs,
                 'generic_signs':real_sigs,'complete_point_basis':real_basis,'functional_on_generic_local_images':real_bits})
    b=int(real_witness or any(any(x['functional_on_generic_local_images']) for x in rows))
    g=int(matrix(GF(2),[[v>>j&1 for j in range(width+len(roots))] for v in joint]).rank())
    m=len(pts);k=m-g;assert k==loc['strict_generic_dimension'];upper=ell-b;assert upper>=g
    return {'status':'PASS','S_finite':S,'polynomial_discriminant':str(delta),
       'derivative_coefficients':[str(pari.lift(beta).polcoef(i)) for i in range(3)],
       'derivative_norm':str(delta**4),'omitted_good_prime_valuations':omitted,'local':rows,
       'generic_dimension':m,'generic_local_dimension':g,'generic_strict_dimension':k,
       'local_point_product_dimension':ell,'reciprocity_constraint_rank_lower_bound':b,
       'Selmer_boundary_dimension_interval':[g,upper],'additional_boundary_capacity_upper_bound':upper-g,
       'generic_reciprocity_checks':len(pts),'additional_strict_dimension':'UNKNOWN',
       'additional_quotient_CT':'UNKNOWN',
       'boundary':'One equation-defined constraint at most; zero means no witness in this test, not absence of other constraints. No additional class or rationality certificate.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for old in r.read(base.OUTPUT)['rows']:
        token=old['token'];path=WORK/f'{token}.json'
        if not path.exists():
            if old['factor']['status']!='PASS':
                r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':'Existing discriminant factorization incomplete; no retry'})
            else:
                with (WORK/f'{token}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['worker_seconds_per_case'])
                        failure='worker failure' if proc.returncode else None
                    except subprocess.TimeoutExpired:failure='bounded worker timeout'
                if failure:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':failure})
        row=r.read(path);assert row['bindings']==bindings();rows.append({'token':token,**row})
        print(token,row['status'],row.get('additional_boundary_capacity_upper_bound'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-strict-boundary.v1','bindings':bindings(),'rows':rows})


def report():
    data=r.read(OUTPUT);labels={x['token']:x for x in r.read(base.MANIFEST)['rows']};rows=[]
    for x in data['rows']:
        y=labels[x['token']];row={k:y[k] for k in ('token','id','family','parameter','generic_rank','retained_rank_lower_bound','observed_quotient_rank')}
        if x['status']!='PASS':row.update(status='UNKNOWN',reason=x['reason'])
        else:
            m=y['generic_rank'];R=y['retained_rank_lower_bound'];k=x['generic_strict_dimension'];g=m-k
            assert m==x['generic_dimension'];h=x['Selmer_boundary_dimension_interval'][1];a=h-g
            total_min=max(k,R-h);extra_min=total_min-k
            assert extra_min==max(0,R-m-a)
            row.update(status='PASS',generic_strict_dimension=k,local_product_dimension=x['local_point_product_dimension'],
                constraint_rank_lower_bound=x['reciprocity_constraint_rank_lower_bound'],boundary_upper_bound=h,
                additional_boundary_capacity_upper_bound=a,necessary_total_rational_strict_dimension=total_min,
                necessary_additional_rational_strict_dimension=extra_min,
                S_class_dimension_lower_bound_from_generic_classes=k,
                S_class_dimension_lower_bound_after_rank_label=total_min,
                rank_upper_bound_formula=f'c_S + {h}',
                additional_Selmer_dimension_formula=f'(c_S - {k}) + e, 0 <= e <= {a}',
                hypothetical_matching_S_class_upper_bound=R-h if R-h>=k else None)
        rows.append(row)
    files=(Path(__file__),OUTPUT,base.MANIFEST)
    r.write_new(REPORT,{'schema':'rank-jump.fresh-strict-boundary-comparison.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Rank-derived necessary dimensions are joined after masked arithmetic. They are not independently measured class dimensions, new ranks, or rational-solubility predictors. Low gains are censored.'})
    for row in rows:print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','report']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='capture':capture()
    elif args.mode=='report':report()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        r.write_new(WORK/f'{args.token}.json',{'bindings':bindings(),**worker(args.token)})
