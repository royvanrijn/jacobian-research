#!/usr/bin/env sage-python
"""GRH generation by an interval-certified triangular explicit-formula test."""
import argparse
import json
from pathlib import Path
import runpy
from sage.all import RealIntervalField
import extend_small_conductor_norm_batch as batch
from audit_bnf_free_s_class_quotient import packed_rank
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,run

ROOT,ART,cert=batch.ROOT,batch.ART,batch.cert
D=ROOT/'artifacts/local/elliptic-curves/small-conductor-smaller-base-v1'
OUT=ART/'small_conductor_smaller_base_v1.json'
BASE=ART/'small_conductor_norm_batch_relations_v1.json'
EXTRA=ART/'small_conductor_special_prime_relations_v1.json'
SOURCE=Path(__file__).resolve()


def interval_record(x):
    return {'lower':str(x.lower()),'upper':str(x.upper())}


def triangle(base,T,bits,terms):
    R=RealIntervalField(bits);L=R(T).log()
    # For F(x)=max(1-x/L,0), expand the archimedean kernels in exp(-(k+1/2)x).
    # I = (pi^2/2 - sum exp(-(k+1/2)L)/(k+1/2)^2)/L.
    # J = pi/2 -4*Catalan/L + sum (-1)^k exp(-(k+1/2)L)/(L*(k+1/2)^2).
    e=L.__neg__().exp();a=(-L/2).exp();positive=R(0);alternating=R(0)
    for k in range(terms):
        term=a/(R(k)+R(1)/2)**2
        positive+=term;alternating+=(-1)**k*term;a*=e
    tail=a/(R(terms)+R(1)/2)**2/(1-e)
    error=R(-tail.upper(),tail.upper())
    I=(R.pi()**2/2-positive+error)/L
    J=R.pi()/2-4*R.catalan_constant()/L+(alternating+error)/L
    total=R(0);prime_ideals=0;prime_powers=0;logs={}
    for col in base['columns']:
        p,f=col['p'],col['f'];q=p**f
        if q>=T:continue
        prime_ideals+=1
        if p not in logs:logs[p]=R(p).log()
        logq=f*logs[p];power=q;m=1
        while power<T:
            total+=logq*(1-m*logq/L)/R(power).sqrt()
            prime_powers+=1;power*=q;m+=1
    logdisc=R(int(base['field_discriminant'])).log()
    margin=2*total-logdisc+3*(R.euler_constant()+(8*R.pi()).log())-3*I+3*J
    return {'cutoff':T,'prime_ideals_contributing':prime_ideals,'prime_power_terms':prime_powers,
            'I_interval':interval_record(I),'J_interval':interval_record(J),
            'prime_sum_interval':interval_record(total),'margin_interval':interval_record(margin),
            'positive_margin_certified':bool(margin.lower()>0)}


def normalized_rows(base,extra):
    canonical={max(r['parity_columns']):set(r['parity_columns']) for r in base['canonical_rational_relations']}
    rows=[]
    for r in base['relations']+extra['relations']:
        s={c for c,e in r['ideal_factorization'] if e%2}
        for c in list(s):
            if c in canonical:s.symmetric_difference_update(canonical[c])
        rows.append(s)
    return canonical,rows


def intersection(base,canonical,rows,T):
    pivots={};inside_rows=[];combinations=[]
    for i,s in enumerate(rows):
        outside=sum(1<<c for c in s if base['columns'][c]['p']>T)
        inside=sum(1<<c for c in s if base['columns'][c]['p']<=T)
        combination=1<<i
        while outside:
            j=outside.bit_length()-1
            if j not in pivots:
                pivots[j]=(outside,inside,combination);break
            o,v,k=pivots[j];outside^=o;inside^=v;combination^=k
        if not outside:
            actual=set()
            for j,row in enumerate(rows):
                if (combination>>j)&1:actual.symmetric_difference_update(row)
            if sum(1<<c for c in actual)!=inside or any(base['columns'][c]['p']>T for c in actual):
                raise ArithmeticError('outside cancellation certificate failed')
            inside_rows.append(inside);combinations.append({'input_row_indices':[j for j in range(len(rows)) if (combination>>j)&1],
                                                           'result_columns':sorted(actual)})
    free=sum(c['p']<=T and i not in canonical for i,c in enumerate(base['columns']))
    rank=packed_rank(inside_rows)
    if rank!=packed_rank([sum(1<<c for c in row) for row in rows])-len(pivots):
        raise ArithmeticError('intersection rank identity failed')
    return {'cutoff':T,'free_columns_before_norm_relations':free,'outside_projection_rank':len(pivots),
            'relations_supported_after_exact_elimination':rank,'supported_combinations':combinations,
            'quotient_dimension':free-rank,'rows_needed_for_target16':max(0,free-rank-16)}


def expected():
    protocol=cert.read(D/'protocol.json')
    frozen=cert.read(D/'checker_sources.json')
    for name,h in {**protocol['sources'],**frozen['sources']}.items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('frozen input differs')
    # This replays the full original factor base and all832 principal relations.
    prior=runpy.run_path(str(ROOT/'elliptic-curves/cas/target_small_conductor_prime_ideals.sage'))
    prior['audit'](check=True)
    base,extra=cert.read(BASE),cert.read(EXTRA)
    R=RealIntervalField(protocol['precision_bits']);logdisc=R(int(base['field_discriminant'])).log()
    if not logdisc.upper()<353:raise ArithmeticError('2025 cubic theorem range failed')
    uniform=int(((R(23)/6)*logdisc**2).upper().ceil())
    canonical,rows=normalized_rows(base,extra)
    if packed_rank([sum(1<<c for c in row) for row in rows])!=832:raise ArithmeticError('old row rank differs')
    records=[]
    for T in protocol['cutoffs']:
        test=triangle(base,T,protocol['precision_bits'],protocol['series_terms'])
        algebra=intersection(base,canonical,rows,T)
        records.append({**test,'generation_certified_under_grh':test['positive_margin_certified'] or T>=uniform,
                        'relation_matrix':algebra})
        print('CUTOFF',T,'MARGIN',test['margin_interval'],'DIMENSION',algebra['quotient_dimension'],flush=True)
    valid=[r for r in records if r['generation_certified_under_grh']]
    best=min(valid,key=lambda r:r['relation_matrix']['quotient_dimension'])
    return {'schema':'elliptic-curves.small-conductor-smaller-base.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,BASE,EXTRA,D/'protocol.json',D/'checker_sources.json']},
            'uniform_2025_grh_bound':uniform,'log_discriminant_interval':interval_record(logdisc),
            'explicit_formula_theorem':'Grenie-Molteni2017 Theorem2.4, restating Belabas-Diaz y Diaz-Friedman: triangular F(x)=max(1-x/log(T),0), F(0)=1, nonnegative Fourier transform. Positive displayed margin proves generation under GRH for nontrivial ideal-class characters.',
            'references':['https://arxiv.org/pdf/1607.02430','https://doi.org/10.1090/mcom/4114'],
            'archimedean_tail_bound':'For lambda=k+1/2, sum_{k>=N} exp(-lambda*L)/lambda^2 <= exp(-(N+1/2)*L)/((N+1/2)^2*(1-exp(-L))); symmetric interval error encloses both positive and alternating tails.',
            'records':records,'best_certified_cutoff':best['cutoff'],
            'conditional_on_grh_class_two_rank_upper_bound':best['relation_matrix']['quotient_dimension'],
            'conditional_on_grh_curve_rank_upper_bound':(best['relation_matrix']['quotient_dimension']+7)//2*2,
            'unconditional_rank_upper_bound':None,'exact_rank':None,
            'claim_boundary':'Certified smaller generation set under GRH and exact cancellation of outside columns in existing principal relations. No new relation search, unconditional generation or exact rank. Negative triangle margin means only that this sufficient test did not certify that cutoff.'}


def launch(check=False):
    stage='check' if check else 'build';p=cert.read(D/'protocol.json')
    if (D/(stage+'.supervisor.json')).exists():raise FileExistsError('preserve run')
    command=['/home/royvanrijn/.local/bin/sage','-python',str(SOURCE),'--check' if check else '--build']
    r=run(command,limits=Limits(p['wall_seconds_per_stage'],p['rss_bytes']),cwd=ROOT,
          log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'))
    print(stage,r['outcome'],r['returncode'],flush=True)
    if r['outcome']!='completed' or r['returncode']!=0:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');p.add_argument('--build',action='store_true');p.add_argument('--launch',action='store_true');a=p.parse_args()
    if a.launch:launch(a.check)
    else:
        data=expected()
        if a.check:
            if cert.read(OUT)!=data:raise ArithmeticError('smaller-base replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve certificate')
            checkpoint(OUT,data)
        print('SMALLER BASE PASS',data['best_certified_cutoff'],data['conditional_on_grh_class_two_rank_upper_bound'],flush=True)
