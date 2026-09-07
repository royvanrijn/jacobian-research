#!/usr/bin/env sage-python
"""Fixed special-prime lattice pilot on 64 untouched factor-base directions."""
import argparse
import json
from math import gcd, prod
from pathlib import Path
import runpy
import time
from sage.all import QQ, pari
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,run
from audit_bnf_free_s_class_quotient import packed_rank

ROOT,ART,cert = batch.ROOT,batch.ART,batch.cert
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-special-primes-v1'
OUT = ART/'small_conductor_special_primes_v1.json'
AUDIT = ART/'small_conductor_special_prime_relations_v1.json'
BASE = ART/'small_conductor_norm_batch_relations_v1.json'
SOURCE = Path(__file__).resolve()


def normalize(row,blocks):
    result=set(row)
    for pivot,support in blocks.items():
        if pivot in result: result.symmetric_difference_update(support)
    return result


def setup():
    data=cert.read(batch.forms.OUT)
    poly=pari.Pol(list(map(int,data['original_monic_cubic_descending'])))
    primes=[int(p) for p,e in cert.read(batch.forms.target.PROOF)['discriminant_factorization']]
    pari.addprimes(primes)
    nf=pari.nfinit([poly,primes])
    if len(pari.nfcertify(nf)) or str(nf.disc())!=data['field_discriminant']:
        raise ArithmeticError('maximal order differs')
    x=pari.Mod('x',poly)
    w=sum(pari(QQ(c))*x**i for i,c in enumerate(data['integral_norm_generator']['w_power_basis']))
    return data,poly,nf,x,w


def reduce_lattice(q,r,hessian):
    A,B,C=map(int,hessian)
    def dot(u,v): return 2*A*u[0]*v[0]+B*(u[0]*v[1]+u[1]*v[0])+2*C*u[1]*v[1]
    u,v=[q,0],[r,1]
    for _ in range(1000):
        if dot(v,v)<dot(u,u): u,v=v,u
        length=dot(u,u); inner=dot(u,v)
        if 2*abs(inner)<=length: break
        k=(2*inner+length)//(2*length)
        v=[v[i]-k*u[i] for i in range(2)]
    else: raise ArithmeticError('lattice reduction cap')
    if abs(u[0]*v[1]-u[1]*v[0])!=q or any((t[0]-r*t[1])%q for t in [u,v]):
        raise ArithmeticError('index-q lattice identity failed')
    return [u,v]


def selection():
    base=cert.read(BASE)
    blocks={max(r['parity_columns']):set(r['parity_columns']) for r in base['canonical_rational_relations']}
    covered=set()
    for r in base['relations']:
        covered.update(normalize([c for c,e in r['ideal_factorization'] if e%2],blocks))
    data,poly,nf,x,w=setup()
    a=int(data['integral_norm_generator']['fixed_a']); M=data['sl2_matrix']
    eligible=[]; seen=set()
    for i,col in enumerate(base['columns']):
        q=col['p']
        if i in covered or i in blocks or col['f']!=1 or col['e']!=1 or q<1000 or q in seen: continue
        if int(nf.disc())%q==0 or int(data['defining_order_index'])%q==0 or a%q==0: continue
        eligible.append(i);seen.add(q)
    targets=[]
    for k in range(64):
        column=eligible[k*(len(eligible)-1)//63]
        col=base['columns'][column];q=col['p']
        roots=[int(pari.lift(r)) for r in pari.polrootsmod(poly,q)]
        root=next(r for r in roots if str(pari.idealhnf(nf,q,x-r))==col['hnf'])
        wbar=0
        for i,c in enumerate(data['integral_norm_generator']['w_power_basis']):
            c=QQ(c);wbar+=(int(c.numerator())*pow(int(c.denominator()),-1,q)*pow(root,i,q))%q
        bm=(a*M[0][0]+wbar*M[1][0])%q
        bn=(a*M[0][1]+wbar*M[1][1])%q
        r=(-bn*pow(bm,-1,q))%q
        if sum(int(c)*r**(3-i) for i,c in enumerate(data['reduced_binary_cubic_descending']))%q:
            raise ArithmeticError('target congruence does not divide norm')
        targets.append({'column':column,'p':q,'hnf':col['hnf'],'theta_root':root,
                        'norm_slope':r,'lattice_basis':reduce_lattice(q,r,data['reduced_hessian'])})
    return {'eligible_distinct_rational_primes':len(eligible),'previously_covered_quotient_columns':len(covered),'targets':targets}


def prepare():
    if (D/'protocol.json').exists(): raise FileExistsError('preserve protocol')
    data=selection()
    sources=[SOURCE,BASE,batch.forms.OUT,Path(batch.__file__).resolve(),ROOT/'elliptic-curves/cas/audit_bnf_free_s_class_quotient.py']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.small-conductor-special-primes-protocol.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sources},
        'selection':data,'coefficient_box':64,'smooth_bound':400000,'workers':1,'rss_bytes':1610612736,
        'stage_seconds':{'worker':180,'check':300,'audit':180},
        'gate':'The preceding296 norm rows are independent but leave quotient dimension27847. Target directions untouched after exact canonical rational-prime elimination.',
        'selection_rule':'Exclude canonical pivots and covered columns; choose one unramified degree-one ideal per rational prime>=1000, excluding index/discriminant/fixed-a divisors. Select64 equally spaced indices in ascending existing-column order, including endpoints.',
        'population':'For every target use the exact Hessian-reduced index-q basis; enumerate coprime(u,v), -64<=u<=64,1<=v<=64; retain only primitive mapped(m,n). Preserve occurrences across targets, including duplicates.',
        'claim_boundary':'A finite ideal-relation pilot on the same curve. Class-rank generation remains GRH-conditional; no exact rank, unconditional upper bound, or new point follows.'})
    print('PREPARED64 TARGETS',data['targets'][0]['p'],data['targets'][-1]['p'],flush=True)


def protocol():
    p=cert.read(D/'protocol.json')
    for name,h in p['sources'].items():
        if cert.hashed(ROOT/name)!=h: raise ArithmeticError('frozen source changed')
    return p


def calculate(check=False):
    p=protocol(); primes=batch.prior.primes_to(p['smooth_bound']); primorial=prod(primes)
    c=list(map(int,cert.read(batch.forms.OUT)['reduced_binary_cubic_descending']))
    chunks=[]; smooth=[]; count=0; started=time.monotonic()
    for index,target in enumerate(p['selection']['targets']):
        v1,v2=target['lattice_basis']; b=p['coefficient_box'];pairs=[]
        for u in range(-b,b+1):
            for v in range(1,b+1):
                if gcd(u,v)!=1: continue
                m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                if gcd(m,n)!=1: continue
                if (m-target['norm_slope']*n)%target['p']: raise ArithmeticError('lattice congruence failed')
                pairs.append((u,v,m,n))
        values=[sum(a*m**(3-i)*n**i for i,a in enumerate(c)) for u,v,m,n in pairs]
        if any(v==0 or v%target['p'] for v in values): raise ArithmeticError('forced norm factor missing')
        if check: remainders=[batch.strip(abs(v),primorial) for v in values]
        else:
            supports=[gcd(abs(v),r) for v,r in zip(values,batch.residues(list(map(abs,values)),primorial))]
            remainders=[batch.strip(abs(v),s) for v,s in zip(values,supports)]
        records=[]
        for (u,v,m,n),value,remainder in zip(pairs,values,remainders):
            record={'u':u,'v':v,'m':m,'n':n,'value':str(value),'remainder':str(remainder)}
            if remainder==1:
                record['factorization']=batch.factor_smooth(abs(value),primes)
                smooth.append({'target_index':index,**record})
            records.append(record)
        path=D/('target_%02d.json'%index);result={'target':target,'records':records}
        if check:
            if cert.read(path)!=result: raise ArithmeticError('scalar replay differs')
        else:
            if path.exists(): raise FileExistsError('preserve target chunk')
            checkpoint(path,result)
        count+=len(records)
        chunks.append({'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path),'pairs':len(records),'smooth':sum(r['remainder']=='1' for r in records)})
        if not check: checkpoint(D/'progress.json',{'completed_targets':index+1,'pairs':count,'smooth':len(smooth),'wall_seconds':time.monotonic()-started})
        if (index+1)%8==0: print('TARGETS',index+1,'PAIRS',count,'SMOOTH',len(smooth),flush=True)
    elapsed=time.monotonic()-started
    result={'schema':'elliptic-curves.small-conductor-special-primes.v1','status':'PASS','protocol':p,
            'chunks':chunks,'candidate_occurrences':count,'smooth_occurrences':len(smooth),
            'smooth_records':smooth,'wall_seconds':elapsed,'claim_boundary':p['claim_boundary']}
    if check:
        result['wall_seconds']=cert.read(OUT)['wall_seconds']
        if cert.read(OUT)!=result: raise ArithmeticError('smooth summary differs')
        checkpoint(D/'scalar_replay.json',{'status':'PASS','source_sha256':cert.hashed(OUT),'wall_seconds':elapsed})
    else:
        if OUT.exists(): raise FileExistsError('preserve norm results')
        checkpoint(OUT,result)
    print('SPECIAL PRIME', 'CHECK' if check else 'WORKER','PASS',count,len(smooth),elapsed,flush=True)


def audit(check=False):
    p=protocol(); summary=cert.read(OUT); scalar=cert.read(D/'scalar_replay.json')
    if scalar['status']!='PASS' or scalar['source_sha256']!=cert.hashed(OUT): raise ArithmeticError('scalar replay required')
    for row in summary['chunks']:
        if cert.hashed(ROOT/row['path'])!=row['sha256']: raise ArithmeticError('chunk differs')
    if selection()!=p['selection']: raise ArithmeticError('target-selection replay differs')
    prior=runpy.run_path(str(ROOT/'elliptic-curves/cas/audit_small_conductor_norm_batch.sage'))
    base=prior['expected'](check=True)
    if base!=cert.read(BASE): raise ArithmeticError('previous matrix replay failed')
    data,poly,nf,x,w=setup();a=int(data['integral_norm_generator']['fixed_a']);M=data['sl2_matrix']
    lookup={(r['p'],r['hnf']):i for i,r in enumerate(base['columns'])}
    canonical=[sum(1<<i for i in r['parity_columns']) for r in base['canonical_rational_relations']]
    old=[sum((e%2)<<c for c,e in r['ideal_factorization']) for r in base['relations']]
    blocks={max(r['parity_columns']):set(r['parity_columns']) for r in base['canonical_rational_relations']}
    old_rank=packed_rank(canonical+old);new=[];relations=[];hit=set();projections=[]
    selected=[t['column'] for t in p['selection']['targets']]
    for record in summary['smooth_records']:
        target=p['selection']['targets'][record['target_index']];m,n=record['m'],record['n'];value=int(record['value'])
        if prod(q**e for q,e in record['factorization'])!=abs(value) or any(q>p['smooth_bound'] or not pari.isprime(q) for q,e in record['factorization']): raise ArithmeticError('norm factorization differs')
        beta=a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
        if pari.nfeltnorm(nf,beta)!=a*a*value: raise ArithmeticError('principal norm identity differs')
        factor=pari.idealfactor(nf,beta)
        if pari.idealhnf(nf,pari.idealfactorback(nf,factor))!=pari.idealhnf(nf,beta): raise ArithmeticError('principal ideal product differs')
        row=[]
        for j in range(factor.nrows()):
            ideal,e=factor[j,0],int(factor[j,1]);column=lookup[(int(ideal[0]),str(pari.idealhnf(nf,ideal)))]
            if e<0: raise ArithmeticError('nonintegral element')
            row.append([column,e])
        if not any(c==target['column'] and e>0 for c,e in row): raise ArithmeticError('intended prime ideal absent')
        sparse=normalize([c for c,e in row if e%2],blocks)
        hit.update(set(selected)&sparse)
        projections.append(sum(1<<i for i,c in enumerate(selected) if c in sparse))
        new.append(sum((e%2)<<c for c,e in row))
        relations.append({'target_index':record['target_index'],'m':m,'n':n,
            'beta_power_basis':[str(pari.polcoef(pari.lift(beta),i)) for i in range(3)],'ideal_factorization':row})
    gain=packed_rank(canonical+old+new)-old_rank;dimension=len(base['columns'])-old_rank-gain
    unique={(m,n) if (n>0 or n==0 and m>0) else (-m,-n) for m,n in [(r['m'],r['n']) for r in relations]}
    sources=[SOURCE,BASE,OUT,D/'scalar_replay.json',D/'protocol.json']
    result={'schema':'elliptic-curves.small-conductor-special-prime-relations.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sources},'relations':relations,
            'smooth_occurrences':len(relations),'unique_smooth_pairs_up_to_sign':len(unique),
            'additional_independent_rows':gain,'target_columns_reached':len(hit),'target_projection_rank':packed_rank(projections),
            'previous_quotient_dimension':len(base['columns'])-old_rank,'new_quotient_dimension':dimension,
            'conditional_on_grh_class_two_rank_upper_bound':dimension,
            'conditional_on_grh_curve_rank_upper_bound':(dimension+7)//2*2,
            'additional_independent_rows_needed_for_target':max(0,dimension-16),
            'unconditional_rank_upper_bound':None,'exact_rank':None,'claim_boundary':p['claim_boundary']}
    if check:
        if cert.read(AUDIT)!=result: raise ArithmeticError('ideal relation replay differs')
    else:
        if AUDIT.exists(): raise FileExistsError('preserve audit')
        checkpoint(AUDIT,result)
    print('TARGETED MATRIX PASS',gain,'NEW ROWS',len(hit),'TARGETS REACHED','DIMENSION',dimension,flush=True)


def launch(stage):
    p=protocol(); path=D/(stage+'.supervisor.json')
    if path.exists(): raise FileExistsError('preserve supervisor record')
    result=run(['/home/royvanrijn/.local/bin/sage','-python',str(SOURCE),stage],
        limits=Limits(p['stage_seconds'][stage],p['rss_bytes']),cwd=ROOT,
        log_path=D/(stage+'.log'),checkpoint_path=path)
    print(stage,result['outcome'],result['returncode'],flush=True)
    if result['outcome']!='completed' or result['returncode']!=0: raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','worker','check','audit','audit-check','launch-worker','launch-check','launch-audit'])
    args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage.startswith('launch-'):launch(args.stage[7:])
    elif args.stage.startswith('audit'):audit(args.stage=='audit-check')
    else:calculate(args.stage=='check')
