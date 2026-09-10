#!/usr/bin/env python3
"""Frozen X1092: conic pairs -> positive-rank base -> certified M19 -> existing V3.

Use a NEW folder. `run` prepares if necessary; `resume` uses the sealed plan.
No current campaign is modified. Failed/unfinished stages are UNKNOWN and are
never silently rerun. The search is opt-in, finite, and makes no model calls.
"""
from __future__ import annotations

import argparse
import fcntl
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from broad_rank_runtime import read, sha, write, normalize_parent, guard
from elkies_pair_policy import (require, pair_order, same_quadratic_class, split_controls,
    shortlist, final_selection, mestre_score, rational_height, short_isomorphic)

CAS=Path(__file__).resolve().parent
ROOT=CAS.parents[1]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
SELF=Path(__file__).resolve()
SOURCE_SUFFIXES={'.py','.sage','.gp','.c','.cpp','.h','.sh'}


def environment(sage):
    from run_broad_rank_search import environment as existing_environment
    return existing_environment(sage)


def prepare(args):
    require(sys.platform.startswith('linux'),'Linux process-tree supervision is required')
    sage=shutil.which(args.sage)
    require(sage and Path('/usr/bin/gp').is_file(),'Sage and /usr/bin/gp are required')
    require(2<=args.pool_size<=256 and 1<=args.pair_limit<=4096,'invalid finite pair window')
    require(4<=args.addresses<=4096 and 0<=args.controls<args.addresses and args.ranked>0,'invalid address/selection window')
    require(args.deep_keep>=args.ranked and args.shallow_bound<args.deep_bound<65537,'invalid score bands')
    require(args.shallow_bound>=5 and args.parameter_bits>=16 and 1<=args.v3_calls<=8192,'invalid score/point limits')
    require(args.job_seconds>0 and args.math_seconds>0,'nonpositive process bound')
    folder=args.folder.resolve();require(not folder.exists(),'use a new folder or resume the sealed plan')
    software=json.loads(subprocess.check_output([sage,'-python','-c',
        'import json,sys,numpy,sage.version; from sage.all import pari; '
        'print(json.dumps(dict(python=sys.executable,sage=sage.version.version,pari=str(pari.version()),numpy=numpy.__version__)))'],
        cwd=ROOT,env=environment(sage),text=True,timeout=60).strip().splitlines()[-1])
    folder.mkdir(parents=True);rt=folder/'runtime/research';bindings={}
    def copy(p):
        dest=rt/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,dest);bindings[str(p.relative_to(ROOT))]=sha(p)
    for directory in (CAS,ROOT/'elliptic-curves/ecsearch',ROOT/'elkies-k3/scripts'):
        for p in sorted(directory.rglob('*')):
            if p.is_file() and p.suffix in SOURCE_SUFFIXES and '__pycache__' not in p.parts:copy(p)
    for p in (ROOT/'elliptic-curves').glob('*.py'):copy(p)
    source=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    original=ART/'curve302_recovered_mw17_parent_v1.json'
    parent=normalize_parent(read(source),read(original)['generic_height_gram'],'x1092-elkies-pairs')
    write(rt/'pair-inputs/parent.json',parent)
    bindings.update({str(p.relative_to(ROOT)):sha(p) for p in (source,original)})
    for name in ('curve302_parent_degree2_multisection_orbits_v1.tsv','curve302_parent_degree2_multisection_lattice_v1.json'):
        copy(ART/name)
    # Only equation exclusions are copied. No withheld points, ranks, or targets.
    inv=ROOT/'elliptic-curves/data/research_curves/database.json'
    exclusions=[]
    if inv.exists():
        exclusions=[r['ainvs'] for r in read(inv)['curves']];bindings[str(inv.relative_to(ROOT))]=sha(inv)
    write(rt/'pair-inputs/exclusions.json',{'models':[[str(a) for a in e] for e in exclusions]})
    config={k:getattr(args,k) for k in ('pool_size','pair_limit','addresses','controls','ranked','deep_keep',
        'shallow_bound','deep_bound','parameter_bits','v3_calls','job_seconds','math_seconds')}
    plan={'schema':'x1092-elkies-pairs.v1','root':str(rt),'sage':sage,'software':software,**config,
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'source_inputs':bindings,'height':125000,'rss_bytes':3*1024**3,'workers':1,
        'same_cover_addresses':['0','1','-1','2','-2'],
        'pair_rule':'Bounded same-squareclass prefix, then low-coefficient odd-intersection pairs. Both limits independent.',
        'score_rule':'Good-fibre Mestre sum through shallow bound; extend survivors and score-independent controls to deep bound. No validation primes >=65537 read.',
        'boundary':'No live or historical campaign is altered. Positive base rank is separate from specialized M19 admission. Search bounds do not prove absence.'}
    write(folder/'plan.json',plan)
    write(folder/'manifest.json',{'plan_sha256':sha(folder/'plan.json'),
        'files':{str(p.relative_to(rt)):sha(p) for p in sorted(rt.rglob('*')) if p.is_file()},
        'executables':{str(p):sha(p) for p in (Path(sage),Path('/usr/bin/gp'),Path(software['python']))}})
    print(json.dumps({'status':'PREPARED_NOT_LAUNCHED','folder':str(folder),'software':software}),flush=True)


def context(folder):
    plan,rt=guard(folder)
    return plan,rt,rt/'pair-results',read(rt/'pair-inputs/parent.json')


def bounded_math(seconds,fn):
    def expired(signum,frame):raise TimeoutError('declared per-attempt mathematical limit')
    old=signal.signal(signal.SIGALRM,expired);signal.alarm(seconds)
    try:return fn()
    finally:signal.alarm(0);signal.signal(signal.SIGALRM,old)


def pool_phase(folder):
    import elkies_pair_geometry as geo
    plan,rt,out,parent=context(folder);E,points,G=geo.load_parent(parent)
    table=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
    lattice=read(ART/'curve302_parent_degree2_multisection_lattice_v1.json')
    require(lattice['status']=='PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT' and sha(table)==lattice['orbits_tsv_sha256'],
            'orbit table lost its complete-census certificate')
    rows=[]
    for line in table.read_text().splitlines()[1:]:
        f=line.split('\t')
        if f[1]=='rational':
            require(int(f[2])==10,'wrong rational bisection norm')
            rows.append({'mask':int(f[0]),'word':list(map(int,f[4].split()))})
    require(len(rows)==lattice['rational_bisections']['translation_orbits'],'incomplete orbit roster')
    rows.sort(key=lambda r:(sum(map(abs,r['word'])),max(map(abs,r['word'])),r['mask']))
    selected=rows[:plan['pool_size']];write(out/'pool-selection.json',{'rows':selected,'uses_specialized_points':False})
    results=[]
    for row in selected:
        path=out/'conics'/f"{row['mask']}.json"
        try:
            c=bounded_math(60,lambda:geo.conic(E,points,G,row))
            write(path,c);result={'mask':row['mask'],'status':'PASS_EXACT_CONIC','sha256':sha(path)}
        except (ValueError,ArithmeticError,TimeoutError,ZeroDivisionError) as e:
            result={'mask':row['mask'],'status':'UNKNOWN_BOUNDED_CONSTRUCTION','reason':str(e)}
        results.append(result);write(out/'pool-progress.json',{'results':results},False)
        print('PAIR_CONIC',json.dumps(result),flush=True)
    write(out/'pool.json',{'results':results,'constructed':sum(r['status']=='PASS_EXACT_CONIC' for r in results)})


def conics(out,geo,E,points,G):
    rows=[]
    for r in read(out/'pool.json')['results']:
        if r['status']!='PASS_EXACT_CONIC':continue
        path=out/'conics'/f"{r['mask']}.json";require(sha(path)==r['sha256'],'conic changed')
        c=read(path);geo.verify_conic(c,E,points,G);rows.append(c)
    return rows


def certify_seed(parent,t,w1,w2,c1,c2,destination):
    from sage.all import QQ
    import elkies_pair_geometry as geo
    import parent_foundry_worker as framework
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import short_curve_has_no_rational_2_torsion_modular_certificate
    t,w1,w2=map(QQ,(t,w1,w2));den=t.denominator()
    try:
        model,base=framework.specialize(parent,str(t))
        extras=[geo.lift(c,t,w) for c,w in ((c1,w1),(c2,w2))]
        pts=base+tuple((F(str(x*den**4)),F(str(y*den**6))) for x,y in extras)
        admission=FinitePointAdmission(model,(),prime_bound=1000)
        outcomes=[admission.consider(p) for p in pts]
        torsion=next((p for p in admission.primes if p<=200 and
            short_curve_has_no_rational_2_torsion_modular_certificate(model,p)),None)
        require(len(admission.points)==19 and torsion is not None,'bounded finite M19 witness unresolved')
        proof=checked_rank(model,pts,admission.primes,torsion)
        require(proof['rank_lower_bound']==19,'M19 finite certificate did not close')
        packet={'family':parent['family'],'parameter':str(t),'curve':list(map(str,model)),
                'points':[list(map(str,p)) for p in pts],'proof':proof,'generic_rank':17,'rank_lower_bound':19}
        framework.verify(packet,parent,str(t))
    except (ValueError,ArithmeticError,ZeroDivisionError) as e:
        write(destination/'gate.json',{'status':'UNKNOWN_SPECIALIZED_M19','reason':str(e),'parameter':str(t)})
        return None
    write(destination/'seed.json',packet)
    write(destination/'gate.json',{'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','rank_lower_bound':19,
          'seed_sha256':sha(destination/'seed.json'),'parameter':str(t),
          'boundary':'Rank lower bound only; finite witness failure would not prove dependence.'})
    return packet


def carrier_phase(folder):
    from sage.all import QQ,PolynomialRing,matrix,ZZ,vector
    import elkies_pair_geometry as geo
    plan,rt,out,parent=context(folder);E,points,G=geo.load_parent(parent)
    cs=conics(out,geo,E,points,G);by={c['mask']:c for c in cs}
    order=pair_order(cs,parent['generic_height_gram'],plan['pair_limit']);write(out/'pair-order.json',{'rows':order})
    attempts=[]
    for i,row in enumerate(order):
        a,b=(by[m] for m in row['masks']);dest=out/'pairs'/f'{i:04d}'
        try:
            if row['kind']=='same-cover':
                ring=E.base_ring().ring();field=E.base_ring();T=geo.function(field,a['T']);W=geo.function(field,a['W'])
                scale=QQ(row['scale']);data=None
                for j,u in enumerate(plan['same_cover_addresses']):
                    try:t,w=T(QQ(u)),W(QQ(u))
                    except (ValueError,ZeroDivisionError):continue
                    packet=bounded_math(90,lambda:certify_seed(parent,t,w,scale*w,a,b,dest/f'probe-{j}'))
                    if packet:
                        data={'kind':'rational','masks':row['masks'],'T':a['T'],'W1':a['W'],
                              'scale':str(scale),'function_field_rank_lower_bound':19,
                              'independence':'Exact generic section identities and one specialization with certified independent 19 images.',
                              'witness_seed':str((dest/f'probe-{j}/seed.json').relative_to(rt)),
                              'witness_seed_sha256':sha(dest/f'probe-{j}/seed.json')}
                        break
                require(data is not None,'no independent same-cover specialization in the fixed probe window')
            else:
                data=bounded_math(90,lambda:geo.carrier(a,b))
                w,v=vector(ZZ,a['word']),vector(ZZ,b['word'])
                gram=matrix(QQ,19,19)
                gram[:17,:17]=4*G
                for j,z in enumerate((w,v)):
                    for k,value in enumerate(2*G*z):gram[k,17+j]=gram[17+j,k]=value
                    gram[17+j,17+j]=16
                gram[17,18]=gram[18,17]=w*G*v
                require(gram.is_positive_definite() and gram.det()==4**17*1092*36,'incorrect biquadratic height Gram')
                data['generic_height_gram']=[list(map(str,r)) for r in gram.rows()]
                data['generic_height_determinant']=str(gram.det())
            data.update(conic_sha256={str(c['mask']):sha(out/'conics'/f"{c['mask']}.json") for c in (a,b)},
                        pair_rule=row,uses_known_fibre=False)
            write(dest/'carrier.json',data);write(out/'carrier.json',data)
            attempts.append({'index':i,**row,'status':'PASS_GENERIC_M19_CARRIER'})
            write(out/'carrier-attempts.json',{'results':attempts},False)
            print('PAIR_CARRIER_FOUND',data['kind'],row['masks'],flush=True);return
        except (ValueError,ArithmeticError,ZeroDivisionError,TimeoutError) as e:
            result={'index':i,**row,'status':'UNKNOWN_BOUNDED_PAIR','reason':str(e)}
            write(dest/'unknown.json',result);attempts.append(result)
            write(out/'carrier-attempts.json',{'results':attempts},False)
            print('PAIR_UNKNOWN',json.dumps(result),flush=True)
    write(out/'carrier-status.json',{'status':'NO_CERTIFIED_CARRIER_IN_FINITE_WINDOW','attempts':len(attempts),
                                   'rank_or_existence_upper_bound':None})


def carrier_addresses(data,count,bitcap):
    from sage.all import QQ,EllipticCurve,PolynomialRing
    import elkies_pair_geometry as geo
    if data['kind']=='rational':
        ring=PolynomialRing(QQ,'u');field=ring.fraction_field()
        T,W=geo.function(field,data['T']),geo.function(field,data['W1'])
        for i in range(1,count+1):
            n=(i+1)//2;a=b=1
            for bit in bin(n)[3:]:
                if bit=='0':b+=a
                else:a+=b
            u=QQ((-1 if i%2==0 else 1)*a)/b
            try:yield i,(T(u),W(u),QQ(data['scale'])*W(u)),None
            except (ValueError,ZeroDivisionError):yield i,None,'EXCEPTIONAL_RATIONAL_CHART'
    else:
        J=EllipticCurve(QQ,data['curve']);P=J(list(map(QQ,data['nontorsion']['point'])));current=J(0)
        for i in range(1,count+1):
            if i%2:current+=P
            point=current if i%2 else -current
            if not point.is_zero() and max(rational_height(c) for c in point.xy())>max(16384,bitcap*16):
                yield i,None,'CARRIER_POINT_HEIGHT_CAP';break
            try:yield i,geo.carrier_value(data,point),None
            except (ValueError,ZeroDivisionError):yield i,None,'EXCEPTIONAL_ELLIPTIC_CHART'


def equation(parent,t):
    from compact_atlas_specialization import polynomial
    t=F(str(t));den=t.denominator
    a=polynomial(parent['A_coefficients_low_to_high'],t)*den**8
    b=polynomial(parent['B_coefficients_low_to_high'],t)*den**12
    require(a.denominator==b.denominator==1,'reduced parent must have integral coefficients')
    require(4*a**3+27*b*b!=0,'singular fibre')
    return int(a),int(b)


def traces(model,primes,cache):
    from sage.all import pari
    a,b=model;result=[]
    for p in primes:
        p=int(p);aa,bb=a,b
        while aa%p**4==0 and bb%p**6==0:aa//=p**4;bb//=p**6
        x,y=aa%p,bb%p
        if (4*x**3+27*y*y)%p==0:result.append([p,None]);continue
        key=(p,x,y)
        if key not in cache:
            n=int(pari.ellinit([0,0,0,x,y],p).ellcard());ap=p+1-n
            # Independent character sums check every distinct tiny-prime curve.
            if p<=31:
                direct=-sum(0 if (v:=(z*z*z+x*z+y)%p)==0 else (1 if pow(v,(p-1)//2,p)==1 else -1) for z in range(p))
                require(ap==direct,'PARI/character-sum trace mismatch')
            cache[key]=ap
        result.append([p,cache[key]])
    return result


def score_phase(folder):
    from sage.all import prime_range
    plan,rt,out,parent=context(folder);data=read(out/'carrier.json');rows=[];attrition=[];seen={}
    # Generic invariants only; turn excluded long models into short equations.
    exclusions=[]
    for aa in read(rt/'pair-inputs/exclusions.json')['models']:
        a1,a2,a3,a4,a6=map(F,aa);b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        exclusions.append((-(b2*b2-24*b4)/48,-(-b2**3+36*b2*b4-216*b6)/864))
    for i,value,error in carrier_addresses(data,plan['addresses'],plan['parameter_bits']):
        if error:attrition.append({'index':i,'reason':error});continue
        t,w1,w2=value
        if rational_height(t)>plan['parameter_bits']:
            attrition.append({'index':i,'reason':'PARAMETER_HEIGHT_CAP','bits':rational_height(t)});continue
        try:a,b=equation(parent,t)
        except (ValueError,ZeroDivisionError) as e:attrition.append({'index':i,'reason':str(e)});continue
        j=F(1728*4*a**3,4*a**3+27*b*b)
        if any(short_isomorphic(a,b,c,d) for c,d in seen.get(str(j),[])):
            attrition.append({'index':i,'reason':'EXACT_Q_ISOMORPHISM_DUPLICATE'});continue
        if any(4*c**3+27*d*d!=0 and F(1728*4*c**3,4*c**3+27*d*d)==j and short_isomorphic(a,b,c,d) for c,d in exclusions):
            attrition.append({'index':i,'reason':'PINNED_EQUATION_EXCLUSION'});continue
        seen.setdefault(str(j),[]).append((a,b))
        rows.append({'id':str(i),'parameter':str(t),'w1':str(w1),'w2':str(w2),
                     'model':[str(a),str(b)],'j':str(j),'j_bits':rational_height(j)})
    write(out/'address-roster.json',{'rows':rows,'attrition':attrition,'requested_addresses':plan['addresses'],
          'boundary':'Fixed generator multiples with exact-isomorphism/height attrition, not uniformly random fibres.'})
    if len(rows)<=plan['controls']:
        write(out/'selection.json',{'rows':[],'status':'INSUFFICIENT_FINITE_ADDRESS_WINDOW'});return
    controls,_=split_controls(rows,plan['controls']);write(out/'controls.json',{'ids':[r['id'] for r in controls]})
    cache={};small=list(prime_range(5,plan['shallow_bound']+1));large=list(prime_range(plan['shallow_bound']+1,plan['deep_bound']+1))
    for r in rows:
        r['shallow_traces']=traces(tuple(map(int,r['model'])),small,cache);r['shallow']=mestre_score(r['shallow_traces'])
    write(out/'shallow.json',{'rows':rows,'bound':plan['shallow_bound']})
    deep=shortlist(rows,controls,plan['deep_keep'])
    for r in deep+controls:
        r['extension_traces']=traces(tuple(map(int,r['model'])),large,cache)
        r['deep']=mestre_score(r['shallow_traces']+r['extension_traces'])
    write(out/'deep.json',{'rows':deep+controls,'bound':plan['deep_bound'],'validation_primes_read':[]})
    write(out/'selection.json',{'rows':final_selection(deep,controls,plan['ranked']),'status':'PASS_FIXED_TWO_STAGE_SELECTION'})


def seed_phase(folder,index):
    from sage.all import QQ,PolynomialRing
    import elkies_pair_geometry as geo
    plan,rt,out,parent=context(folder);E,points,G=geo.load_parent(parent);data=read(out/'carrier.json')
    row=read(out/'selection.json')['rows'][index];cs=[read(out/'conics'/f'{m}.json') for m in data['masks']]
    for c in cs:
        require(sha(out/'conics'/f"{c['mask']}.json")==data['conic_sha256'][str(c['mask'])],'carrier conic changed')
        geo.verify_conic(c,E,points,G)
    t,w1,w2=map(QQ,(row['parameter'],row['w1'],row['w2']))
    ring=E.base_ring().ring();require(w1*w1==ring(cs[0]['q'])(t) and w2*w2==ring(cs[1]['q'])(t),'invalid two-sheet fibre witness')
    job=out/'cases'/f'{index:03d}';packet=certify_seed(parent,t,w1,w2,*cs,job)
    if packet:
        write(job/'request.json',{'parent':'pair-inputs/parent.json','parent_sha256':sha(rt/'pair-inputs/parent.json'),
              'parameter':str(t),'packet':str((job/'seed.json').relative_to(rt)),
              'packet_sha256':sha(job/'seed.json'),'allowance':plan['v3_calls'],'height':plan['height'],'bank_index':index,
              'carrier_sha256':sha(out/'carrier.json'),'selection_sha256':sha(out/'selection.json')})


def self_test():
    """Elkies's published pair is a control, never an X1092 selection input."""
    from sage.all import QQ,PolynomialRing
    import elkies_pair_geometry as geo
    R=PolynomialRing(QQ,'u');K=R.fraction_field();u=R.gen()
    q1=R([289444,38636,4225]);q2=R([22473889,-3269604,54756])
    T=K(q1[0]-u*u)/(130*u-q1[1]);W=65*T+u
    first={'mask':1,'q':geo.polynomial_record(q1),'T':geo.function_record(T),'W':geo.function_record(W)}
    second={'mask':2,'q':geo.polynomial_record(q2)}
    result=geo.carrier(first,second)
    require(result['nontorsion']['rank_lower_bound']==1,'published pair failed positive-rank control')
    print(json.dumps({'status':'PASS_PUBLISHED_PAIR_POINTED_MAPS_AND_NONTORSION','curve':result['curve'],
                      'nontorsion':result['nontorsion']}))


def run_stage(folder,name,phase=None,index=0,seconds=None,backend=False):
    from research_runtime.supervisor import Limits,run as supervise
    plan,rt,out,parent=context(folder);job=out/'phases'/name;receipt=job/'supervisor.json'
    if receipt.exists():
        r=read(receipt)
        if r['outcome']!='completed' or r['returncode']!=0:return False
        require((job/'seal.json').exists(), 'completed stage has no output seal; retain as UNKNOWN')
        seal=read(job/'seal.json')
        require(sha(receipt)==seal['receipt_sha256'], 'supervisor receipt changed')
        for name,digest in seal['files'].items():
            require(sha(out/name)==digest, 'completed stage output changed: '+name)
        return True
    require(not (job/'started.json').exists(),'interrupted unreceipted stage is UNKNOWN; use a fresh campaign, not automatic replay')
    write(job/'started.json',{'phase':phase,'index':index,'time':time.time()})
    script=rt/'elliptic-curves/cas'/('parent_foundry_worker.py' if backend else SELF.name)
    command=[plan['sage'],'-python',str(script)]
    if backend:command+=['--job',str(out/'cases'/f'{index:03d}')]
    else:command+=['_phase','--folder',str(folder),'--phase',phase,'--index',str(index)]
    r=supervise(command,limits=Limits(wall_seconds=seconds or plan['math_seconds'],rss_bytes=plan['rss_bytes']),
                log_path=job/'worker.log',checkpoint_path=receipt,cwd=rt,env=environment(plan['sage']))
    if r['outcome']=='completed' and r['returncode']==0:
        write(job/'seal.json',{'receipt_sha256':sha(receipt),'files':{
            str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*'))
            if p.is_file() and p.relative_to(out).parts[0]!='phases'}})
    print('PAIR_STAGE',name,r['outcome'],r['returncode'],flush=True)
    return r['outcome']=='completed' and r['returncode']==0


def execute(folder):
    folder=folder.resolve();plan,rt,out,parent=context(folder)
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for phase in ('self-test','pool','carrier','score'):
            if (folder/'STOP').exists():return
            if not run_stage(folder,phase,phase):
                write(folder/'REPORT.json',{'status':'UNKNOWN_FAILED_OR_CENSORED_STAGE','phase':phase},False);return
            if phase=='carrier' and not (out/'carrier.json').exists():
                write(folder/'REPORT.json',read(out/'carrier-status.json'),False);return
        rows=read(out/'selection.json')['rows'];results=[]
        for i,row in enumerate(rows):
            if (folder/'STOP').exists():break
            if not run_stage(folder,f'seed-{i:03d}','seed',i,180):
                results.append({'index':i,'status':'UNKNOWN_SEED_STAGE'});continue
            job=out/'cases'/f'{i:03d}';gate=read(job/'gate.json')
            if gate['status']!='PASS_TWO_FINITE_IMPLEMENTATIONS':
                results.append({'index':i,'control':row['control'],**gate});continue
            ok=run_stage(folder,f'v3-{i:03d}',index=i,seconds=plan['job_seconds'],backend=True)
            result={'index':i,'control':row['control'],'initial_rank':19,'rank_lower_bound':19,
                    'status':'UNKNOWN_V3_STAGE','seed_sha256':sha(job/'seed.json')}
            if ok:
                packet=read(job/'packet.json');proof=read(job/'verified.json');r=read(job/'result.json')
                require(proof['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS' and
                        sha(job/'packet.json')==proof['packet_sha256']==r['packet_sha256'],'V3 packet/replay hash differs')
                require(packet['points'][:19]==read(job/'seed.json')['points'],'V3 lost the M19 seed prefix')
                require(r['status']=='PASS_CERTIFIED_PARENT_EVALUATION' and
                        r['request_sha256']==sha(job/'request.json') and r['parameter']==row['parameter'],
                        'V3 request/fibre binding differs')
                require(r['initial_rank']==19 and r['rank_lower_bound']>=19 and r['calls']<=plan['v3_calls'],'V3 rank/budget contract differs')
                result.update(status='PASS_V3_REPLAY',rank_lower_bound=r['rank_lower_bound'],calls=r['calls'],
                              exposure_complete=r['exposure_complete'],packet_sha256=sha(job/'packet.json'))
            results.append(result)
            write(folder/'REPORT.json',{'status':'RUNNING','carrier_kind':read(out/'carrier.json')['kind'],'results':results},False)
            if result['rank_lower_bound']>=32:break
        status=('STOPPED' if (folder/'STOP').exists() else 'TARGET_LOWER_BOUND_REACHED'
                if any(r.get('rank_lower_bound',0)>=32 for r in results) else 'COMPLETE_BOUNDED_CAMPAIGN')
        write(folder/'REPORT.json',{'status':status,
            'carrier_kind':read(out/'carrier.json')['kind'],'results':results,'selected':len(rows),
            'boundary':'Only sealed M19/V3 packets support fibre rank lower bounds. Censored and unattempted cases are not null ranks.'},False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','run','resume','status','self-test','_phase'])
    parser.add_argument('--folder',type=Path,default=ROOT/'artifacts/local/elliptic-curves/x1092-elkies-pairs-v1')
    parser.add_argument('--sage',default='sage')
    for name,default in [('pool-size',64),('pair-limit',128),('addresses',128),('controls',2),('ranked',8),
        ('deep-keep',32),('shallow-bound',997),('deep-bound',32749),('parameter-bits',1024),
        ('v3-calls',100),('job-seconds',1800),('math-seconds',1800)]:
        parser.add_argument('--'+name,type=int,default=default)
    parser.add_argument('--phase',choices=['self-test','pool','carrier','score','seed']);parser.add_argument('--index',type=int,default=0)
    args=parser.parse_args()
    if args.command=='self-test':self_test();return
    if args.command=='_phase':
        if args.phase=='self-test':self_test()
        elif args.phase=='seed':seed_phase(args.folder,args.index)
        else:globals()[args.phase+'_phase'](args.folder)
        return
    if args.command=='status':
        p=args.folder/'REPORT.json';print(json.dumps(read(p) if p.exists() else {'status':'NO_COMPLETED_REPORT'},indent=2));return
    if args.command in ('prepare','run') and not args.folder.exists():prepare(args)
    elif args.command=='prepare':raise ValueError('campaign already exists; use resume')
    if args.command in ('run','resume'):execute(args.folder)


if __name__=='__main__':main()
