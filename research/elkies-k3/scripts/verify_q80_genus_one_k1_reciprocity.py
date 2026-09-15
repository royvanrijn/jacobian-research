#!/usr/bin/env python3
"""Independent composite-code joins and generic C++ degree-three replay for Q80 k1."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path
import resource
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
BOUNDARY='artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
PRODUCER='elkies-k3/scripts/certify_q80_genus_one_k1_reciprocity.sage'
CPP='elkies-k3/scripts/certify_q80_cubic_norm_census.cpp'
REPLAY='elkies-k3/scripts/replay_q80_cubic_norm_census.cpp'
TEST='tests/test_q80_genus_one_k1_reciprocity.py'
P=131
def require(b,message):
    if not b:raise ValueError(message)
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a
def coeff(a):return trim([Fraction(v).numerator*pow(Fraction(v).denominator,-1,P)%P for v in a])
def add(a,b):return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%P for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([v*c%P for v in a])
def mul(a,b):
    if not a or not b:return []
    c=[0]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):c[i+j]+=v*w
    return trim([v%P for v in c])
def divrem(a,b):
    a=trim(a);b=trim(b);require(bool(b),'nonzero polynomial divisor')
    q=[0]*max(0,len(a)-len(b)+1);bi=pow(b[-1],-1,P)
    while len(a)>=len(b):
        k=len(a)-len(b);v=a[-1]*bi%P;q[k]=v
        for j,w in enumerate(b):a[k+j]=(a[k+j]-v*w)%P
        a=trim(a)
    return trim(q),a
def gcd(a,b):
    while b:a,b=b,divrem(a,b)[1]
    return scale(a,pow(a[-1],-1,P)) if a else []
def power_mod(a,n,f):
    r=[1]
    while n:
        if n&1:r=divrem(mul(r,a),f)[1]
        a=divrem(mul(a,a),f)[1];n>>=1
    return r
def ev(a,t):
    v=0
    for c in reversed(a):v=(v*t+c)%P
    return v
def check_fibres(record,A,B):
    delta=add(scale(mul(mul(A,A),A),4),scale(mul(B,B),27));prod=[1]
    factors=record['discriminant_factors']
    require([len(f)-1 for f,e in factors]==[2,3,19],'complete bad-fibre degree partition')
    for f,e in factors:
        require(e==1,'simple discriminant factors')
        n=len(f)-1;xp=power_mod([0,1],P,f)
        require(len(gcd(add(xp,[0,P-1]),f))==1,'no degree-one subfactor')
        h=[0,1]
        for _ in range(n):h=power_mod(h,P,f)
        require(h==[0,1],'prime-degree irreducibility')
        prod=mul(prod,f)
    require(prod==scale(delta,pow(delta[-1],-1,P)),'exact complete discriminant factorization')
    require(all(ev(delta,t) for t in range(P)) and (4*A[8]**3+27*B[12]**2)%P,'all rational projective fibres smooth')
    require(all((t**3+t+3)%P for t in range(P)),'irreducible cubic field modulus')
    require(record['all_projective_rational_fibres_smooth'] and record['cubic_modulus_irreducible'],'fibre and field assertions')
def rational_data(rat):
    agreement={}
    for r in rat:
        for e,c in zip(r['rational_roots'],r['rational_codes']):
            require(c!=0 and c not in agreement,'distinct nonzero rational agreement codes')
            agreement[c]=(r['t'],e)
    require(len(agreement)==110,'all rational agreement root points')
    split=[r for r in rat if len(r['rational_roots'])==3]
    require(len(split)==16,'full rational splitting roster')
    return agreement,split
def permitted_agreement(codes,agreement,forbidden_bases):
    require(codes[0]^codes[1]^codes[2]==0,'cubic-root character product')
    if codes[2] or codes[0] not in agreement:return None
    point=agreement[codes[0]]
    return point if point[0] not in forbidden_bases else None
def check_composite(record,rat,quad):
    agreement,split=rational_data(rat);ropts=[];rindex=defaultdict(list)
    for i,r in enumerate(split):
        for option in itertools.permutations(zip(r['rational_roots'],r['rational_codes'])):
            codes=[c for e,c in option];row=(i,r['t'],codes);ropts.append(row);rindex[codes[2]].append(row)
    counts={'rational_triple_choices':len(list(itertools.combinations(split,3)))*6**3};candidates=[];unused=0
    for i,a,ca in ropts:
        for j,b,cb in ropts:
            if i>=j:continue
            for k,c,cc in rindex[ca[2]^cb[2]]:
                if j>=k:continue
                unused+=1;codes=[ca[n]^cb[n]^cc[n] for n in range(3)]
                point=permitted_agreement(codes,agreement,[a,b,c])
                if point is not None:candidates.append({'type':'111','agreement':list(point),'d_bases':[a,b,c]})
    if unused:counts['rational_triples_unused_zero']=unused
    qindex=defaultdict(list);qtotal=0
    for row in quad:
        if len(row['roots'])!=3:continue
        for option in set(itertools.permutations(zip(row['roots'],row['norm_codes']))):
            codes=[c for e,c in option];qindex[codes[2]].append((row['t'],codes));qtotal+=1
    counts['rational_quadratic_choices']=len(ropts)*qtotal;unused=0
    for i,t,ca in ropts:
        for u,cb in qindex[ca[2]]:
            unused+=1;codes=[ca[n]^cb[n] for n in range(3)]
            point=permitted_agreement(codes,agreement,[t])
            if point is not None:candidates.append({'type':'12','agreement':list(point),'d_bases':[t,u]})
    if unused:counts['rational_quadratic_unused_zero']=unused
    require(counts['rational_triple_choices']==120960 and counts['rational_quadratic_choices']==838944,'complete composite degree ranges')
    require(not candidates,'no coupled composite norm solution')
    expected={'counts':counts,'candidate_rows':[],'agreement_points':110,'full_split_rational_sites':[r['t'] for r in split]}
    require(record==expected,'independent norm-code joins reproduce complete counters')
    return agreement,counts
def make_input(source,agreement):
    model=source['weierstrass_model'];A=coeff(model['A_coefficients_low_to_high']);B=coeff(model['B_coefficients_low_to_high'])
    lines=[]
    def emit(a):lines.append(str(len(a))+' '+' '.join(map(str,a)))
    emit(A);emit(B);lines.append('17')
    require(len(source['sections']['records'])==17,'seventeen literal source sections')
    for r in source['sections']['records']:
        x=r['X'];n=coeff(x['numerator_coefficients_low_to_high']);d=coeff(x['denominator_coefficients_low_to_high'])
        g=gcd(n,d);n,rn=divrem(n,g);d,rd=divrem(d,g);require(not rn and not rd,'exact rational coordinate cancellation')
        emit(n);emit(d)
    lines.append('110')
    for c,(t,e) in agreement.items():lines.append(str(c)+' '+str(-2 if t is None else t))
    return '\n'.join(lines)+'\n'
def check_cubic_result(result,agreement):
    require(result['status']=='COMPLETE' and result['orbits']==749320,'complete cubic orbit census')
    require((result['smooth_split'],result['nodal'])==(125595,1),'all split cubic and nodal orbits')
    require(result['max_factor_trials']<=64,'factorization stayed within its declared bound')
    rows=result['zero_norm_rows'];require(len(rows)==8 and len({r['t'] for r in rows})==8,'complete distinct zero-norm rows')
    require(sum(not r['smooth'] for r in rows)==1,'nodal zero-norm row retained')
    for row in rows:
        require(0 in row['norm_codes'] and row['norm_codes'][0]^row['norm_codes'][1]^row['norm_codes'][2]==0,'unused-root and product norm conditions')
        require(row['agreement_code']==-1,'no rational agreement code')
        for j,c in enumerate(row['norm_codes']):
            if c==0:require(row['norm_codes'][(j+1)%3] not in agreement,'zero unused norm does not supply agreement')
    require(result.get('candidate_rows',[])==[] and result.get('candidates',0)==0,'no cubic norm candidates')
def verify(out):
    started=time.process_time();wall=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    paths={SOURCE,PRODUCER,CPP,REPLAY,TEST,str(Path(__file__).relative_to(ROOT)),
           *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
           *[BOUNDARY+'/'+n for n in ['input.json','result.json','character-gate.json','independent-replay.json']]}
    require(set(packet['bindings'])==paths,'exact frozen binding set')
    for p,h in {**packet['bindings'],**packet['preserved_preflight']}.items():require(digest(ROOT/p)==h,'source/preflight binding '+p)
    require((packet['prime'],packet['cubic_orbits'],packet['factorization_trial_limit'])==(131,749320,64),'frozen finite scope')
    require(packet['cubic_modulus_low_to_high']==[3,1,0,1],'cubic field presentation')
    require((packet['driver_cpu_seconds'],packet['census_cpu_seconds'],packet['driver_memory_bytes'],packet['census_memory_bytes'])==(30,90,2*1024**3,512*1024**2),'declared process limits')
    require(result['status']=='PASS' and result['input_sha256']==digest(out/'input.json'),'bound result')
    require(set(result['records'])=={'fibre-hypotheses.json','composite-norms.json','cubic-input.txt','cubic-census.json','cubic-row-checks.json'},'complete finite records')
    for p,h in result['records'].items():require(digest(out/p)==h,'result binding '+p)
    source=read(ROOT/SOURCE);A=coeff(source['weierstrass_model']['A_coefficients_low_to_high']);B=coeff(source['weierstrass_model']['B_coefficients_low_to_high'])
    check_fibres(read(out/'fibre-hypotheses.json'),A,B)
    rat=read(ROOT/DEGREE/'rational-fibres.json')['rows'];quad=read(ROOT/DEGREE/'quadratic-fibres.json')['rows']
    agreement,counts=check_composite(read(out/'composite-norms.json'),rat,quad)
    input_text=make_input(source,agreement);require(input_text==(out/'cubic-input.txt').read_text(),'literal independently reduced coefficient input')
    cubic=read(out/'cubic-census.json');check_cubic_result(cubic,agreement)
    with tempfile.TemporaryDirectory(prefix='q80-k1-independent-') as tmp:
        exe=Path(tmp)/'replay'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/REPLAY),'-o',str(exe)],check=True,timeout=30)
        completed=subprocess.run([str(exe)],input=input_text,text=True,capture_output=True,check=True,timeout=95)
    replay=json.loads(completed.stdout);check_cubic_result(replay,agreement)
    for key in ['status','orbits','smooth_split','nodal','zero_norm_rows']:require(cubic[key]==replay[key],'independent cubic census '+key)
    row_check=read(out/'cubic-row-checks.json')
    require(row_check=={'sage_checked_rows':[{'t':r['t'],'norm_codes':r['norm_codes']} for r in cubic['zero_norm_rows']]},'all special row cross-checks')
    require(result['finite_candidate_count']==0 and result['written_integrality_argument_required'],'finite gate versus written integrality scope')
    require(result['excluded_allocations']==['k1,j0','k1,j1','k1,j2','k1,j3'] and result['genus1_allocations_remaining']==5,'exact chart allocation scope')
    require(not result['positive_mw17_target_complete'] and not result['formal_verification'] and not result['old_norm8_replay_upgraded'],'positive endpoint and assurance boundary')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
            'checker_sha256':digest(Path(__file__)),'independent_census_source_sha256':digest(ROOT/REPLAY),
            'composite_counts':counts,'cubic_orbits':749320,'smooth_split_cubic_orbits':125595,
            'nodal_cubic_orbits':1,'zero_norm_rows':8,'coupled_norm_candidates':0,
            'python_cpu_seconds':time.process_time()-started,'census_cpu_seconds':replay['cpu_seconds'],
            'wall_seconds':time.monotonic()-wall,'census_factor_trials':replay['max_factor_trials'],
            'driver_cpu_limit':30,'census_cpu_limit':90,'driver_memory_bytes':2*1024**3,'census_memory_bytes':512*1024**2,
            'written_integrality_argument_required':True,'prior_complete_censuses_rerun':False}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,default=DEFAULT)
    parser.add_argument('--receipt',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,95));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    record=verify(args.directory)
    if args.receipt:
        with args.receipt.open('x') as f:json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(record,sort_keys=True),flush=True)
