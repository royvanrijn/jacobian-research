#!/usr/bin/env python3
"""Eight-prime cost and exact sample gate for a11952 extended projective trace cache."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import benchmark_r17_extended_prime_traces as scalar
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/extended-projective-trace-cache-benchmark-v1';OUT=ART/'extended_projective_trace_cache_benchmark_v1.json'
GATE=ART/'outer_known29_retention_comparison_v1.json'

def sources():return {**scalar.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),spec.ATLAS,GATE,ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed cache benchmark')
    gate=cert.read(GATE)
    if gate['status']!='PASS' or gate['short_score_position']!=435 or gate['extended_score_position']!=2 or gate['same_saved_outer_candidates']!=967:raise ArithmeticError('exact retention-reordering gate required')
    primes=[p for p in _primes_up_to(32749) if p>=4099];cases=[primes[i*(len(primes)-1)//7] for i in range(8)]
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']=='11952');model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.extended-projective-trace-cache-benchmark.v1','sources':sources(),'family':'11952','model':model,'model_hash':digest(model),'primes':cases,'full_selection_prime_roster':primes,'projective_rows':sum(p+1 for p in cases),'full_selection_projective_rows':sum(p+1 for p in primes),'gp_sha256':cert.hashed(scalar.GP),'gp_seconds_per_prime':20,'rss_bytes':536870912,'outer_wall_seconds':300,'direct_residues':'0,floor(p/3),floor(2p/3),p-1,infinity','projected_single_worker_seconds_gate':1800,'gate':'The known native11952 control changes from hypothetical position435 under the short score to2 under extended S1 among the same967 saved outer candidates. Extending many more short-retained survivors with a separate GP call per curve repeats finite-field fibres. Before any broader survivor campaign, measure whether caching each projective residue once per selection prime is affordable. The family is chosen by this exact calibration issue, not generic rank alone.','scope':'Only eight deterministic prime-table cases, spaced across4099through32749, one worker. Compute allp+1 raw homogeneous finite-field fibres in each case with PARI ellcard, retain singular markers, check every prime/residue/Hasse frame and independently verify five complete character sums per prime. Cost extrapolation is descriptive, not a certified runtime bound. No full extended cache, new parameter search, new selector, point search or automatic escalation occurs in this benchmark.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(scalar.GP):raise ArithmeticError('fixed projective-cache inputs changed')
    return p

def evaluate(c,t,p):
    z=0
    for a in reversed(c):z=(z*t+a)%p
    return z

def program(model,p):
    a=[int(v)%p for v in model['A_coefficients_low_to_high']];b=[int(v)%p for v in model['B_coefficients_low_to_high']]
    return f'A={a};B={b};p={p};ev(v,t)={{my(z=0);forstep(j=#v,1,-1,z=(z*t+v[j])%p);z}};gettime();for(i=0,p,a=if(i==p,A[#A],ev(A,i));b=if(i==p,B[#B],ev(B,i));if((4*a^3+27*b^2)%p==0,print("B|",i),print("T|",i,"|",p+1-ellcard(ellinit([Mod(a,p),Mod(b,p)])))));print("MS|",gettime());print("DONE");quit\n'

def table(prime,p,create):
    folder=D/str(prime);path=folder/'raw.json';code=program(p['model'],prime)
    if create:
        if path.exists():raise FileExistsError('preserve prime-table benchmark call')
        c=capture([str(scalar.GP),'-q','-s','256000000'],input_text=code,limits=Limits(p['gp_seconds_per_prime'],p['rss_bytes']),log_path=folder/'gp.log',separate_stderr=True,check=False)
        checkpoint(path,{'program':code,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['program']!=code or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('table call failed or censored')
    lines=raw['stdout'].splitlines();a=[int(v)%prime for v in p['model']['A_coefficients_low_to_high']];b=[int(v)%prime for v in p['model']['B_coefficients_low_to_high']];values=[]
    if len(lines)!=prime+3 or lines[-1]!='DONE' or not lines[-2].startswith('MS|'):raise ArithmeticError('complete projective table framing differs')
    for i,line in enumerate(lines[:-2]):
        aa=a[-1] if i==prime else evaluate(a,i,prime);bb=b[-1] if i==prime else evaluate(b,i,prime);bad=(4*aa**3+27*bb**2)%prime==0;v=line.split('|')
        if len(v)!=(2 if bad else 3) or v[0]!=('B' if bad else 'T') or int(v[1])!=i:raise ArithmeticError('residue/discriminant framing differs')
        value=None if bad else int(v[2])
        if value is not None and value*value>4*prime:raise ArithmeticError('projective Hasse bound differs')
        values.append(value)
    checks=[]
    for i in [0,prime//3,2*prime//3,prime-1,prime]:
        aa=a[-1] if i==prime else evaluate(a,i,prime);bb=b[-1] if i==prime else evaluate(b,i,prime);direct=scalar.direct([0,0,0,aa,bb],prime)
        if direct!=values[i]:raise ArithmeticError('independent table character sum differs')
        checks.append([i,direct])
    saved={'input':{'family':p['family'],'model_hash':p['model_hash'],'prime':prime},'traces':[v or 0 for v in values],'good':[v is not None for v in values]};tp=folder/'table.json'
    if create:
        if tp.exists():raise FileExistsError('preserve benchmark table')
        checkpoint(tp,saved)
    elif cert.read(tp)!=saved:raise ArithmeticError('saved projective table differs')
    return {'prime':prime,'projective_rows':prime+1,'cpu_ms':int(lines[-2][3:]),'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path),'table_sha256':cert.hashed(tp),'direct_checks':checks}

def run(check=False):
    p=protocol()
    if not check and (D/'result.json').exists():raise FileExistsError('preserve fixed eight-prime benchmark')
    data={'status':'RUNNING','protocol_hash':digest(p),'rows':[]}
    if not check:checkpoint(D/'result.json',data)
    for prime in p['primes']:
        data['rows'].append(table(prime,p,not check))
        if not check:checkpoint(D/'result.json',data)
        print('PROJECTIVE TRACE CACHE PRIME',prime,data['rows'][-1]['wall_seconds'],'seconds',flush=True)
    data['status']='PASS';projected=max(r['wall_seconds']/r['projective_rows'] for r in data['rows'])*p['full_selection_projective_rows']
    data.update(projected_single_worker_seconds=projected,cost_gate_passed=projected<=p['projected_single_worker_seconds_gate'],character_sum_checks=40)
    if check:
        if cert.read(D/'result.json')!=data or cert.read(OUT)!={'schema':'elliptic-curves.extended-projective-trace-cache-benchmark-result.v1',**data,'sources':sources(),'claim_boundary':p['scope']}:raise ArithmeticError('projective cache benchmark replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve benchmark certificate')
        checkpoint(D/'result.json',data);checkpoint(OUT,{'schema':'elliptic-curves.extended-projective-trace-cache-benchmark-result.v1',**data,'sources':sources(),'claim_boundary':p['scope']})
    print('PROJECTED SINGLE WORKER',projected,'SECONDS; COST GATE',data['cost_gate_passed'],flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','run','check']);v=a.parse_args();prepare() if v.stage=='prepare' else run(v.stage=='check')
