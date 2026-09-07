"""Optional fixed-map PARI backend with exact shared pointed-chart witnesses.

This is a separate backend, never a synthetic GMP transcript. A completed
finite-box claim trusts the pinned PARI execution; replay checks every retained
square and rational map without invoking either search engine.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import subprocess,time
from pointed_quartic_search import CoordinatePolicy,point_record,sources as pointed_sources
from search_observability import transform,multiply,prepare_chart
ROOT=Path(__file__).resolve().parents[2]


def sources():
    paths=[Path(__file__).resolve(),Path(__file__).with_name('search_observability.py')]
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}}


def validate_map(search,mapping):
    x,y=search.to_short(search.centre);A=search.short_model[3]
    raw=(-3*x*x-4*A,-8*y,-6*x,F(0),F(1))
    if list(map(str,raw))!=mapping['raw_coefficients']:raise ArithmeticError('raw pointed quartic differs')
    M=tuple(map(F,mapping['matrix']));first=tuple(map(F,mapping['first_matrix']));second=tuple(map(F,mapping['second_matrix']))
    if M!=multiply(first,second) or mapping['coordinate_policy']!={'kind':'raw','matrix':mapping['matrix']}:raise ArithmeticError('horizontal maps differ')
    if search.policy.record()!=CoordinatePolicy.parse(mapping['coordinate_policy']).record():raise ArithmeticError('search coordinate differs')
    P,Q=[list(map(F,mapping[k])) for k in ('reduced_P','reduced_Q')]
    if len(P)!=5 or len(Q)!=3:raise ArithmeticError('reduced polynomial degrees differ')
    disc=tuple(4*P[j]+sum(Q[k]*Q[j-k] for k in range(3) if 0<=j-k<3) for j in range(5));ratio=F(mapping['square_ratio'])
    if list(map(str,disc))!=mapping['discriminant_quartic'] or ratio<=0 or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator or transform(raw,M)!=tuple(ratio*v for v in disc):raise ArithmeticError('reduced quartic identity failed')
    prepare_chart(search.chart_record())
    return P,Q


def program(mapping,height):
    if type(height) is not int or not 1<=height<=1000000:raise ValueError('invalid finite box')
    poly=lambda a:'+'.join(f'({F(v)})*x^{j}' for j,v in enumerate(a))
    return 'C=['+poly(mapping['reduced_P'])+','+poly(mapping['reduced_Q'])+'];gettime();R=hyperellratpoints(C,'+str(height)+');print("MS|",gettime());for(i=1,#R,print("X|",R[i][1]));print("DONE|",#R);quit\n'


def witnesses(search,stdout,status,height):
    hits=[];found=set();milliseconds=None
    if status=='bounded_search_complete':
        lines=stdout.splitlines();xlines=[s[2:] for s in lines if s.startswith('X|')];done=[s[5:] for s in lines if s.startswith('DONE|')];ms=[s[3:] for s in lines if s.startswith('MS|')]
        if len(done)!=1 or len(ms)!=1 or int(done[0])!=len(xlines):raise ArithmeticError('incomplete PARI output framing')
        milliseconds=int(ms[0]);xs=sorted(set(map(F,xlines)))
        for t in xs:
            n,d=t.numerator,t.denominator
            if max(abs(n),d)>height:raise ArithmeticError('PARI hit outside fixed box')
            f=sum(v*n**j*d**(4-j) for j,v in enumerate(search.coefficients))
            if f<0 or isqrt(f)**2!=f:raise ArithmeticError('returned coordinate is not an exact square')
            hits.append((n,d,isqrt(f)))
    elif status not in ('bounded_search_timeout','backend_failure'):raise ArithmeticError('unknown backend status')
    lead=search.coefficients[4]
    if lead>=0 and isqrt(lead)**2==lead:hits.append((1,0,isqrt(lead)))
    for n,d,root in hits:
        for sign in {root,-root}:
            p=search.map_hit(n,d,sign)
            if p is not None:found.add(p)
    return hits,tuple(sorted(found)),milliseconds


def execute(search,mapping,height,seconds,expected_gp_hash):
    validate_map(search,mapping);gp=Path('/usr/bin/gp');actual=sha256(gp.read_bytes()).hexdigest()
    if actual!=expected_gp_hash:raise ArithmeticError('PARI binary changed')
    command=program(mapping,height);started=time.monotonic()
    try:
        p=subprocess.run([str(gp),'-q','-s','256000000'],input=command,text=True,capture_output=True,timeout=seconds)
        stdout,stderr,returncode=p.stdout,p.stderr,p.returncode
        status='bounded_search_complete' if returncode==0 and '***' not in stderr else 'backend_failure'
    except subprocess.TimeoutExpired as e:
        decode=lambda x:x.decode('utf8',errors='replace') if isinstance(x,bytes) else (x or '')
        stdout,stderr,returncode=decode(e.stdout),decode(e.stderr),None;status='bounded_search_timeout'
    wall=time.monotonic()-started;hits,points,ms=witnesses(search,stdout,status,height)
    record={**search.chart_record(),'backend':'pari_fixed_pointed_v1','source_hashes':sources(),'mw_state_key':search.state.key,'height_bound':height,'timeout_seconds':seconds,'status':status,'gp_binary_sha256':actual,'program':command,'stdout':stdout,'stderr':stderr,'returncode':returncode,'wall_seconds':wall,'search_cpu_ms':ms,'primitive_square_hits':[[str(v) for v in h] for h in hits],'finite_curve_points':[point_record(p) for p in points],'infinity_checked':True,'claim_boundary':'Retained exact square and curve-map witnesses. Completed box coverage trusts the pinned PARI invocation; incomplete invocations have no denominator-prefix claim. Infinity is checked exactly in Python. No rank upper bound.'}
    return record,points


def replay(search,mapping,record):
    validate_map(search,mapping)
    if record['backend']!='pari_fixed_pointed_v1' or record['source_hashes']!=sources() or record['mw_state_key']!=search.state.key:raise ArithmeticError('backend or state binding changed')
    if any(record[k]!=v for k,v in search.chart_record().items()):raise ArithmeticError('pointed chart changed')
    if record['program']!=program(mapping,record['height_bound']) or not record['infinity_checked']:raise ArithmeticError('fixed PARI program changed')
    if record['status']=='bounded_search_complete' and (record['returncode']!=0 or '***' in record['stderr']):raise ArithmeticError('failed PARI call asserted complete')
    hits,points,ms=witnesses(search,record['stdout'],record['status'],record['height_bound'])
    if record['primitive_square_hits']!=[[str(v) for v in h] for h in hits] or record['finite_curve_points']!=[point_record(p) for p in points] or record['search_cpu_ms']!=ms:raise ArithmeticError('PARI square/point witnesses changed')
    return points
