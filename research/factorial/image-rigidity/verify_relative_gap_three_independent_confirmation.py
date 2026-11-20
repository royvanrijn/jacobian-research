#!/usr/bin/env python3
"""Authenticate an executed COMPLETE finite replay; never launch that job."""
if not __debug__:
    raise SystemExit('Unoptimized Python is required for exact verification.')
from collections import Counter
from hashlib import sha256
from math import comb,factorial,isqrt
from pathlib import Path
import argparse,copy,importlib.util,json,struct

HERE=Path(__file__).resolve().parent
SOURCE_PIN='12d554bd8fe2dc2a933d3be47b6ba791c8c1d689126588be1df59425184e5a38'
RECEIPT_PIN='b5336126eb70bbb25258b215af1a3e3234ee58e6737ddb337a80dc9b4412a865'
PRIMES=(1000033,1000037,1000039,1000003)
FILES={
 'proof':'RELATIVE_GAP_THREE_INDEPENDENT_FINITE_CONFIRMATION.md',
 'checker':'verify_relative_gap_three_independent_confirmation.py',
 'input':'relative_gap_three_independent_confirmation_inputs.json',
 'finite_checker':'independent_verify_relative_gap_three_finite_box.py',
 'finite_receipt':'independent_relative_gap_three_finite_box.json'}
CONTROLS={'literal_integer_coefficient_vectors':54,'wrong_actual_shape_rejected':True,
          'common_factor_detected':True,'corrupted_Bezout_identity_rejected':True}
def digest(p):return sha256(p.read_bytes()).hexdigest()

def validate_input(v):
    assert v=={
      'schema_version':1,'gap_range':[9,63],'minimum_lower_exponent':3,
      'relative_gap_factor':3,'minimum_terminal':2,'absolute_terminal_onset':4000,
      'linear_terminal_factor':128,'support_pairs':532,'terminal_pairs':3694040,
      'maximum_terminal':10751,'maximum_factorial_argument':903084,
      'full_Bezout_attempts':3694045,'primes':list(PRIMES),
      'successful_primes':{'1000033':3694035,'1000037':5},
      'fallbacks':[[5,31,1168],[9,33,2932],[10,50,7079],[5,57,7734],[4,60,7495]],
      'source_sha256':SOURCE_PIN,'completed_receipt_sha256':RECEIPT_PIN}

def pin_sources(packet,ledger):
    assert ledger['schema_version']==1
    records=ledger['files']
    assert len(records)==len(FILES) and {r['role']for r in records}==set(FILES)
    out={}
    for row in records:
        assert row['scope']=='packet' and row['path']==FILES[row['role']]
        path=packet/row['path']
        assert path.is_file() and digest(path)==row['sha256'],'Mandatory source pin mismatch: '+row['role']
        out[row['role']]=path.resolve()
    assert out['checker']==Path(__file__).resolve()
    assert digest(out['finite_checker'])==SOURCE_PIN
    assert digest(out['finite_receipt'])==RECEIPT_PIN
    return out

def domain():
    cases=[(a,d,max(4000,128*(a+d))-1)for d in range(9,64)for a in range(3,d//3+1)]
    assert len(cases)==len(set(cases))==532
    assert sum(last-1 for a,d,last in cases)==3694040
    assert max(last for a,d,last in cases)==10751
    assert max((a+d)*last for a,d,last in cases)==903084
    for p in PRIMES:assert p>903084 and all(p%j for j in range(2,isqrt(p)+1))
    return cases

def verify_completed(r,inp,stream=True):
    assert r['status']=='PASS' and r['complete_finite_box']is True
    assert r['tested_supports']==inp['support_pairs']==532
    assert r['tested_terminal_pairs']==inp['terminal_pairs']==3694040
    assert r['verified_full_Bezout_identities']==inp['full_Bezout_attempts']
    assert r['maximum_factorial_argument_for_full_domain']==903084
    assert r['source_sha256']==SOURCE_PIN
    assert r['primary_prime']==PRIMES[0] and r['alternative_primes']==list(PRIMES[1:])
    assert r['actual_shapes']==['a*N+1','a*(N-1)+1']
    assert r['coefficient_construction']=='Separate binomial/sliding-rising recurrences; preceding actual polynomial reused at next terminal; no factorial tables'
    assert r['versions']=={'python':'3.12.3','python_flint':'0.9.0','FLINT':'3.6.0'}
    assert r['controls']==CONTROLS
    expected=domain();rows=r['rows'];assert len(rows)==len(expected)
    total=attempts_total=stream_cases=0;observed=[];global_counts=Counter();fallback_records=[]
    for row,(a,d,last)in zip(rows,expected):
        assert (row['lower_exponent'],row['gap'],row['last_terminal'])==(a,d,last)
        assert row['first_terminal']==2 and row['tested_cases']==last-1
        assert row['maximum_factorial_argument']==(a+d)*last<=903084
        assert row['maximum_upper_degree']==last and row['maximum_previous_degree']==last-1
        by_terminal={};extra=0;counts=Counter({str(PRIMES[0]):last-1})
        for f in row['fallbacks']:
            N=f['terminal'];attempts=f['attempts']
            assert 2<=N<=last and N not in by_terminal
            assert 2<=len(attempts)<=len(PRIMES)
            assert [p for p,g in attempts]==list(PRIMES[:len(attempts)])
            assert all(isinstance(g,int) and 1<=g<=N-1 for p,g in attempts[:-1])
            assert attempts[-1][1]==0
            chosen=attempts[-1][0];by_terminal[N]=chosen
            counts[str(PRIMES[0])]-=1;counts[str(chosen)]+=1;extra+=len(attempts)-1
            observed.append([a,d,N]);fallback_records.append({'a':a,'d':d,'N':N,'attempts':attempts})
        assert {k:v for k,v in counts.items()if v}==row['successful_primes']
        assert row['verified_full_Bezout_identities']==last-1+extra
        assert sum(row['successful_primes'].values())==last-1
        if stream:
            h=sha256()
            for N in range(2,last+1):h.update(struct.pack('<IIII',a,d,N,by_terminal.get(N,PRIMES[0])))
            assert h.hexdigest()==row['case_stream_sha256']
            stream_cases+=last-1
        total+=last-1;attempts_total+=last-1+extra;global_counts.update(row['successful_primes'])
    assert observed==inp['fallbacks']
    assert total==3694040 and attempts_total==3694045
    assert dict(global_counts)==inp['successful_primes']
    return {'supports':len(rows),'terminal_pairs':total,'terminal_stream_entries':stream_cases,
            'full_Bezout_attempts':attempts_total,'successful_primes':dict(global_counts),
            'fallbacks':fallback_records}

def literal_controls(source):
    assert tuple(source.PRIMES)==PRIMES and source.domain()==domain()
    frozen=source.controls();assert frozen==CONTROLS
    count=0
    for a,d in[(1,1),(2,6),(3,9),(3,63),(4,12),(7,21),(10,30),(20,60),(21,63)]:
        tab=source.make_tables(d,6*(a+d),PRIMES[0])
        for m in range(1,7):
            literal=[comb(m,k)*(factorial(a*m+d*k)//factorial(a*m))%PRIMES[0]for k in range(m+1)]
            assert source.coefficients(a,d,m,PRIMES[0],tab)==literal;count+=1
    assert count==54
    return count,frozen

def factorial_pair(a,d,N,p):
    # Different from the executed inverse/sliding-rising construction.
    from flint import nmod_poly
    limit=(a+d)*N;assert limit<p
    f=[1]*(limit+1)
    for j in range(1,limit+1):f[j]=f[j-1]*j%p
    inv=[1]*(limit+1);inv[-1]=pow(f[-1],-1,p)
    for j in range(limit,0,-1):inv[j-1]=inv[j]*j%p
    def poly(m):
        pref=f[m]*inv[a*m]%p
        cs=[pref*f[a*m+d*k]*inv[k]*inv[m-k]%p for k in range(m+1)]
        q=nmod_poly(cs,p);assert q.degree()==m and q[0]==1 and all(cs)
        return q
    return poly(N),poly(N-1)

def replay_fallbacks(records):
    from flint import nmod_poly
    out=[]
    for row in records:
        for p,degree in row['attempts']:
            hi,lo=factorial_pair(row['a'],row['d'],row['N'],p)
            q,s,t=hi.xgcd(lo)
            assert q.degree()==degree and s*hi+t*lo==q
            if degree==0:assert q==nmod_poly([1],p)
            out.append({'a':row['a'],'d':row['d'],'N':row['N'],'prime':p,
                        'gcd_degree':degree,'full_Bezout_identity_verified':True})
    assert len(out)==10
    return out

def negative_controls(r,inp):
    negatives=0
    for label in['partial','missing_support','case_count','attempt_count','shape','source','fallback','stream']:
        bad=copy.deepcopy(r)
        if label=='partial':bad['complete_finite_box']=False
        elif label=='missing_support':bad['rows'].pop()
        elif label=='case_count':bad['tested_terminal_pairs']-=1
        elif label=='attempt_count':bad['verified_full_Bezout_identities']-=1
        elif label=='shape':bad['actual_shapes']=['a*N+1','a*N+1']
        elif label=='source':bad['source_sha256']='0'*64
        elif label=='fallback':next(q for q in bad['rows']if q['fallbacks'])['fallbacks'][0]['attempts'][-1][0]=PRIMES[0]
        else:bad['rows'][0]['case_stream_sha256']='0'*64
        try:verify_completed(bad,inp,stream=(label=='stream'))
        except AssertionError:negatives+=1
        else:raise AssertionError('Corrupted artifact accepted: '+label)
    assert negatives==8
    return negatives

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet-root',type=Path,default=HERE)
    ap.add_argument('--output',type=Path,default=HERE/'relative_gap_three_independent_confirmation_verification.json')
    args=ap.parse_args();packet=args.packet_root.resolve()
    ledger_path=packet/'relative_gap_three_independent_confirmation_source_ledger.json'
    ledger=json.loads(ledger_path.read_text());paths=pin_sources(packet,ledger)
    inp=json.loads(paths['input'].read_text());validate_input(inp)
    raw=paths['finite_receipt'].read_bytes();completed=json.loads(raw)
    counts=verify_completed(completed,inp)
    spec=importlib.util.spec_from_file_location('frozen_independent_actual_source',paths['finite_checker'])
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    literal,frozen_controls=literal_controls(source)
    fallback_attempts=replay_fallbacks(counts['fallbacks'])
    negatives=negative_controls(completed,inp)
    assert paths['finite_receipt'].read_bytes()==raw and digest(paths['finite_checker'])==SOURCE_PIN
    result={'schema_version':1,'status':'PASS','scope':'Authenticated completed FINITE actual domain only. The all-order theorem imports separate analytical premises. No full-domain restart, external review or formal verification.',
      'source_ledger_sha256':digest(ledger_path),'actual_executed_source_sha256':SOURCE_PIN,
      'actual_completed_receipt_sha256':RECEIPT_PIN,'completed_supports':counts['supports'],
      'completed_terminal_pairs':counts['terminal_pairs'],'reconstructed_terminal_stream_entries':counts['terminal_stream_entries'],
      'verified_full_Bezout_attempt_count':counts['full_Bezout_attempts'],
      'successful_primes':counts['successful_primes'],'fresh_literal_coefficient_vectors':literal,
      'frozen_source_controls_freshly_replayed':frozen_controls,
      'fresh_full_fallback_attempt_replays':fallback_attempts,'artifact_negative_controls':negatives,
      'full_domain_restarted':False,'assurance_scope':'Additional internal, differently implemented finite confirmation only.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','supports':counts['supports'],'terminal_pairs':counts['terminal_pairs'],
      'full_Bezout_attempts':counts['full_Bezout_attempts'],'fresh_literals':literal,
      'fresh_fallback_attempts':len(fallback_attempts),'artifact_negatives':negatives}))
if __name__=='__main__':main()
