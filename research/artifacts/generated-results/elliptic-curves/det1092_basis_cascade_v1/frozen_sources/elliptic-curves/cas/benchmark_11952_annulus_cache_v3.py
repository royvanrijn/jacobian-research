#!/usr/bin/env python3
"""Explicit cache-reader height bounds, old-output regression and scalar checks."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import benchmark_retained_extended_cache as old
import benchmark_r17_extended_prime_traces as scalar
import extend_retained_r17_prime_scores as scoring
import extend_outer131072_r17 as models
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
ROOT=old.ROOT;CAS=old.CAS;ART=old.ART;D=ROOT/'artifacts/local/elliptic-curves/11952-annulus-cache-v3'
CPP=CAS/'newfamily/score_retained_projective_cache_v3.cpp';TEST=ROOT/'elliptic-curves/tests/test_retained_projective_cache_scorer_v3.py';BINARY=D/'scorer'
OUT=ART/'annulus_11952_cache_reader_v3.json'
FAILED=ROOT/'artifacts/local/elliptic-curves/11952-new-annulus-scores-v1/cache-raw.json'
PAIRS=[(131072,131073),(-131072,131073),(524287,524288),(-524287,524288),(524288,524287),(-524288,524287),(1,4099*127),(-1,4099*127),(524288,1),(-524288,1)]
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),CPP,TEST,old.OUT,old.TABLE,old.INPUT,old.D/'raw.json',FAILED,models.spec.ATLAS,Path(scalar.__file__),Path(scoring.__file__),Path(models.__file__),CAS/'research_runtime/supervisor.py',CAS/'research_runtime/store.py']}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve cache-reader gate')
    f=next(r for r in cert.read(models.spec.ATLAS)['families'] if r['family']=='11952');proof=cert.read(old.OUT);failure=cert.read(FAILED)
    if proof['status']!='PASS' or cert.hashed(old.TABLE)!=proof['cache_sha256'] or cert.hashed(old.INPUT)!=proof['candidate_sha256']:raise ArithmeticError('prior exact cache fixture differs')
    if failure['supervision']['returncode']!=1 or failure['stdout'] or failure['stderr'].strip()!='invalid primitive candidate':raise ArithmeticError('preserved out-of-domain failure differs')
    rows=[{'parameter':str(cert.F(n,d)),'numerator':n,'denominator':d,'model':models.model_at(f,str(cert.F(n,d)))} for n,d in PAIRS]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.11952-annulus-cache-reader.v3','sources':sources(),'rows':rows,'maximum_height':524288,'compile_seconds':60,'test_seconds':30,'cached_seconds':30,'scalar_seconds_per_row':20,'rss_bytes':2147483648,'gate':'The old reader rejected the already frozen new32768 bank before producing scores because its candidate validator hardcodes131072. Its projective lookup only uses signed residues and invertible denominators. Add an explicit validated height argument without changing quantized lookup, table format, ordering or summation. The default remains131072. Integer safety caps the new argument at2^30; this experiment passes only524288.','scope':'Compile one new binary, run exact signed/infinity and framing unit regressions, compare all967 old fixture scores byte-for-byte with recorded old output, then compare ten fixed boundary parameters against fresh scalar4099..65521 traces and independent sums at4099 and32771. No parameter enumeration, point search, rank claim or implicit retry of the failed frozen pipeline.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('fixed cache-reader gate inputs changed')
    return p

def execute(name,cmd,seconds,check=False,input_text=None):
    p=protocol();path=D/(name+'.json')
    if not check:
        if path.exists():raise FileExistsError('preserve bounded gate invocation')
        c=capture(cmd,input_text=input_text,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'input_text':input_text,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=cmd or raw['input_text']!=input_text or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('fixed cache-reader gate failed/censored')
    return raw

def run(check=False):
    p=protocol()
    execute('compile',['g++','-O3','-std=c++17',str(CPP),'-o',str(BINARY)],60,check)
    test=execute('tests',['/usr/bin/python3',str(TEST)],30,check)
    if 'Ran 3 tests' not in test['stderr'] or not test['stderr'].rstrip().endswith('OK'):raise ArithmeticError('three boundary regressions required')
    original=cert.read(old.D/'raw.json');regression=execute('old967',[str(BINARY),str(old.TABLE),str(old.INPUT)],30,check)
    if regression['stderr'] or regression['stdout']!=original['stdout'] or len(regression['stdout'].splitlines())!=968:raise ArithmeticError('old967 exact score output differs')
    path=D/'boundary-candidates.txt';text='R17-CANDIDATES-V1 10\n'+''.join(f'{n} {d}\n' for n,d in PAIRS)
    if check:
        if path.read_text()!=text:raise ArithmeticError('boundary fixture changed')
    else:
        if path.exists():raise FileExistsError('preserve boundary bank')
        path.write_text(text)
    raw=execute('boundary-cache',[str(BINARY),str(old.TABLE),str(path),'524288'],30,check);lines=raw['stdout'].splitlines();rows=[]
    if raw['stderr'] or len(lines)!=11 or lines[-1]!='S 10 2948':raise ArithmeticError('complete boundary cache frame required')
    for i,r in enumerate(p['rows']):
        code=scalar.program(r['model']);trace=execute('scalar-'+str(i),[str(scalar.GP),'-q','-s','256000000'],20,check,code)
        if trace['stderr']:raise ArithmeticError('scalar boundary stderr')
        values,ms=scalar.parse(trace['stdout'],r['model']);s=scoring.sums(values)
        if lines[i]!=f"R {i} {s['extension_selection_units']} {s['extension_good']}":raise ArithmeticError('new-height cache/scalar disagreement')
        direct=[]
        for q in (4099,32771):
            ap=next(a for prime,a in values if prime==q)
            if scalar.direct(r['model'],q)!=ap:raise ArithmeticError('independent boundary character sum differs')
            direct.append([q,ap])
        rows.append({**r,'scores':s,'direct_checks':direct,'raw_sha256':cert.hashed(D/('scalar-'+str(i)+'.json'))})
    result={'schema':'elliptic-curves.11952-annulus-cache-reader-proof.v3','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'binary_sha256':cert.hashed(BINARY),'old_fixture_rows':967,'boundary_rows':rows,'unit_tests':3,'new_height_bound':524288,'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('cache-reader proof replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve cache-reader proof')
        checkpoint(OUT,result)
    print('ANNULUS CACHE V3: 3 TESTS,967 OLD ROWS,10 NEW BOUNDARY SCALARS PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check']);a=p.parse_args()
    prepare() if a.stage=='prepare' else run(a.stage=='check')
