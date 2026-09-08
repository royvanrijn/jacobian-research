"""Build/replay a bounded reusable trace table and a fixed timing comparison.

Sixteen first good presentation primes <=397; one worker, 120 seconds.
Reuse the shared exact character-sum engine and independently check every
table entry with PARI. No fibre search, parameter promotion or rank claim.
"""
import argparse,json,runpy,signal,sys,time
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prime_range,pari

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
BASE=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=BASE/'curve302_parent_trace_tables_v1.json'
REPORT=ROOT/'artifacts/local/elliptic-curves/curve302-search-preparation/trace-benchmark.json'
PARENT=BASE/'curve302_recovered_mw17_parent_v1.json'
PROOF=BASE/'curve302_recovered_mw17_parent_proof_v1.json'
sys.path.insert(0,str(CAS))
from curve302_parent_trace_cache import traces_at,score,legacy_text
from research_runtime.finite_fields import family_traces
from research_runtime.store import FiniteFieldFacts
from research_runtime.memory_store import MemoryFactStore


def build():
    start=time.perf_counter()
    proof=json.loads(PROOF.read_text())
    assert proof['status']=='PASS_FULL_ARITHMETIC_MW17_PARENT'
    assert proof['generic_arithmetic_MW_rank']==17 and proof['rational_NS_rank']==19
    assert proof['input_sha256'][str(PARENT.relative_to(ROOT))]==sha256(PARENT.read_bytes()).hexdigest()
    E,_,_=runpy.run_path(str(CAS/'load_curve302_recovered_parent.sage'))['load_curve302_recovered_parent']()
    R=E.base_ring().ring()
    A,B=R(-E.c4()/48),R(-E.c6()/864)
    assert A.degree()==8 and B.degree()==12
    tables=[];skipped=[]
    facts=FiniteFieldFacts(MemoryFactStore())
    for pp in prime_range(5,398):
        p=int(pp);Rp=PolynomialRing(GF(p),'t')
        try:a,b=Rp(A),Rp(B)
        except (ZeroDivisionError,TypeError):
            skipped.append({'prime':p,'reason':'coefficient denominator'});continue
        delta=-16*(4*a**3+27*b**2)
        if delta.degree()!=24 or not delta.is_squarefree():
            skipped.append({'prime':p,'reason':'presentation lacks squarefree degree24 reduction'});continue
        raw=family_traces(list(A),list(B),p,facts=facts)
        table={'prime':p,'band':'discovery' if len(tables)<8 else 'held_out',
               'A':list(map(int,a.list())),'B':list(map(int,b.list())),
               'traces':[row['trace'] for row in raw['fibres']]}
        tables.append(table)
        if len(tables)==16:break
    assert len(tables)==16
    built=time.perf_counter()-start
    # Independent implementation: PARI cardinality, including the infinity
    # fibre in the degree(8,12) Weierstrass chart. Unknown singular entries
    # are checked against their discriminant, not assigned fake traces.
    for row in tables:
        p=row['prime'];a,b=map(PolynomialRing(GF(p),'t'),[row['A'],row['B']])
        for residue,trace in enumerate(row['traces']):
            aa,bb=(a(residue),b(residue)) if residue<p else (a[8],b[12])
            if not (4*aa**3+27*bb**2):assert trace is None
            else:assert trace==p+1-int(pari.ellinit([0,0,0,int(aa),int(bb)],p).ellcard())
    paths=[Path(__file__),CAS/'curve302_parent_trace_cache.py',CAS/'research_runtime/finite_fields.py',PARENT,PROOF]
    result={'schema':'curve302.projective-traces.v1','parent_path':str(PARENT.relative_to(ROOT)),
            'parent_sha256':sha256(PARENT.read_bytes()).hexdigest(),
            'sources':{str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in paths},
            'capacity_gate':'certified arithmetic MW17, rational Picard19; same surface, no new candidate',
            'selection':'first16 primes>=5 and<=397 with integral coefficients and squarefree degree24 discriminant',
            'skipped_primes':skipped,'tables':tables,'independently_verified_entries':sum(len(row['traces']) for row in tables),
            'infinity_index':'p','singular_entry':None,'rank_claim':None,
            'score_boundary':'Scores are scheduling heuristics; missing entries and low scores exclude no rational fibre rank.'}
    # Timing workload is fixed in advance and never used to promote fibres.
    # Coprime parameters with |numerator|,denominator<=32, excluding zero;
    # first512 in deterministic order. Add all sixteen bad-denominator
    # controls and projective infinity. Direct comparison uses the same
    # character-sum algorithm without amortized family precomputation.
    parameters=sorted({Fraction(m,n) for n in range(1,33) for m in range(-32,33) if m},key=lambda q:(q.denominator, q.numerator))[:512]
    parameters += [Fraction(1,row['prime']) for row in tables]+['infinity']
    before=time.perf_counter();cached=[traces_at(result,q) for q in parameters];lookup=time.perf_counter()-before
    def ev(coeff,t,p):
        value=0
        for c in reversed(coeff):value=(value*t+c)%p
        return value
    chars={}
    for row in tables:
        p=row['prime'];ch=[-1]*p;ch[0]=0
        for x in range(1,p):ch[x*x%p]=1
        chars[p]=ch
    before=time.perf_counter();direct=[]
    for q in parameters:
        rowvals=[]
        for row in tables:
            p=row['prime'];den=0 if q=='infinity' else q.denominator%p
            v=p if not den else q.numerator*pow(den,-1,p)%p
            aa,bb=(row['A'][8],row['B'][12]) if v==p else (ev(row['A'],v,p),ev(row['B'],v,p))
            rowvals.append(None if (4*aa**3+27*bb**2)%p==0 else -sum(chars[p][(x*x*x+aa*x+bb)%p] for x in range(p)))
        direct.append(rowvals)
    direct_seconds=time.perf_counter()-before
    assert direct==cached
    report={'status':'PASS','parameters':len(parameters),'parameter_sha256':sha256(json.dumps(list(map(str,parameters))).encode()).hexdigest(),
            'prime_count':16,'projective_entries':result['independently_verified_entries'],
            'cold_build_seconds':built,'warm_lookup_seconds':lookup,'uncached_character_sum_seconds':direct_seconds,
            'warm_lookup_speedup':direct_seconds/lookup,'total_seconds':time.perf_counter()-start,
            'boundary':'Single observed timing; excludes file load, not a point-search speedup or rank-success predictor. No parameter selection or promotion.'}
    REPORT.parent.mkdir(parents=True,exist_ok=True);REPORT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args();signal.alarm(120)
    result=build()
    if args.build:
        assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n')
    assert result==json.loads(OUT.read_text())
    exported=OUT.with_suffix('.txt')
    if args.build:
        assert not exported.exists();exported.write_text(legacy_text(result))
    assert exported.read_text()==legacy_text(result)
    print('PASS exact cached traces, independent PARI replay, finite and infinity parameter routing')
