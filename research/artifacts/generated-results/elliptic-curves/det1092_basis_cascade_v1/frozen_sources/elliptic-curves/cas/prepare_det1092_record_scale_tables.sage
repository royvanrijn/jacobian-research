#!/usr/bin/env sage-python
"""Complete projective traces <=997, independently checked at every residue."""
import sys
from pathlib import Path
from sage.all import pari
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import det1092_record_scale_selection as campaign
from research_runtime.finite_fields import family_traces
from research_runtime.store import FiniteFieldFacts,checkpoint
from research_runtime.memory_store import MemoryFactStore
def main():
    p=campaign.protocol();A,B=campaign.coefficients();out=campaign.D/'tables.json';assert not out.exists()
    tables=[];facts=FiniteFieldFacts(MemoryFactStore())
    for q in p['first_primes']:
        raw=family_traces(A,B,q,facts=facts);rows=raw['fibres']
        for r in rows:
            a,b=r['a'],r['b']
            if r['singular']:assert r['trace'] is None and (4*a**3+27*b*b)%q==0
            else:assert r['trace']==q+1-int(pari.ellinit([0,0,0,a,b],q).ellcard())
        tables.append(dict(prime=q,traces=[r['trace'] for r in rows],ambiguous=[r['a']==r['b']==0 for r in rows]))
        checkpoint(campaign.D/'table-progress.json',dict(status='RUNNING',completed_primes=len(tables)))
    # Exact regression for the historical lost-good-prime branch.
    for q in [5,13]:
        for k in [1,2]:
            row={'model':['0','0','0',str(q**(4*k)),str(q**(6*k))]}
            a,b,scale,bad=campaign.scalar.reduction(row,q)
            assert (a,b,scale,bad)==(1,1,k,False)
            assert campaign.scalar.direct(a,b,q)==q+1-int(pari.ellinit([0,0,0,a,b],q).ellcard())
    checkpoint(out,dict(status='PASS',protocol_sha256=campaign.cert.hashed(campaign.D/'protocol.json'),
        source_sha256=campaign.cert.hashed(Path(__file__)),primes=p['first_primes'],tables=tables,
        independently_checked_entries=sum(len(t['traces']) for t in tables),
        scope='Exact projective smooth traces. Singular A=B=0 cells are flagged for per-parameter p4/p6 minimization; they are not assigned a permanently omitted good-prime contribution. Independent PARI replay at every residue, including infinity.'))
    print('PASS',len(tables),'complete independently checked prime tables',flush=True)
if __name__=='__main__':main()
