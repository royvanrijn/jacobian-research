#!/usr/bin/env python3
"""Exact same-population short/extended rank contrast for the known29 control."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_outer_known29_score as audit
from research_runtime.store import checkpoint
ROOT=audit.ROOT;ART=audit.ART;OUT=ART/'outer_known29_retention_comparison_v1.json'

def expected():
    proof=cert.read(audit.OUT)
    if proof!=json.loads(json.dumps(audit.expected())):raise ArithmeticError('exact control scores differ')
    parent=audit.extension.parent;p=parent.protocol();saved=cert.read(audit.extension.D/'result.json');c=proof['control'];rows=[r for r in saved['rows'] if r['family']=='11952']
    shortkey=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator'])
    key=shortkey({'score_units':c['short_score_units'],'good_primes':c['short_good_primes'],'denominator':c['denominator'],'numerator':c['numerator']})
    position=1+sum(shortkey(r)<key for r in rows);positive=next(r for r in p['rows'] if r['family']=='11952' and r['sign']==1)
    sources=[Path(__file__).resolve(),audit.OUT,Path(audit.__file__).resolve(),parent.D/'protocol.json',audit.extension.D/'result.json']
    return {'schema':'elliptic-curves.outer-known29-retention-comparison.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sources},'family':'11952','parameter':'89074/31895','same_saved_outer_candidates':len(rows),'short_score_position':position,'extended_score_position':proof['hypothetical_extended_score_position'],'selected_positive_denominator_shard':positive['shard'],'control_denominator_shard':(c['denominator']-1)%1024,'denominator_modulus':1024,'claim_boundary':'Retrospective exact ordering of one known29 control against the same967 saved outer candidates. The control is not in the sampled denominator residue; these positions are hypothetical insertion ranks conditional on the retained candidate population, not full parameter-population quantiles. Short and extended S1 rankings differ substantially on this exact case. No prospective sensitivity, rank density, global discard probability or new selector is proved. The completed48-curve point trial and its score/point inputs remain unchanged.'}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');v=a.parse_args();r=expected()
    if v.check:
        if cert.read(OUT)!=r:raise ArithmeticError('same-population retention contrast differs')
    else:
        if OUT.exists():raise FileExistsError('preserve retention comparison')
        checkpoint(OUT,r)
    print('SAME967 POPULATION:KNOWN29 SHORT POSITION',r['short_score_position'],'EXTENDED',r['extended_score_position'],flush=True)
