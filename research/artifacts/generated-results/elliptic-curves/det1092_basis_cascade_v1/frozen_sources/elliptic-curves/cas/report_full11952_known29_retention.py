#!/usr/bin/env python3
"""Post-freeze control retention audit; never changes the prospective cohort."""
import argparse,json,struct
from pathlib import Path
import certify_compact_r17_candidates as cert
import score_full11952_retained as score
import full11952_64_r17_pari_batch as batch
import audit_outer_known29_score as control
from research_runtime.store import checkpoint
OUT=score.ART/'full11952_known29_retention_v1.json'

def expected():
    p=score.protocol();frozen=batch.protocol();selection=cert.read(score.OUT)
    reference=cert.read(control.OUT);c=reference['control']
    ledger=cert.read(batch.CONTROL/'ledger.json')
    bank=cert.read(score.D/'bank.json');extended=cert.read(score.D/'extended/result.json')
    if ledger['status']!='PASS' or len(ledger['rows'])!=11 or len(frozen['rows'])!=64 or reference['status']!='PASS':raise ArithmeticError('completed score pipeline and frozen point roster required')
    if cert.hashed(score.META)!=bank['metadata_sha256'] or cert.hashed(score.EXT)!=extended['extended_vector_sha256']:raise ArithmeticError('score bank binding differs')
    n,d=c['numerator'],c['denominator'];shortkey=(c['short_score_units'],c['short_good_primes'],-d,-n);longkey=(c['combined_selection_units'],c['combined_good'],-d,-n)
    shortahead=longahead=0;found=[]
    with score.EXT.open('rb') as f:
        if f.read(16)!=b'R17EXT01'+struct.pack('<II',p['rows'],p['extension_primes']):raise ArithmeticError('extended frame differs')
        for i,(nn,dd,s,g) in score.iter_bank(p):
            raw=f.read(score.EXTRA.size)
            if len(raw)!=score.EXTRA.size:raise ArithmeticError('truncated score bank')
            e,h=score.EXTRA.unpack(raw)
            shortahead+=(s,g,-dd,-nn)>shortkey;longahead+=(s+e,g+h,-dd,-nn)>longkey
            if (nn,dd)==(n,d):
                if (s,g,e,h)!=(c['short_score_units'],c['short_good_primes'],c['extension_selection_units'],c['extension_good']):raise ArithmeticError('independently scored control differs')
                found.append(i)
        if f.read(1):raise ArithmeticError('trailing extended bytes')
    if len(found)!=1:raise ArithmeticError('one retained known control required')
    i=found[0];selected=[j+1 for j,r in enumerate(selection['selected']) if r['retained_index']==i]
    if len(selected)!=1:raise ArithmeticError('known control must appear once in frozen selection')
    parent=score.scan.protocol();slice_row=parent['rows'][i//512]
    paths=[Path(__file__).resolve(),control.OUT,score.OUT,score.META,score.EXT,score.D/'bank.json',score.D/'short/result.json',score.D/'extended/result.json',batch.D/'protocol.json',batch.D/'maps-ledger.json',batch.CONTROL/'ledger.json']
    return {'schema':'elliptic-curves.full11952-known29-retention.v1','status':'PASS','sources':{str(q.relative_to(score.ROOT)):cert.hashed(q) for q in paths},'parameter':str(cert.F(n,d)),'public_id':12,'retained_rows':p['rows'],'retained_index':i,'slice_id':slice_row['id'],'short_position_within_slice':i%512+1,'short_position_among_retained':shortahead+1,'extended_position_among_retained':longahead+1,'position_in_frozen64':selected[0],'claim_boundary':'Retrospective exact placement of the independently scored known29 control after prospective selection and all point maps were frozen. It was retained and selected without a public label entering score ordering. The curve is already known and cannot be a new inventory addition. These finite retained-population positions establish neither full-population extended-score optimality nor prospective sensitivity, rank density or universal novelty. The prospective roster, point budgets and execution remain unchanged.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('retention audit differs')
    else:
        if OUT.exists():raise FileExistsError('preserve control retention audit')
        checkpoint(OUT,d)
    print('KNOWN29 RETAINED',d['retained_index'],'SHORT POSITION',d['short_position_among_retained'],'EXTENDED POSITION',d['extended_position_among_retained'],'FROZEN64 POSITION',d['position_in_frozen64'],flush=True)
