#!/usr/bin/env python3
"""Replay the fresh4096 equation roster and deterministic64 selection."""
import json
import full11952_late_band_selection as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import digest

def main():
    p=batch.protocol();d=cert.read(batch.D/'result.json');selected=cert.read(batch.OUT)
    if d['status']!='COMPLETE_FROZEN4096' or len(d['rows'])!=4096 or d['protocol_sha256']!=cert.hashed(batch.D/'protocol.json'):raise ArithmeticError('complete fixed scores required')
    projected=sum(r['wall_seconds'] for r in d['rows'][:8])*4096/8
    if projected!=d['projected_serial_seconds'] or projected>p['cost_gate_maximum_projected_serial_seconds']:raise ArithmeticError('runtime gate replay differs')
    old=cert.read(batch.OLD);previous=old['previous_equations']+[{'address':batch.OLD.name+':'+r['id'],'curve':r['curve']} for r in old['curves']]
    catalogue=[{'id':r['id'],'curve':r['ainvs']} for r in cert.read(batch.CAT)['curves']]
    if previous!=p['previous_equations'] or catalogue!=p['catalogue_equations']:raise ArithmeticError('equation exclusion roster differs')
    seen={};rows=[];skips=[]
    def jvalue(model):
        v=cert.weierstrass_invariants(model);return v['c4']**3/v['discriminant']
    for r in previous+catalogue:
        model=tuple(map(cert.F,r['curve']));seen.setdefault(jvalue(model),[]).append((model,str(r.get('address',r.get('id')))))
    prefix=batch.parent.ranked_prefix(batch.parent.protocol())
    if len(prefix)!=32768 or digest(prefix)!=p['old_prefix_sha256']:raise ArithmeticError('fixed prefix differs')
    family=next(f for f in cert.read(batch.parent.extension.spec.ATLAS)['families'] if f['family']=='11952')
    for position,(u,g,nd,nn,i,s,sg,e,eg) in enumerate(prefix):
        n,den=-nn,-nd;t=str(cert.F(n,den));model=tuple(map(cert.F,batch.parent.extension.model_at(family,t)));j=jvalue(model)
        matches=[name for m,name in seen.get(j,[]) if cert.isomorphic(m,model)]
        if matches:skips.append({'retained_index':i,'parameter':t,'matches':matches});continue
        row={'id':f'11952-{i:07}','family':'11952','retained_index':i,'parameter':t,'numerator':n,'denominator':den,
            'score_units':s,'good_primes':sg,'extension_selection_units':e,'extension_good':eg,
            'combined_selection_units':u,'combined_good':g,'model':list(map(str,model)),'old_prefix_position':position+1}
        rows.append(row);seen.setdefault(j,[]).append((model,row['id']))
        if len(rows)==4096:break
    if rows!=p['rows'] or skips!=p['skipped']:raise ArithmeticError('complete fresh equation roster differs')
    combined=[]
    for row,score in zip(rows,d['rows']):
        if row['id']!=score['id'] or score['combined_late_units']!=row['combined_selection_units']+score['scores']['validation_units'] or score['combined_late_good']!=row['combined_good']+score['scores']['validation_good']:raise ArithmeticError('combined scores differ')
        combined.append({**row,**score})
    combined.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
    expected={'schema':'elliptic-curves.full11952-late-band-finalists.v1','status':'PASS_FROZEN64_SELECTION','sources':batch.sources(),
        'protocol_sha256':cert.hashed(batch.D/'protocol.json'),'scores_sha256':cert.hashed(batch.D/'result.json'),
        'score_replay_sha256':cert.hashed(batch.D/'check.supervisor.json'),'selected':combined[:64],'claim_boundary':p['scope']}
    if selected!=expected:raise ArithmeticError('fixed64 selection differs')
    print('FRESH4096 EXCLUSIONS AND FINAL64 ORDER REPLAY PASS',flush=True)
if __name__=='__main__':main()
