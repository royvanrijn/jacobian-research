#!/usr/bin/env python3
"""Retrospective one-control score audit after the fixed outer48 trial closes."""
import argparse,json
from math import log
from pathlib import Path
import certify_compact_r17_candidates as cert
import extend_outer131072_r17 as extension
import benchmark_periodic_nagao_scanner as short
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=extension.ROOT;ART=extension.ART;CAS=extension.CAS;D=extension.parent.LOCAL/'outer-known29-score-audit-v1'
OUT=ART/'outer_known29_score_audit_v1.json';HEIGHTS=ART/'public_compact_parameter_heights_v1.json';TRIAL=ART/'outer48_experiment_v1.json'

def sources():
    return {**extension.engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),HEIGHTS,TRIAL,extension.D/'result.json',CAS/'extend_retained_r17_prime_scores.py',extension.spec.ATLAS]}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve retrospective score protocol')
    trial=cert.read(TRIAL);height=cert.read(HEIGHTS);control=next(r for r in height['rows'] if r['id']==12)
    if trial['status']!='PASS' or trial['completed_point_boxes']!=2160 or trial['rank_lower_bound_counts']!={'17':48}:raise ArithmeticError('completed fixed outer48 trial required')
    if control['status']!='PASS_EXACT_COMPACT_Q_ISOMORPHISM' or control['family']!='11952' or control['compact_parameter']!='89074/31895':raise ArithmeticError('pinned native control differs')
    f=next(r for r in cert.read(extension.spec.ATLAS)['families'] if r['family']=='11952');model=extension.model_at(f,control['compact_parameter'])
    if model!=control['compact_specialization'] or not cert.isomorphic(model,control['public_curve']):raise ArithmeticError('control model transport differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.outer-known29-score-audit.v1','sources':sources(),'family':'11952','parameter':control['compact_parameter'],'model':model,'public_id':12,'gp_sha256':cert.hashed(extension.engine.GP),'wall_seconds':120,'gp_call_wall_seconds':20,'rss_bytes':1073741824,'direct_validation_indices':[0,1,len(extension.engine.PRIMES)//2,len(extension.engine.PRIMES)-1], 'gate':'The completed48-curve outer trial returns no direction beyond the generic17 subgroup, while the known native11952 rank29 control in this envelope has a successful generic17-only detector. Retrospectively locate that one control under the unchanged short and extended scores among the saved outer candidates to distinguish a possible selection miss from an undiagnosed point-exposure or incidence issue.','scope':'Exactly one known curve, one GP trace call through65521 and four direct character-sum checks. Reuse the canonical short-prime tables. This is retrospective labelled calibration after all selection, maps and points are immutable; no new curve, candidate selector, point search, validation-based tie breaking, density or sensitivity theorem follows.'})

def expected(create=False):
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(extension.engine.GP):raise ArithmeticError('frozen score audit inputs changed')
    n,d=cert.F(p['parameter']).as_integer_ratio();units=0;good=0;hashes={}
    for prime in range(5,4094):
        if not short._is_prime(prime):continue
        path=short.old.PARENT/p['family']/'trace-tables'/f'{prime}.json';t=cert.read(path);hashes[str(path.relative_to(ROOT))]=cert.hashed(path);x=n*pow(d,-1,prime)%prime if d%prime else prime
        if t['good'][x]:a=t['traces'][x];units+=round((2-a)/(prime+1-a)*log(prime)*10**12);good+=1
    short.exact_scores([{'numerator':n,'denominator':d,'score_units':units,'good_primes':good}],{'family':p['family']})
    rawpath=D/'raw.json';program=extension.engine.program(p['model'])
    if create:
        if rawpath.exists():raise FileExistsError('preserve one control trace call')
        q=capture([str(extension.engine.GP),'-q','-s','256000000'],input_text=program,limits=Limits(p['gp_call_wall_seconds'],p['rss_bytes']),log_path=D/'gp.log',separate_stderr=True,check=False)
        checkpoint(rawpath,{'program':program,'stdout':q.stdout,'stderr':q.stderr,'supervision':q.supervision})
    raw=cert.read(rawpath)
    if raw['program']!=program or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('control trace failed or censored')
    values,ms=extension.engine.parse(raw['stdout'],p['model']);checks=[]
    for i in p['direct_validation_indices']:
        prime,value=values[i];direct=extension.engine.direct(p['model'],prime)
        if value!=direct:raise ArithmeticError('control direct character sum differs')
        checks.append([prime,value])
    scored={**extension.scoring.sums(values),'combined_selection_units':units+extension.scoring.sums(values)['extension_selection_units'],'combined_good':good+extension.scoring.sums(values)['extension_good'],'numerator':n,'denominator':d}
    saved=cert.read(extension.D/'result.json');rows=[r for r in saved['rows'] if r['family']==p['family']]
    key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator'])
    rank=1+sum(key(r)<key(scored) for r in rows);chosen=[r for r in rows if r['retained_index'] in saved['selection'][p['family']]]
    return {'schema':'elliptic-curves.outer-known29-score-audit-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'raw_sha256':cert.hashed(rawpath),'short_table_hashes':hashes,'control':{**scored,'parameter':p['parameter'],'short_score_units':units,'short_good_primes':good},'direct_checks':checks,'outer_family_candidates':len(rows),'hypothetical_extended_score_position':rank,'selected_family_count':len(chosen),'outranks_selected_count':sum(key(scored)<key(r) for r in chosen),'best_outer_selection_units':max(r['combined_selection_units'] for r in rows),'worst_selected_selection_units':min(r['combined_selection_units'] for r in chosen),'claim_boundary':p['scope']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','build','check']);v=a.parse_args()
    if v.stage=='prepare':prepare()
    else:
        if v.stage=='build' and OUT.exists():raise FileExistsError('preserve one control score report')
        r=expected(v.stage=='build')
        if v.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(r)):raise ArithmeticError('retrospective control score differs')
        else:checkpoint(OUT,r)
        print('KNOWN29 UNCHANGED-SCORE POSITION',r['hypothetical_extended_score_position'],'among',r['outer_family_candidates'],'outer candidates; outranks',r['outranks_selected_count'],'of8 selected',flush=True)
