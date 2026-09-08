#!/usr/bin/env python3
"""Retrospective exact visibility of public curve542 points in frozen MW16 charts."""
import argparse
from hashlib import sha256
from collections import Counter
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import search_observability as visibility
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1/a1-fibration-04/candidate-00/result.json'

def transport(ainvs,points,model):
    original=tuple(map(cert.F,ainvs));a1,a2,a3,a4,a6=original;c2=a2+a1*a1/4
    a,b=cert.weierstrass_invariants(original),cert.weierstrass_invariants(model)
    u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    if u is None or not cert.isomorphic(original,model):raise ArithmeticError('public equation transport failed')
    out=[]
    for x,y in points:
        x,y=cert.F(x),cert.F(y);p=((x+c2/3)/u**2,(y+(a1*x+a3)/2)/u**3)
        if not cert.is_on_weierstrass_curve(model,p):raise ArithmeticError('public point transport failed')
        out.append(p)
    return str(u),out

def summaries(charts,points):
    counts=Counter();minima=[];discrepancies=[]
    for j,(x,y) in enumerate(points):
        choices=[]
        for sign in (1,-1):
            for i,chart in enumerate(charts):
                v=visibility.point_visibility(chart,(x,sign*y));counts[v['status']]+=1
                row={'public_point_index':j,'sign':sign,'chart_index':i,**v}
                if v['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(row)
                if v.get('minimum_affine_height') is not None:choices.append(row)
        minima.append(min(choices,key=lambda r:(r['minimum_affine_height'],r['chart_index'],r['sign'])) if choices else {'public_point_index':j,'no_finite_parameter':True})
    return {'status_counts':dict(counts),'minimum_by_public_point':minima,'coverage_discrepancies':discrepancies}

def sources():
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),Path(visibility.__file__).resolve(),
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def build(output):
    if output.exists():raise FileExistsError('preserve visibility audit')
    data=cert.read(INPUT)
    if data['status']!='COMPLETE_DECLARED_PILOT' or data['rank_lower_bound']!=25 or len(data['charts'])!=43:raise ArithmeticError('frozen initial measurement changed')
    public=next(r for r in cert.read(cert.DATABASE)['curves'] if r['id']==542)
    model=tuple(map(cert.F,data['curve']));u,points=transport(public['ainvs'],public['points'],model)
    basis=[tuple(map(cert.F,p)) for p in data['final_state']['state']['reductions']['points']]
    witness=None
    for j,p in enumerate(points):
        try:proof=cert.checked_rank(model,basis+[p])
        except ArithmeticError:continue
        witness={'public_point_index':j,'points':[list(map(str,q)) for q in basis+[p]],'rank_certificate':proof};break
    if witness is None:raise ArithmeticError('rank26 witness not certified within finite prime bound')
    fields=('short_model','short_model_x_shift','base_point','pointed_chart','horizontal_matrix','ordinate_scale','coefficients',
            'height_bound','denominator_start','denominator_end','completed_denominator','finite_curve_points','infinity_checked','primitive_square_hits')
    charts=[{'input':{'curve':r['search']['input']['curve']},**{k:r['search'][k] for k in fields}} for r in data['charts']]
    result={'schema':'elliptic-curves.curve542-initial-visibility.v1','sources':sources(),'retrospective_only':True,
        'input_path':str(INPUT.relative_to(ROOT)),'input_sha256':cert.hashed(INPUT),'catalogue_sha256':cert.hashed(cert.DATABASE),
        'public_model':public['ainvs'],'public_points':public['points'],'curve':data['curve'],'scale_u':u,
        'transported_points':[list(map(str,p)) for p in points],'rank26_witness':witness,'charts':charts,
        'full_chart_hashes':[sha256(json.dumps(r['search'],sort_keys=True,separators=(',',':')).encode()).hexdigest() for r in data['charts']],
        'claim_boundary':'Public oracle points are used only after the fixed initial batch and its exact replay. No selector or search is rerun or changed. Signed representative visibility is not visibility of every translate or a rank upper bound.'}
    result.update(summaries(charts,points));checkpoint(output,result)
    print('CURVE542 RANK26 WITNESS',witness['public_point_index'],'BEST VISIBILITY',result['minimum_by_public_point'][witness['public_point_index']],flush=True)
    print('VISIBILITY COUNTS',result['status_counts'],'discrepancies',len(result['coverage_discrepancies']),flush=True)

def check(path):
    d=cert.read(path)
    if d['sources']!=sources():raise ArithmeticError('audit source changed')
    model=tuple(map(cert.F,d['curve']));u,points=transport(d['public_model'],d['public_points'],model)
    if u!=d['scale_u'] or [list(map(str,p)) for p in points]!=d['transported_points']:raise ArithmeticError('public transport changed')
    w=d['rank26_witness'];q=[tuple(map(cert.F,p)) for p in w['points']];old=w['rank_certificate']
    if len(q)!=26 or q[-1]!=points[w['public_point_index']]:raise ArithmeticError('escape witness changed')
    actual=cert.checked_rank(model,q,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
    if json.dumps(actual,sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('rank26 certificate changed')
    result=summaries(d['charts'],points)
    if any(result[k]!=d[k] for k in result):raise ArithmeticError('visibility audit changed')
    print('REPLAYED CURVE542 RANK26 AND INITIAL VISIBILITY',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.output)
