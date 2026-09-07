#!/usr/bin/env python3
"""Independently replay the cascade evidence and publish a compact report."""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v1'
OUT=ART/'adaptive_visibility_cascade_v1.json'
from research_runtime.store import checkpoint
from half_lattice_pointed_sieve import linear_combination
from alternate_quartic_covers import short_add,alternate_cover
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend
import audit_recorded_point_mod2_rank_v3 as mod2
import audit_retained_cloud_modl as modl
from search_observability import transform
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))

def matrices():
    vis=read(ART/'curve302_residual_visibility_geometry_v1.json')
    cloud=read(ART/'curve302_focused_curve302_generic17_mod2_v1.json');model=tuple(map(F,cloud['curve']))
    base=tuple(tuple(map(F,p)) for p in cloud['independent_points'][:17]);targets={}
    for e in vis['directions']:
        roots={short_add(model,tuple(map(F,c['exact_target_point_short_model'])),linear_combination(model,base,[-int(x) for x in c['target_translation_m17_word']])) for c in e['exact_pointed_quartic_rows']}
        if len(roots)!=1:raise ArithmeticError('target labels changed between original charts')
        targets[e['id']]=roots.pop()
    reports=[]
    for dirname,filename in [('visibility-cascade-matrix-v1','curve302_visibility_cascade_matrix_v1.json'),('visibility-cascade-support-matrix-v1','curve302_visibility_cascade_support_matrix_v1.json')]:
        directory=D.parent/dirname;data=read(ART/filename);protocol=read(directory/'protocol.json')
        if sha(directory/'protocol.json')!=data['protocol_sha256']:raise ArithmeticError('matrix protocol changed')
        for path,h in data['inputs'].items():
            if sha(ROOT/path)!=h:raise ArithmeticError('matrix input changed '+path)
        subgroup={}
        for path in data['inputs']:
            if path.endswith('_mod2_v1.json'):
                c=read(ROOT/path);subgroup[c['rank_lower_bound']]=tuple(tuple(map(F,p)) for p in c['independent_points'])
        subgroup[17]=base
        checked=0;last={};counts={}
        for c in data['cells']:
            name=c['direction'];rank=c['rank'];basis=subgroup[rank];target=targets[name]
            full=read(directory/f'{name}-M{rank}.json')
            if {k:v for k,v in full.items() if k!='trials'}!=c:raise ArithmeticError('matrix summary differs from trials')
            for r in full['trials']:
                centre=linear_combination(model,basis,r['centre_word']);moved=short_add(model,target,linear_combination(model,basis,r['translation_word']))
                if r['status']=='EXACT_SUBGROUP_CONTAINMENT_ENDPOINT':
                    if moved is not None and moved!=centre:raise ArithmeticError('false containment endpoint')
                else:
                    if list(map(str,centre))!=r['centre_point'] or list(map(str,moved))!=r['translated_target']:raise ArithmeticError('exact subgroup translation differs')
                    n,d=map(int,r['coordinate']);h=max(abs(n),abs(d))
                    if gcd(n,d)!=1 or str(h)!=r['height']:raise ArithmeticError('nonprimitive coordinate or wrong height')
                    a,b,cc,dd=map(F,r['coordinate_matrix']);t,w=alternate_cover(model,centre).curve_point_to_cover(moved)
                    if t*(cc*n+dd*d)!=a*n+b*d:raise ArithmeticError('horizontal coordinate failed')
                    coef=tuple(map(F,r['quartic']));root=F(r['square_root'])
                    if sum(v*n**i*d**(4-i) for i,v in enumerate(coef))!=root*root:raise ArithmeticError('quartic square failed')
                    x,y=centre;raw=(-3*x*x-4*model[3],-8*y,-6*x,F(0),F(1));trans=transform(raw,(a,b,cc,dd))
                    ratio=next(v/q for v,q in zip(trans,coef) if q)
                    if trans!=tuple(ratio*q for q in coef) or ratio<=0 or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator:raise ArithmeticError('quartic twist or transport changed')
                checked+=1
            if c['winner']:
                h=int(c['winner']['height'])
                if name in last and h>last[name]:raise ArithmeticError('persistent minimum increased')
                last[name]=h
            if c['cheap_representatives']<counts.get(name,0):raise ArithmeticError('cheap witness count decreased')
            counts[name]=c['cheap_representatives']
        reports.append({'artifact':rel(ART/filename),'sha256':sha(ART/filename),'cells':len(data['cells']),'exact_translations_and_quartic_witnesses_checked':checked,'status':'PASS'})
    return reports

def build():
    p=read(D/'protocol.json');cp=read(D/'continuation-v1/protocol.json');roster=read(D/'roster.json')
    for protocol in (p,cp):
        for key in ('sources','inputs'):
            for path,h in protocol.get(key,{}).items():
                if sha(ROOT/path)!=h:raise ArithmeticError('frozen input/source changed '+path)
    rows=[];files={rel(D/'protocol.json'):sha(D/'protocol.json'),rel(D/'continuation-v1/protocol.json'):sha(D/'continuation-v1/protocol.json')}
    for row in roster:
        d=D/row['id'];terminal=read(D/'continuation-v1'/row['id']/'terminal.json')
        if terminal['prefix_sha256']!=sha(d/'terminal.json'):raise ArithmeticError('continuation prefix changed')
        seed=read(d/'seed.json');model=tuple(map(F,seed['curve']));priorbasis=seed['points'];bits=[];seconds=[];allhits=0;chartcount=0
        for s in terminal['stages']:
            cloudpath=ROOT/s['mod2'];oddpath=ROOT/s['modl'];wd=cloudpath.parent;cloud=read(cloudpath);odd=read(oddpath)
            if sha(cloudpath)!=s['mod2_sha256'] or sha(oddpath)!=s['modl_sha256']:raise ArithmeticError('rank certificate hash changed')
            if cloud['input_sha256']!=sha(wd/'result.json') or odd['input_sha256']!=sha(cloudpath):raise ArithmeticError('certificate input changed')
            mod2.check(cloudpath);modl.check(oddpath)
            select=read(wd/'selection.json');basis=tuple(tuple(map(F,pt)) for pt in select['basis'])
            if select['basis']!=priorbasis or select['rank']!=s['before']:raise ArithmeticError('adaptive seed link failed')
            result=read(wd/'result.json')
            if cloud['independent_points'][:len(priorbasis)]!=priorbasis:raise ArithmeticError('subgroup inclusion lost')
            for chart,c in zip(result['charts'],select['centres']):
                if chart['centre']!=c:raise ArithmeticError('chart differs from frozen selection')
                point=linear_combination(model,basis,c['representative']);record=chart['search'];mapping=chart['mapping']
                search=PointedQuarticSearch(curve=model,subgroup=(),centre={'point':point},coordinate_policy=mapping['coordinate_policy'])
                backend.validate_map(search,mapping)
                if tuple(map(F,record['coefficients']))!=search.coefficients or record['source_hashes']!=backend.sources():raise ArithmeticError('backend geometry/source changed')
                hits,points,ms=backend.witnesses(search,record['stdout'],record['status'],p['height'])
                if hits!=[tuple(map(int,h)) for h in record['primitive_square_hits']] or [dict(zip(('x','y'),map(str,pt))) for pt in points]!=record['finite_curve_points']:raise ArithmeticError('independent witness/map replay failed')
                if record['program']!=backend.program(mapping,p['height']) or record['height_bound']!=p['height'] or record['timeout_seconds']!=p['seconds_per_chart']:raise ArithmeticError('chart budget changed')
                allhits+=len(hits);chartcount+=1;bits.append(record['maximum_coefficient_bits']);seconds.append(record['wall_seconds'])
            if len(result['charts'])!=len(select['centres']) or len(result['charts'])!=s['charts']:raise ArithmeticError('chart coverage mismatch')
            files.update({rel(x):sha(x) for x in wd.glob('*.json')});priorbasis=cloud['independent_points']
        last=read(ROOT/terminal['stages'][-1]['mod2']);odd=read(ROOT/terminal['stages'][-1]['modl'])
        rows.append({'id':row['id'],'stratum':row.get('stratum','calibration'),'height_bin':row.get('height_bin'),
            'rank_path':[17]+[s['after'] for s in terminal['stages']],'stages':terminal['stages'],'charts':chartcount,
            'completed_charts':sum(s['completed'] for s in terminal['stages']),'exact_square_hits':allhits,
            'quartic_bits_min_median_max':[min(bits),sorted(bits)[len(bits)//2],max(bits)],'search_seconds':sum(seconds),
            'final_mod2_certificate':last,'final_modl_certificate':odd,
            'stop':'complete_no_gain_wave' if terminal['stages'][-1]['before']==terminal['stages'][-1]['after'] else 'budget_or_rank_endpoint'})
        files[rel(d/'terminal.json')]=sha(d/'terminal.json');files[rel(D/'continuation-v1'/row['id']/'terminal.json')]=sha(D/'continuation-v1'/row['id']/'terminal.json')
    matrix_checks=matrices()
    return {'schema':'adaptive-visibility-cascade-report.v1','status':'PASS_BOUNDED_AUTONOMOUS_CASCADE_AND_FINITE_ATLAS_REPLAY',
        'source_sha256':sha(Path(__file__)),'initial_protocol':p,'budget_continuation_protocol':cp,'rows':rows,'matrix_checks':matrix_checks,'checkpoints':files,
        'boundary':'The selector was calibrated on known curve302 history, then executed without target points or residual labels. One explicit budget extension followed four gaining waves; selection never changed. Controls are ten score-stratified siblings, not an unbiased population sample. Matrix minima are finite-atlas bounds. No exact-rank upper bound or rank32 discovery is inferred.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();result=build()
    if a.check:
        if read(OUT)!=result:raise ArithmeticError('report differs')
    else:
        if OUT.exists():raise FileExistsError('preserve published report')
        checkpoint(OUT,result)
    print('CASCADE REPORT PASS',[(r['id'],r['rank_path']) for r in result['rows']],flush=True)
