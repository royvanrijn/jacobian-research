#!/usr/bin/env python3
"""Immutable, retrospective rank-jump panel. No search, CAS or cache mutation.

capture reads pinned git blobs and hash-bound completed worker transcripts.
replay reads only the portable captured input and recomputes the diagnostics.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/generated-results/elliptic-curves'
BASE = '7471645748ae07b8d942f6957e314378ca8aed1f'
INPUT = OUT / 'rank_jump_retrospective_inputs_v1.json'
RESULT = OUT / 'rank_jump_retrospective_report_v1.json'


def digest(b): return hashlib.sha256(b).hexdigest()
def write_new(p, obj):
    with p.open('x') as f: json.dump(obj, f, indent=2, sort_keys=True); f.write('\n')
def read(p): return json.loads(p.read_bytes())
def pack(row): return sum(int(x) << i for i, x in enumerate(row))
def basis(vectors):
    piv = {}
    for v in vectors:
        while v:
            k = v.bit_length()-1
            if k not in piv: piv[k] = v; break
            v ^= piv[k]
    return piv

def reduce(v, piv):
    for k in sorted(piv, reverse=True):
        if v >> k & 1: v ^= piv[k]
    return v

def rank(vs): return len(basis(vs))
def transpose(rows):
    return [pack([r[i] for r in rows]) for i in range(len(rows[0]))]

def short(model, points):
    a1,a2,a3,a4,a6 = map(F,model)
    b2=a1*a1+4*a2; b4=a1*a3+2*a4; b6=a3*a3+4*a6
    A=b4/2-b2*b2/48; B=b6/4-b2*b4/24+b2**3/864
    pts=[]
    for x,y in points:
        x,y=F(x),F(y)
        assert y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6
        pts.append([str(x+b2/12),str(y+(a1*x+a3)/2)])
    return ['0','0','0',str(A),str(B)],pts

def primes(n):
    return [p for p in range(3,n+1,2) if all(p%d for d in range(3,isqrt(p)+1,2))]
def mod(a,p):
    a=F(a); return a.numerator*pow(a.denominator,-1,p)%p

@lru_cache(None)
def roots_at(A,B,p):
    try: a,b=mod(A,p),mod(B,p)
    except ValueError: return None
    if (4*a**3+27*b*b)%p==0: return None
    return tuple(x for x in range(p) if (x**3+a*x+b)%p==0)

def point_signature(model, point, blocks):
    """Exact split-root Kummer characters, including reduction to O and 2-torsion."""
    A=model[3]; x,y=map(F,point); ans=0; col=0
    for p,roots in blocks:
        if x.denominator%p==0:
            # A good integral model: nonintegral rational points reduce to O.
            col+=len(roots); continue
        xp=mod(x,p)
        for r in roots:
            v=(xp-r)%p
            if v==0: v=(3*r*r+mod(A,p))%p
            assert v
            ans |= (pow(v,(p-1)//2,p)==p-1) << col
            col+=1
    return ans

def galois(model):
    A,B=map(F,model[3:]); disc=-4*A**3-27*B**2
    assert disc
    square=disc>0 and isqrt(disc.numerator)**2==disc.numerator and isqrt(disc.denominator)**2==disc.denominator
    witness=next((p for p in primes(1999) if roots_at(str(A),str(B),p)==()),None)
    return {'cubic_discriminant':str(disc),'discriminant_rational_square':square,
            'irreducibility_prime':witness,'galois_group':'C3' if witness and square else 'S3' if witness else 'UNKNOWN',
            'layer':'incidence','claim':'Equation-only E[2] representation test, not a rank predictor.'}


def capture():
    bindings={}; rows={}; charts={}
    def blob(path):
        data=subprocess.check_output(['git','show',BASE+':'+path],cwd=ROOT)
        bindings[path]={'sha256':digest(data),'git_commit':BASE}
        return json.loads(data)
    def art(name): return blob('artifacts/generated-results/elliptic-curves/'+name)
    def worker(w):
        path=w['path']; raw=(ROOT/path).read_bytes()
        if digest(raw)!=w['sha256']: raise ValueError('Live transcript differs from pinned certificate: '+path)
        bindings[path]={'sha256':digest(raw),'role':'completed_transcript_projection'}
        return json.loads(raw)
    def add_recent(r,source,role='prospective_discovery',with_charts=False):
        fam=r.get('family','published-R17'); key=fam+':'+r['parameter']
        if key in rows and not with_charts:return key
        w=r.get('discovery_witness',r.get('input'))
        d=worker(w) if w and with_charts else None
        gr=16 if fam.startswith('a1') else 17
        gp=r.get('generic_points',d.get('generic_points') if d else None)
        if gp is None and d:
            gp=d['initial_state']['state']['reductions']['points']
        if gp is None: gp=r['points'][:gr] # source certificates bind this ordered prefix to transported sections
        assert len(gp)==gr
        points=r['points']; proof=r['rank_certificate']
        rows[key]={'id':key,'family':fam,'parameter':r['parameter'],'origin':role,'source':source,
                   'model':r['curve'],'generic_points':gp,'points':points,'generic_rank':gr,
                   'rank_lower_bound':proof['rank_lower_bound'],'prime_list':[s['prime'] for s in proof['signatures']],
                   'generic_rank_scope':'specialized marked subgroup; full finite Kummer rank required below',
                   'score_units':r.get('score_units'),'cohort':r.get('cohort'),
                   'chart_count':r.get('completed_charts',r.get('chart_count'))}
        if d:
            rows[key]['chart_count']=len(d['charts'])
            rows[key]['cohort']=str(Path(w['path']).parts[3])
            pop=ROOT/Path(w['path']).parent.parent/'population.json'
            if pop.exists():
                raw=pop.read_bytes(); popd=json.loads(raw)
                candidates=popd.get('finalists',popd.get('selected',[]))
                for x in candidates:
                    if x.get('parameter')==r['parameter']:
                        rows[key]['score_units']=x.get('score_units')
                        bindings[str(pop.relative_to(ROOT))]={'sha256':digest(raw),'role':'completed_selection_metadata'}
            if with_charts:
                charts[key]={'source':w,'charts':[{'index':i,'centre':c['centre'],
                  'recorded_rank_lower_bound':c['rank_lower_bound'],
                  'points':[[p['x'],p['y']] for p in c['search']['finite_curve_points']],
                  'completed_denominator':c['search']['completed_denominator'],
                  'denominator_end':c['search']['denominator_end'],
                  'coefficient_bits':c['search']['maximum_coefficient_bits']}
                  for i,c in enumerate(d['charts'])]}
        return key
    index=art('new_high_rank_curve_index_v2.json')
    sources={n:art(n) for n in index['source_certificate_hashes']}
    detailed={'07ca9:505/794','a1-fibration-05:307/206','a1-fibration-04:-1647/91','published-R17:-2300/843'}
    for r in index['curves']:
        sr=sources[r['source_certificate']]['curves'][r['source_curve_index']]
        assert sr['parameter']==r['parameter']
        add_recent(sr,r['source_certificate'],with_charts=r['family']+':'+r['parameter'] in detailed)
    r=art('compact_r17_wide_rank26_initial_proof_v1.json')
    r['discovery_witness']={'path':r['input_path'],'sha256':r['input_sha256']}
    add_recent(r,'compact_r17_wide_rank26_initial_proof_v1.json',with_charts=True)
    for name in ('prospective_mw16_results_v1.json','prospective_mw16_wide_results_v1.json'):
        for r in art(name)['curves']:
            key=r['family']+':'+r['parameter']
            add_recent(r,name,role='retrospective_search_control',with_charts=key in ('a1-fibration-05:-3158/1291','a1-fibration-04:-2177/2397'))
    measures=art('compact_r17_initial_measurements_v1.json')['rows']
    high=next(r for r in measures if r['parameter']=='-2300/843')
    # Predeclared deterministic match: same cohort, zero observed gain, closest coefficient bits then score.
    lows=[r for r in measures if len(r['points'])==17 and r['cohort']==high['cohort']]
    low=min(lows,key=lambda r:(abs(r['short_coefficient_bits']-high['short_coefficient_bits']),abs(r['score_units']-high['score_units']),r['parameter']))
    add_recent(high,'compact_r17_initial_measurements_v1.json',with_charts=True)
    rk=add_recent(low,'compact_r17_initial_measurements_v1.json',role='retrospective_search_control',with_charts=True)
    for r in art('exceptional_soluble_selmer_panel_v1.json')['curves']:
        gp=r['generic_points']; pts=gp+[p['point'] for p in r['candidate_points']]
        model,pts=short(r['model'],pts)
        key='ICARM-'+str(r['curve_id'])
        rows[key]={'id':key,'family':r['frame'],'parameter':None,'origin':'historic_oracle',
          'source':'exceptional_soluble_selmer_panel_v1.json','model':model,
          'generic_points':pts[:len(gp)],'points':pts,'generic_rank':r['generic_rank'],
          'rank_lower_bound':r['generic_rank']+r['known_soluble_residual_dimension_lower_bound'],
          'prime_list':[b['prime'] for b in r['exact_kummer_signature']['blocks']],
          'cover_count':len(r['basis_covers']),'cover_labels':[c['label'] for c in r['basis_covers']],
          'ideal_data':'UNKNOWN: no full class-group certificate','full_selmer':'UNKNOWN'}
    for r in art('elkies_2026_known_kummer_quotients_suite_v1.json')['runs']:
        key=r['case_id']; model,_=short(r['global_minimal_model'],[])
        rows[key]={'id':key,'family':'published-R17','parameter':r['parameter'],
          'origin':'historic_oracle' if r['exceptional_point_count'] else 'retrospective_search_control',
          'source':'elkies_2026_known_kummer_quotients_suite_v1.json','model':model,'generic_rank':17,
          'rank_lower_bound':r['known_kummer_rank'],
          'supplied_global_kummer_fingerprints':[p['local_squareclass_row'] for p in r['points']],
          'signature_scope':'Exact row algebra replay of prior Kummer certificate; points not re-evaluated.'}
    deep=art('elkies_2026_deep_cover_exceptional_quotients_v1.json')
    covers=[{k:r[k] for k in ('parameter','known_rank_lower_bound','cumulative_captured_exceptional_rank','tested_cover_universe')} for r in deep['fibres']]
    pairs=[{'high':'a1-fibration-05:307/206','low':'a1-fibration-05:-3158/1291','match':'same family and 43-chart policy, adjacent score cohorts; scale mismatch retained'},
           {'high':'a1-fibration-04:-1647/91','low':'a1-fibration-04:-2177/2397','match':'same family, same H4096 cohort, same 43-chart policy; coefficient size not matched'},
           {'high':'published-R17:-2300/843','low':rk,'match':'same cohort; nearest coefficient bits among observed-zero controls, then nearest score'}]
    write_new(INPUT,{'schema':'rank-jump.retrospective-inputs.v1','baseline':BASE,'bindings':bindings,
       'rows':list(rows.values()),'chart_inputs':charts,'pairs':pairs,'low_degree_cover_capture':covers,
       'limits':{'point_searches':0,'parameter_sweeps':0,'full_descents':0,'prime_bound':1999,
                 'chart_scope':'only retained points from eight completed initial transcripts'},
       'oracle_policy':'All diagnostics retrospective. No point or outcome enters a prospective selector.'})


def characterize(row):
    model=row['model']; g=row['generic_rank']
    gal=galois(model)
    assert gal['irreducibility_prime'] is not None, 'Independence proof requires E(Q)[2]=0'
    if 'supplied_global_kummer_fingerprints' in row:
        sig=[pack(r) for r in row['supplied_global_kummer_fingerprints']]; blocks=[]
    else:
        _,pts=short(model,row['generic_points']+row['points'])
        blocks=[(p,roots_at(model[3],model[4],p)) for p in row['prime_list']]
        blocks=[(p,rs) for p,rs in blocks if rs]
        sig=[point_signature(model,p,blocks) for p in pts]
    gr=rank(sig[:g]); total=rank(sig); assert gr==g,(row['id'],gr,g)
    assert total>=row['rank_lower_bound'],(row['id'],total,row['rank_lower_bound'])
    selected=[];span=[]
    for i,v in enumerate(sig):
        if rank(span+[v])>len(span):selected.append(i);span.append(v)
    q=total-g
    out={k:row.get(k) for k in ('id','family','parameter','origin','source','score_units','cohort','chart_count')}
    out.update({'generic_subgroup_rank_exact':g,'generic_kummer_dimension_exact':gr,'epsilon_M':0,
      'certified_independent_subgroup_rank_exact':total,'curve_rank_lower_bound':total,'curve_rank_upper_bound':None,
      'independent_input_indices':selected,
      'independent_input_order':'generic points followed by supplied points, or the historic signature-only source order',
      'certified_independent_quotient_rank_exact':q,'full_quotient_rank_interval':[q,None],
      'soluble_residual_dimension_lower_bound':q,'galois':gal,
      'short_coefficient_bits':max(max(F(x).numerator.bit_length(),F(x).denominator.bit_length()) for x in model[3:]),
      'parameter_height':max(abs(F(row['parameter']).numerator),F(row['parameter']).denominator) if row['parameter'] else None,
      'block_report':{'witnessed_soluble_block_dimension':q,
        'ct_rows_against_every_selmer_class':'ZERO_BY_RATIONAL_WITNESSES' if q else 'NO_RESIDUAL_WITNESS',
        'ct_complement':'UNKNOWN','ideal_class_blocks':'UNKNOWN','common_auxiliary_curve':'UNKNOWN',
        'common_descent_algebra':'the same irreducible 2-division cubic for every direction on this fibre',
        'causal_decomposition':'UNKNOWN; an independent basis alone defines no preferred arithmetic blocks'}})
    return out,blocks,sig


def chart_map(row,blocks,chart_input):
    g=row['generic_rank']; model=row['model']
    gp=[point_signature(model,p,blocks) for p in row['generic_points']]
    piv={}; chosen=[]
    def insert(v):
        rem,word=solve(v)
        if rem:
            j=len(chosen);chosen.append(v);piv[rem.bit_length()-1]=(rem,word^(1<<j))
            return True
        return False
    def solve(v):
        word=0
        for k in sorted(piv,reverse=True):
            if v>>k&1:v^=piv[k][0];word^=piv[k][1]
        return v,word
    for v in gp:assert insert(v)
    for p in row['points']:insert(point_signature(model,p,blocks))
    q=len(chosen)-g; allm=[]; records=[]; exact_points={}
    for c in chart_input['charts']:
        ms=[]; unresolved=0
        _,points=short(model,c['points'])
        for p in points:
            v=point_signature(model,p,blocks); rem,w=solve(v)
            if rem:unresolved+=1;continue
            m=w>>g;ms.append(m)
            # Keep one actual representative per finite quotient signature, up to sign.
            if m:exact_points.setdefault(m,p)
        vals=sorted(set(ms)-{0});before=rank(allm);allm+=vals
        records.append({'index':c['index'],'generic_parity_mask':c['centre']['mask'],
          'centre_metric_norm':c['centre']['metric_norm'],'finite_quotient_masks':vals,
          'returned_point_count':len(points),'distinct_nonzero_signatures':len(vals),
          'standalone_quotient_rank':rank(vals),'prefix_quotient_rank':rank(allm),
          'new_prefix_rank':rank(allm)-before,'outside_reference_finite_span':unresolved,
          'recorded_rank_lower_bound':c['recorded_rank_lower_bound'],
          'denominator_coverage':[c['completed_denominator'],c['denominator_end']],
          'coefficient_bits':c['coefficient_bits']})
    # Deterministic greedy chart span cover, not a claim of minimum chart count.
    union=[];chosen_charts=[];remaining=records[:]
    while rank(union)<rank(allm):
        c=max(remaining,key=lambda c:(rank(union+c['finite_quotient_masks'])-rank(union),-c['index']))
        chosen_charts.append(c['index']);union+=c['finite_quotient_masks'];remaining.remove(c)
    overlaps=[{'left':a['index'],'right':b['index'],
      'intersection_dimension':a['standalone_quotient_rank']+b['standalone_quotient_rank']-rank(a['finite_quotient_masks']+b['finite_quotient_masks']),
      'parity_hamming_distance':(a['generic_parity_mask']^b['generic_parity_mask']).bit_count()}
      for i,a in enumerate(records) for b in records[:i]]
    return {'id':row['id'],'layer':'visibility','reference_quotient_dimension':q,
      'signature_method':'finite split-root characters modulo the generic image; injective on the certified basis span only',
      'limitation':'A matching finite signature is not a certified global relation for an arbitrary returned point. Ranks of images give lower bounds; counts of masks need not count global classes.',
      'charts':records,'chart_pair_intersections':overlaps,'greedy_span_cover_chart_indices':chosen_charts,
      'prefix_charts_for_full_observed_span':next((c['index']+1 for c in records if c['prefix_quotient_rank']==q),None) if q else 0,
      'maximum_single_chart_rank':max(c['standalone_quotient_rank'] for c in records),
      'distinct_nonzero_signatures':len(set(allm)),
      'finite_mask_representatives':{str(k):v for k,v in exact_points.items()},
      'signature_union_rank':rank(allm)}


def replay(check=False):
    inp=read(INPUT); panel=[];maps=[]
    for row in inp['rows']:
        out,blocks,sigs=characterize(row);panel.append(out)
        if row['id'] in inp['chart_inputs']:
            maps.append(chart_map(row,blocks,inp['chart_inputs'][row['id']]))
        print(row['id'],'subgroup',out['generic_subgroup_rank_exact'],'+',out['certified_independent_quotient_rank_exact'],out['galois']['galois_group'],flush=True)
    result={'schema':'rank-jump.retrospective-report.v1','input_sha256':digest(INPUT.read_bytes()),
      'script_sha256':digest(Path(__file__).read_bytes()),'baseline':inp['baseline'],
      'panel':panel,'chart_maps':maps,'pairs':inp['pairs'],'low_degree_cover_capture':inp['low_degree_cover_capture'],
      'summary':{'rows':len(panel),'galois_counts':dict(Counter(r['galois']['galois_group'] for r in panel)),
        'all_generic_images_certified':True,'all_zero_gains_censored':True},
      'claim_boundary':'No new curve/rank discovery, rank upper bound, full descent or prospective selection. Replays exact membership, independent split-root signatures and binary incidence; prior signature-only historic controls replay row algebra.'}
    if check:
        if read(RESULT)!=result:raise ValueError('report mismatch')
        print('PASS deterministic replay')
    else:write_new(RESULT,result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('capture','build','check'));a=p.parse_args()
    capture() if a.mode=='capture' else replay(a.mode=='check')
