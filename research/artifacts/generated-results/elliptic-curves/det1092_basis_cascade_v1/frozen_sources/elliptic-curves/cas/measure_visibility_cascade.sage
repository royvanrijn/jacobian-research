#!/usr/bin/env sage-python
"""Exact retrospective finite-atlas cascade matrix, with persistent witnesses.

The atlas is the original 32 vetted charts per target plus centres of completed
historical recovery arms once available. It is NOT every parity of M/2M.
Old exact witnesses persist; each stage additionally tries a CVP translate.
"""
import argparse,json,hashlib,time
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from sage.all import matrix,ZZ,QQ,pari,vector
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/visibility-cascade-matrix-v1'
VIS=ART/'curve302_residual_visibility_geometry_v1.json'
CHAIN=ART/'curve302_residual_strict_bootstrap_chain_v1.json'
def load(n):return SourceFileLoader('matrix_'+n.replace('.','_'),str(CAS/n)).load_module()
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
from research_runtime.store import checkpoint
from half_lattice_pointed_sieve import linear_combination
from alternate_quartic_covers import short_add

def run():
    if D.exists():raise FileExistsError('preserve matrix experiment')
    D.mkdir(parents=True)
    v=read(VIS);chain=read(CHAIN)
    stages=[(17,ART/'curve302_focused_curve302_generic17_mod2_v1.json'),(19,ART/'curve302_focused_curve302_generic17_mod2_v1.json'),(22,ART/'curve302_recovered_followup_wave_01_mod2_v1.json'),(24,ART/'curve302_recovered_followup_wave_03_mod2_v1.json')]
    arms=chain['recovery_arms']
    stages += [(a['rank_after'],ROOT/next(k for k in a['immutable_checkpoints'] if k.endswith('_mod2_v1.json'))) for a in arms]
    inputs={rel(p):sha(p) for p in {VIS,CHAIN,Path(__file__),*(p for _,p in stages)}}
    checkpoint(D/'protocol.json',{'inputs':inputs,'ranks':[r for r,_ in stages],
        'atlas':'For each of 14 directions retain its 32 original vetted M17 charts; adjoin all seven historical recovery centres when their seed rank is available. Retain earlier witnesses and try one DD/MPFR-agreeing CVP translate per chart/stage. Every winner is a finite-atlas bound, not a global minimum.',
        'cheap_height':125000,'metric':'rounded 384-bit canonical height matrix; no exact canonical CVP optimality claim',
        'endpoint_rule':'O or centre reached by an exact subgroup translation establishes containment; record separately, do not count as new visibility.'})
    cvp=load('audit_curve302_residual_visibility_geometry.sage');helpers=load('audit_curve302_m25_remaining_strict_cvp.sage');geo=load('prospective_half_lattice_v3.sage');mapper=load('factor_free_pari_mapping.sage')
    model=tuple(map(F,read(stages[0][1])['curve']));base17=tuple(tuple(map(F,p)) for p in read(stages[0][1])['independent_points'][:17])
    targets={};atlases={};kept={};cells=[];mapping_cache={}
    for entry in v['directions']:
        name=entry['id'];first=entry['exact_pointed_quartic_rows'][0]
        target=tuple(map(F,first['exact_target_point_short_model']))
        back=linear_combination(model,base17,[-int(x) for x in first['target_translation_m17_word']])
        targets[name]=short_add(model,target,back)
        atlases[name]=[{'word':c['centre_m17_word'],'introduced':17,'orbit':c['degree_two_translation_orbit'],'origin':'original vetted M17 chart'} for c in entry['exact_pointed_quartic_rows']]
        kept[name]=[]
    for rank,path in stages:
        start=time.monotonic();cloud=read(path);basis=tuple(tuple(map(F,p)) for p in cloud['independent_points'][:rank])
        if rank!=len(basis):raise ArithmeticError('stage basis dimension')
        if rank>17 and basis[:len(previous)]!=previous:raise ArithmeticError('historical subgroup chain is not nested by prefix')
        previous=basis
        allpts=(*basis,*targets.values());hg,asym=geo.canonical_height_gram(model,allpts);h=matrix(ZZ,geo.rounded_gram(hg,1000000));g=h[:rank,:rank]
        change=matrix(ZZ,pari(g).qflllgram()).transpose();inv=change.inverse();reduced=change*g*change.transpose()
        dd=cvp.RoundedMetricCVP(reduced,'dd',192);mp=cvp.RoundedMetricCVP(reduced,'mpfr',192)
        for ti,(name,target) in enumerate(targets.items()):
            atlas=list(atlases[name])
            for arm in arms:
                if arm['rank_before']<=rank:
                    atlas.append({'word':arm['exact_centre']['representative'],'introduced':arm['rank_before'],'orbit':{'orbit_mask':arm['degree_two_orbit_mask'],'generic_shell':arm['generic_shell']},'origin':'historical recovery centre'})
            projection=vector(QQ,g.solve_right(vector(ZZ,[2*h[i,rank+ti] for i in range(rank)])))*inv
            trials=[]
            for item in atlas:
                word=item['word']+[0]*(rank-len(item['word']));rword=vector(ZZ,word)*inv;p=vector(ZZ,[int(x)%2 for x in rword]);u=(rword-p)/2
                z,dist=cvp.cross_precision_closest(dd,mp,(projection-p)/2)
                translation=list(map(int,(u-z)*change));moved=short_add(model,target,linear_combination(model,basis,translation))
                centre=linear_combination(model,basis,word)
                record={'centre_word':word,'translation_word':translation,'orbit':item['orbit'],'introduced':item['introduced'],'parity':sum((x%2)<<i for i,x in enumerate(word)),
                        'extension_bits':sum((x%2)<<i for i,x in enumerate(word[17:])),'centre_rounded_norm':int(vector(ZZ,word)*g*vector(ZZ,word)),
                        'cvp_distance':str(dist),'origin':item['origin']}
                if moved is None or moved==centre:
                    record.update(status='EXACT_SUBGROUP_CONTAINMENT_ENDPOINT',height=None,coordinate=None)
                else:
                    key=tuple(centre)
                    if key not in mapping_cache:mapping_cache[key]=mapper.mapping(model,basis,{'representative':word})
                    class FixedMapper:
                        @staticmethod
                        def mapping(*args):return mapping_cache[key]
                    exact=helpers.coordinate(FixedMapper,model,basis,word,moved)
                    record.update(status='EXACT_QUARTIC_SQUARE',height=exact['reduced_coordinate_height'],coordinate=exact['reduced_coordinate'],
                        square_root=exact['reduced_quartic_square_root_absolute'],quartic=exact['mapping']['discriminant_quartic'],coordinate_matrix=exact['mapping']['matrix'],
                        target_point=list(map(str,target)),translated_target=list(map(str,moved)),centre_point=list(map(str,centre)),
                        quartic_max_coefficient_bits=max(abs(F(x).numerator).bit_length() for x in exact['mapping']['discriminant_quartic']))
                trials.append(record)
            # Count distinct chart/translate witnesses, not repeated stage observations.
            retained={json.dumps([r['centre_word']+[0]*(rank-len(r['centre_word'])),r['translation_word']+[0]*(rank-len(r['translation_word']))]):r for r in kept[name]+trials}
            kept[name]=list(retained.values());finite=[r for r in kept[name] if r['height'] is not None]
            winner=min(finite,key=lambda r:int(r['height'])) if finite else None
            cell={'direction':name,'rank':rank,'new_trials':len(trials),'retained_witness_count':len(kept[name]),'cheap_representatives':sum(int(r['height'])<=125000 for r in finite),
                  'contained_by_exact_endpoint':any(r['status']=='EXACT_SUBGROUP_CONTAINMENT_ENDPOINT' for r in kept[name]),'winner':winner,'trials':trials}
            cells.append(cell);checkpoint(D/f'{name}-M{rank}.json',cell)
        print('MATRIX stage',rank,'directions',len(targets),'seconds',round(time.monotonic()-start,2),flush=True)
    payload={'status':'PASS_FINITE_ATLAS_CASCADE_MATRIX','protocol_sha256':sha(D/'protocol.json'),'inputs':inputs,
             'boundary':'Minimum over the retained finite atlas and vetted translations only. Exact rational witnesses; rounded metric decisions. No global optimality, no point search.',
             'cells':[{k:v for k,v in c.items() if k!='trials'} for c in cells]}
    output=ART/'curve302_visibility_cascade_matrix_v1.json'
    if output.exists():raise FileExistsError('preserve generated matrix')
    checkpoint(output,payload)

if __name__=='__main__':run()
