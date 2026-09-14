#!/usr/bin/env python3
"""Collect retained independent point packets; no point search or model selection.

Source certificates are retained evidence, not relabelled as fresh rank replays.
Exact short-model transport and membership are checked here. Generic-prefix
claims require explicit generic points or a packet's generic_rank convention.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import gzip
import json
from pathlib import Path
import sqlite3

from analyze_seed_amplification_history import LEGACY
from half_lattice_pointed_sieve import integral_short_scale, linear_combination

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
OUT = ART/'finite_cancellation_predictor_v1'

def digest(x): return sha256(x).hexdigest()
def canonical(x): return json.dumps(x, sort_keys=True, separators=(',', ':')).encode()
def write(path, x):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(x, indent=2, sort_keys=True)+'\n')
def short(model, points):
    a1,a2,a3,a4,a6 = map(F,model)
    c2=a2+a1*a1/4;c4=a4+a1*a3/2;c6=a6+a3*a3/4
    shift=c2/3
    curve=(F(0),F(0),F(0),c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27)
    u,curve=integral_short_scale(curve)
    result=[]
    for P in points:
        if len(P)==3:
            if F(P[2])!=1:raise ArithmeticError('unsupported non-affine projective source')
            P=P[:2]
        x,y=map(F,P)
        if y*y+a1*x*y+a3*y != x**3+a2*x*x+a4*x+a6:
            raise ArithmeticError('source point is off its displayed curve')
        xx=(x+shift)*u*u;yy=(y+(a1*x+a3)/2)*u**3
        if yy*yy != xx**3+curve[3]*xx+curve[4]:raise ArithmeticError('short transport failed')
        result.append([str(xx),str(yy)])
    return list(map(str,curve)),result
def pointkey(p): return str(F(p[0]))+','+str(abs(F(p[1])))

class Corpus:
    def __init__(self):
        self.inputs={};self.gaps=[];self.packets=[];self.seen=set();self.registry=[];self.generic_banks=[];self.pending=[]
    def read(self,p,expected=None):
        p=Path(p)
        if not p.exists():
            self.gaps.append({'path':str(p.relative_to(ROOT)), 'reason':'MISSING_RETAINED_SOURCE'});return None
        b=p.read_bytes();h=digest(b)
        if expected and h!=expected:raise ArithmeticError(f'changed input: {p}')
        self.inputs[str(p.relative_to(ROOT))]=h
        return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
    def add(self,p,source,family=None,generic=None):
        model=p.get('curve',p.get('ainvs',p.get('model')))
        points=p.get('points')
        proof=p.get('rank_certificate',p.get('proof'))
        if not isinstance(model,list) or len(model)!=5 or not isinstance(points,list) or not points or not isinstance(proof,dict):return
        if any(not isinstance(z,list) or len(z)!=2 for z in points):return
        rank=proof.get('rank_lower_bound',p.get('rank_lower_bound'))
        if rank!=len(points):
            self.gaps.append({'source':source,'reason':'POINT_LIST_NOT_FULL_INDEPENDENT_CERTIFICATE','rank':rank,'points':len(points)});return
        gen=p.get('generic_points',generic)
        n=p.get('generic_rank')
        if gen is None and isinstance(n,int):gen=points[:n]
        if gen is None and p.get('discovery_witness'):
            w=p['discovery_witness']
            try:raw=self.read(ROOT/w['path'],w.get('sha256'))
            except ArithmeticError:
                # A live search can have advanced after an interim export.
                # Preserve the broken historical binding; use current bytes
                # solely as a generic-subset candidate, checked below against
                # the frozen independent points. Do not retag the certificate.
                self.gaps.append({'source':source,'path':w['path'],'expected_sha256':w.get('sha256'),
                    'reason':'STALE_INTERIM_WITNESS_CURRENT_GENERIC_SUBSET_ONLY'})
                raw=self.read(ROOT/w['path'])
            if raw is not None and raw.get('curve')==model:gen=raw.get('generic_points')
        if gen is None:
            self.pending.append((p,source,family));return
        if not isinstance(gen,list) or not gen or not isinstance(gen[0],list):return
        cm,cp=short(model,points);gm,gp=short(model,gen)
        assert cm==gm
        # Explicit generic points must be contained in this independent list.
        ps={pointkey(z) for z in cp};gs={pointkey(z) for z in gp}
        if len(gs)!=len(gp) or not gs<=ps:
            self.gaps.append({'source':source,'reason':'GENERIC_SUBGROUP_NOT_A_CERTIFIED_SUBSET'});return
        targets=[z for z in cp if pointkey(z) not in gs]
        self.generic_banks.append({'curve':cm,'generic_points':gp,'generic_rank':len(gp),'source':source})
        if not targets:return
        key=digest(canonical([cm,gp,cp]))
        if key in self.seen:
            if family:
                for previous in self.packets:
                    if previous['id']==key[:20] and previous['family']=='UNCLASSIFIED':previous['family']=family
            return
        self.seen.add(key)
        A,B=map(F,cm[3:]);j=1728*4*A**3/(4*A**3+27*B*B)
        self.packets.append({'id':key[:20],'family':family or p.get('family','UNCLASSIFIED'),
            'curve':cm,'generic_points':gp,'targets':targets,'source':source,
            'generic_rank':len(gp),'rank_lower_bound':rank,'j_group':digest(str(j).encode()),
            'proof_sha256':digest(canonical(proof)),
            'evidence':'RETAINED_INDEPENDENT_CERTIFICATE; exact model transport and membership checked'})
    def walk(self,d,source,family=None):
        if isinstance(d,dict):
            fam=d.get('family',family)
            model=d.get('curve',d.get('elliptic_a_invariants'))
            if model and d.get('generic_points'):
                cm,gp=short(model,d['generic_points'])
                self.generic_banks.append({'curve':cm,'generic_points':gp,'generic_rank':len(gp),'source':source})
            self.add(d,source,fam)
            for k,v in d.items():
                if k not in ['bindings','evidence','points','generic_points','rank_certificate','proof','signatures','baseline','height_gram']:
                    if isinstance(v,(list,dict)):self.walk(v,source+'/'+k,fam)
        elif isinstance(d,list):
            for i,v in enumerate(d):
                if isinstance(v,dict):self.walk(v,source+'/'+str(i),family)
    def file(self,p,expected=None):
        d=self.read(p,expected)
        if d is not None:self.walk(d,str(Path(p).relative_to(ROOT)))
        return d

def extract():
    c=Corpus();c.read(OUT/'protocol.json')
    # Full cohorts, including certified gains below the rank22 inventory cutoff.
    for name in LEGACY:c.file(ART/(name+'_results_v1.json'))
    c.file(ART/'r17_60_panel_results_v1.json')
    folder=LOCAL/'broad-rank-v1/runtime/research'
    for p in sorted((folder/'broad-cases').glob('*/state.json')):
        s=c.read(p)
        if s.get('rank',0)>17 and s.get('packet'):
            p=folder/s['packet'];d=c.read(p,s.get('packet_sha256'))
            if d:c.walk(d,str(p.relative_to(ROOT)),s['parent_id'])
    # Adaptive foundry endpoints are descriptive population, not unbiased trials.
    db=LOCAL/'high-rank-foundry-v3/ledger.sqlite';b=db.read_bytes();c.inputs[str(db.relative_to(ROOT))]=digest(b)
    con=sqlite3.connect(f'file:{db}?mode=ro',uri=True)
    rows=[json.loads(r[0]) for r in con.execute('select record from curves order by id')];con.close()
    if digest(db.read_bytes())!=digest(b):raise ArithmeticError('foundry ledger changed')
    foundry_roots=sorted(LOCAL.glob('high-rank-foundry-v*/runtime/research'),reverse=True)
    for r in rows:
        if r.get('rank') and r['rank']>17 and r.get('packet'):
            p=next((f/r['packet'] for f in foundry_roots if (f/r['packet']).exists()),foundry_roots[0]/r['packet'])
            d=c.read(p,r.get('packet_sha256'))
            if d:c.walk(d,str(p.relative_to(ROOT)),r['family'])
    p=LOCAL/'parent-foundry-v4/state.json';state=c.read(p)
    roots=sorted(LOCAL.glob('parent-foundry-v*/runtime/research'),reverse=True)
    for r in state['fibres'].values():
        if r.get('rank') and r['rank']>r['generic_rank'] and r.get('packet'):
            p=next((f/r['packet'] for f in roots if (f/r['packet']).exists()),roots[0]/r['packet'])
            d=c.read(p,r.get('packet_sha256'))
            if d:c.walk(d,str(p.relative_to(ROOT)),r['family'])
    # Current registry and every rank source; sources can supply explicit prefixes.
    dbp=ROOT/'elliptic-curves/data/research_curves/database.json';registry=c.read(dbp)['curves']
    for name in sorted({r['rank_source_certificate'] for r in registry}):c.file(ROOT/name)
    # Two nonstandard exact certificate schemas in the current registry.
    p=ART/'det1092_funnel_first_seeds_v1/rank21/standalone-proof.json';d=c.read(p)
    assert d['inherited_rank']==17 and d['matrix_rank']==len(d['points'])==21
    c.add({**d,'proof':{'rank_lower_bound':21,'retained_finite_group_certificate':digest(canonical(d))},
           'generic_points':d['points'][:17]},str(p.relative_to(ROOT)),'det1092')
    p=ART/'det1092_bifibration_conic_sources_v1/seed-generic-frame.json';d=c.read(p)
    assert d['inherited_rank']==len(d['basis'])==17
    cm,gp=short(d['curve'],d['basis'])
    c.generic_banks.append({'curve':cm,'generic_points':gp,'generic_rank':17,'source':str(p.relative_to(ROOT))})
    # Reuse the explicitly retained native17 seed packets behind early fresh-six
    # followups, whose final packets omit generic_rank metadata.
    for folder in sorted(LOCAL.glob('fresh6-*')):
        for p in sorted(folder.glob('*/seed-M17.json')):c.file(p)
        for p in sorted(folder.glob('*preparation/seed-M18.json')):c.file(p)
    pending=c.pending;c.pending=[]
    bank_index=defaultdict(list)
    for g in c.generic_banks:bank_index[tuple(g['curve'])].append(g)
    for p,source,fam in pending:
        model=p.get('curve',p.get('ainvs',p.get('model')));cm,cp=short(model,p['points']);keys={pointkey(z) for z in cp}
        candidates=[g for g in bank_index[tuple(cm)] if {pointkey(z) for z in g['generic_points']}<=keys]
        if candidates:
            g=min(candidates,key=lambda x:(x['generic_rank'],x['source']))
            z={**p,'curve':cm,'points':cp,'generic_points':g['generic_points']};c.add(z,source,fam)
    # Join displayed registry points to already documented generic subgroups on
    # the exact same normalized model. Do not assume that its first17 are generic.
    bymodel=defaultdict(list)
    for p in c.generic_banks:bymodel[tuple(p['curve'])].append(p)
    for r in registry:
        cm,cp=short(r['ainvs'],r['points']);keys={pointkey(z) for z in cp}
        candidates=[p for p in bymodel[tuple(cm)] if {pointkey(z) for z in p['generic_points']}<=keys]
        if not candidates:
            c.gaps.append({'id':r['id'],'source':r['rank_source_certificate'],'reason':'REGISTRY_GENERIC_JOIN_UNAVAILABLE'});continue
        parent=min(candidates,key=lambda p:(p['generic_rank'],p['source']))
        p={'curve':cm,'points':cp,'generic_points':parent['generic_points'],'rank_certificate':r['rank_certificate']}
        c.add(p,str(dbp.relative_to(ROOT))+'/'+r['id'],r['family'])
        c.registry.append(r['id'])
    # Public302: use the documented generic embedding; the complement is a
    # fixed unimodular completion, not an assumed prefix of its public D basis.
    import icarm_curve302 as curve302
    p=ART/'curve302_recovered_mw17_parent_v1.json';parent=c.read(p)
    model=curve302.short_coefficients();points=curve302.SHORT_POINTS
    emb=parent['basis_embedding_in_public_D']
    gen=[list(map(str,linear_combination(model,points,[row[i] for row in emb]))) for i in range(17)]
    ix=[0,2,3,5,8,12,13,16,18,21,25,28,29,30]
    # Verify the unimodular completion independently over the integers.
    import sympy as sp
    M=sp.Matrix(emb).row_join(sp.eye(31)[:,ix]);assert abs(M.det())==1
    certpath=ART/'icarm_curve302_rank31_v1.json.gz';cert=c.read(certpath)
    c.inputs['elliptic-curves/cas/icarm_curve302.py']=digest((ROOT/'elliptic-curves/cas/icarm_curve302.py').read_bytes())
    cp=[list(map(str,points[i])) for i in ix]
    c.add({'curve':list(map(str,model)),'points':gen+cp,'generic_points':gen,
           'proof':{'rank_lower_bound':31,'public_certificate_sha256':c.inputs[str(certpath.relative_to(ROOT))],
                    'unimodular_completion_determinant':str(M.det())}},'public302-completed-generic-basis','Curve302-development')
    # Three other public rank29 controls with retained generic17 identifications.
    import elkies_klagsbrun_rank29 as ek29
    import icarm_curve356 as c356
    controlpath=ROOT/'elliptic-curves/data/half_lattice_rank29_control_inputs_v1.json'
    control=c.read(controlpath)['cases']
    q12path=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json'
    q12=c.read(q12path)
    projection=c.read(ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json')
    rec385=next(r for r in projection['records'] if r['id']==385)
    model385,points385=short(rec385['ainvs'],rec385['points'])
    allpublic=[(list(map(str,ek29.short_weierstrass_coefficients())),[list(map(str,p)) for p in ek29.published_short_points()]),
               (list(map(str,c356.short_coefficients())),[list(map(str,p)) for p in c356.SHORT_POINTS]),(model385,points385)]
    for i,(ctl,(model,public)) in enumerate(zip(control,allpublic)):
        cm,cp=short(model,public);gm,gp=short(ctl['short_model'],ctl['generic_points']);assert cm==gm
        if i==0:
            ix=[int(label[1:])-1 for label in q12['displayed_exceptional_quotient']['free_basis_modulo_specialized_generic']]
            M=sp.Matrix(q12['specialized_generic_subgroup']['coordinate_matrix_rows_in_ordered_29_public_points']).row_join(sp.eye(29)[:,ix])
            assert abs(M.det())==1
            points=gp+[cp[j] for j in ix]
        else:
            assert gp==cp[:17];points=cp
        c.add({'curve':cm,'points':points,'generic_points':gp,'proof':{'rank_lower_bound':29,
            'retained_public_rank29_control':ctl['label'],'generic_identification_source':str(controlpath.relative_to(ROOT))}},
            'public-rank29/'+ctl['label'],'11952' if i==0 else '074d9')
    # Merge exact point duplicates only for identical ordered generic banks.
    groups={}
    for p in c.packets:
        k=digest(canonical([p['curve'],p['generic_points']]))[:20]
        if k not in groups:groups[k]={**p,'id':k,'targets':{},'sources':[], 'families':set()}
        g=groups[k];g['sources'].append({'source':p['source'],'proof_sha256':p['proof_sha256']})
        g['families'].add(p['family'])
        for t in p['targets']:g['targets'].setdefault(pointkey(t),t)
    records=[]
    for g in groups.values():
        g['targets']=list(g['targets'].values());g['families']=sorted(g['families'])
        g['family']=next((f for f in g['families'] if f!='UNCLASSIFIED'),'UNCLASSIFIED')
        if 'Curve302-development' in g['families']:g['family']='Curve302-development'
        records.append(g)
    records.sort(key=lambda x:x['id'])
    OUT.mkdir(parents=True,exist_ok=True)
    with gzip.GzipFile(str(OUT/'corpus.json.gz'),'wb',mtime=0) as f:f.write(canonical(records))
    write(OUT/'census.json',{'status':'PASS_RETAINED_PACKET_CENSUS','packets':len(c.packets),
        'curve_generic_banks':len(records),'j_groups':len({r['j_group'] for r in records}),
        'point_observations':sum(len(r['targets']) for r in records),
        'families':dict(Counter(r['family'] for r in records)),
        'registry_endpoints':len(registry),'registry_generic_joins':len(c.registry),
        'registry_joined_ids':c.registry,'inputs':c.inputs,'gaps':c.gaps,
        'corpus_sha256':digest((OUT/'corpus.json.gz').read_bytes()),
        'boundary':'Complete for enumerated available endpoint sources, with explicit join/missing gaps. Union points and repeated bases are not mutually independent observations. No missing points or generic subgroup reconstructed.'})
    print(json.dumps({k:v for k,v in json.loads((OUT/'census.json').read_text()).items() if k not in ['inputs','gaps','registry_joined_ids']}))

if __name__=='__main__':extract()
