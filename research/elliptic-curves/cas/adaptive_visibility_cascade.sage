#!/usr/bin/env sage-python
"""Bounded generic-only low-shell cascade; freeze precedes all discovery."""
import argparse
import csv
import hashlib
import json
import sys
import time
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
D = ROOT / 'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v1'
ORBITS = ART / 'curve302_parent_degree2_multisection_orbits_v1.tsv'
READS = set()

def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.resolve().relative_to(ROOT))
def load(name): return SourceFileLoader('cascade_'+name.replace('.','_'), str(CAS/name)).load_module()

from research_runtime.store import checkpoint

def sources():
    import pari_pointed_backend as backend
    names = ['factor_free_pari_mapping.sage','prospective_half_lattice_v3.sage',
             'audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py',
             'research_runtime/search_state.py','research_runtime/mw_state.py',
             'research_runtime/pointed_orbit_compression.py']
    return {**backend.sources(), rel(Path(__file__)):sha(Path(__file__)),
            **{rel(CAS/n):sha(CAS/n) for n in names}}

def freeze():
    if D.exists(): raise FileExistsError('preserve frozen cascade')
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    reduced=read(ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json')
    generic={k:parent[k] for k in ('a_invariants','basis_weierstrass_coordinates','generic_height_gram')}
    reduced_generic={k:reduced[k] for k in ('a_invariants','basis_weierstrass_coordinates')}
    selection=read(ROOT/'artifacts/local/elliptic-curves/det1092-record-scale-selection-v1/selection-result.json')
    # Predetermined score strata, never a post-search rank or response.
    controls=[]
    for height in (10,11):
        for stratum,count in [('strong',3),('moderate',1),('lower_fixed',1)]:
            eligible=sorted((r for r in selection['selected'] if r['height_bin']==height and r['stratum']==stratum),key=lambda r:r['id'])
            controls += [{k:r[k] for k in ('id','parameter','model','height_bin','stratum')} for r in eligible[:count]]
    assert len(controls)==10
    D.mkdir(parents=True)
    checkpoint(D/'generic.json',generic);checkpoint(D/'reduced-generic.json',reduced_generic)
    checkpoint(D/'roster.json', [{'id':'calibration302','parameter':'0','presentation':'normalized'}]+
               [{**r,'presentation':'reduced'} for r in controls])
    checkpoint(D/'protocol.json',{
        'schema':'adaptive-visibility-cascade.v1','sources':sources(),
        'inputs':{rel(p):sha(p) for p in (D/'generic.json',D/'reduced-generic.json',D/'roster.json',ORBITS)},
        'selection_calibration':'Designed after the known chain; execution and ranking use no exceptional points, residual labels or preferred orbit IDs.',
        'candidate_rule':'Enumerate every exported generic orbit of minimum norm 8 or 10. Rank canonical words separately in each shell by descending rounded specialized height norm, then orbit mask. Retain all candidates in this order. Each wave takes the next 25 unsearched canonical words per shell. Independently take the first 16 anchors per shell and their first 64 binary extension masks in the newly admitted basis coordinates, CVP reduce, and take the 48 deepest unsearched parities, with parity as tie break. No novelty floor. Previously selected charts persist.',
        'charts_per_wave_cap':98,'max_waves':4,'height':125000,'seconds_per_chart':10,
        'stop':'Stop after a complete no-gain wave, four waves, or 32 certified independent points; timeouts are censored and are not negative evidence.',
        'metric':'384-bit PARI canonical heights rounded at 10^6; floating CVP with exact integer norm/parity recheck, no optimality assertion.',
        'gp_sha256':sha(Path('/usr/bin/gp')),
        'reduction_prime_bound':1000,'scope':'One calibration fibre and ten score-stratified determinant-1092 controls; identical policy and budgets.'})
    print('FROZEN generic-only cascade: 11 fibres, 98 charts/wave, 4 waves maximum',flush=True)

def guard():
    def audit(event,args):
        if event!='open' or not args or not isinstance(args[0],(str,bytes)):return
        p=Path(args[0]).resolve()
        if not p.is_relative_to(ROOT/'artifacts'):return
        if p.is_relative_to(D) or p==ORBITS or p.is_relative_to(ROOT/'artifacts/local/elliptic-curves/pointed-sieve-build'):
            READS.add(rel(p));return
        raise PermissionError('cascade rejected artifact: '+rel(p))
    sys.addaudithook(audit)

def protocol():
    p=read(D/'protocol.json')
    if p['sources']!=sources() or any(sha(ROOT/k)!=v for k,v in p['inputs'].items()):
        raise ArithmeticError('frozen cascade source/input changed')
    return p

def seed(row):
    from sage.all import QQ, PolynomialRing
    ring=PolynomialRing(QQ,'t');t=QQ(row['parameter'])
    parent=read(D/('generic.json' if row['presentation']=='normalized' else 'reduced-generic.json'))
    def ev(v):return ring(v['numerator'])(t)/ring(v['denominator'])(t)
    aa=[ev(v) for v in parent['a_invariants']]
    if row['presentation']=='normalized':
        a1,a2,a3,a4,a6=aa;b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        model=[0,0,0,-27*(b2*b2-24*b4),-54*(-b2**3+36*b2*b4-216*b6)]
        points=[(36*ev(x)+3*b2,108*(2*ev(y)+a1*ev(x)+a3)) for x,y in parent['basis_weierstrass_coordinates']]
    else:
        q=t.denominator();model=[0,0,0,aa[3]*q**8,aa[4]*q**12]
        assert list(map(str,model))==row['model']
        points=[(ev(x)*q**4,ev(y)*q**6) for x,y in parent['basis_weierstrass_coordinates']]
    return tuple(F(str(x)) for x in model),tuple(tuple(F(str(x)) for x in pt) for pt in points)

def candidates(model,basis,seen):
    from sage.all import matrix,ZZ,pari
    import numpy as np
    geo=load('prospective_half_lattice_v3.sage')
    hg,asym=geo.canonical_height_gram(model,basis);g=matrix(ZZ,geo.rounded_gram(hg,1000000))
    if not g.is_positive_definite():raise ArithmeticError('nonpositive height decision metric')
    rank=len(basis);buckets={8:[],10:[]}
    with ORBITS.open() as stream:
        for r in csv.DictReader(stream,delimiter='\t'):
            norm=int(r['minimum_norm'])
            if norm in buckets:
                buckets[norm].append({'orbit':int(r['orbit_mask']),'shell':norm,'representative':list(map(int,r['parent_MW17_w'].split()))})
    chosen=[];anchors=[]
    for shell,rows in buckets.items():
        w=np.asarray([r['representative'] for r in rows],dtype=np.int64);gg=np.asarray(g[:17,:17].rows(),dtype=np.int64)
        if 17**2*int(abs(w).max())**2*int(abs(gg).max())>=2**62:raise ArithmeticError('integer vector norm fast path overflow')
        norms=np.einsum('ij,jk,ik->i',w,gg,w,optimize=True)
        for r,n in zip(rows,norms):r['metric_norm']=int(n)
        rows.sort(key=lambda r:(-r['metric_norm'],r['orbit']))
        anchors+=rows[:16]
        eligible=[r for r in rows if ('canonical',r['orbit']) not in seen]
        for r in eligible[:25]:
            chosen.append({**r,'representative':r['representative']+[0]*(rank-17),'lane':'canonical'})
    if rank>17:
        u=matrix(ZZ,pari(g).qflllgram()).transpose();inv=u.inverse();oracle=geo.CosetOracle((u*g*u.transpose()).rows())
        extensions=[]
        for anchor in anchors:
            for ext in range(min(64,2**(rank-17))):
                residue=[x%2 for x in anchor['representative']]+[(ext>>i)&1 for i in range(rank-17)]
                mask=sum(x<<i for i,x in enumerate(residue))
                if ('parity',mask) in seen:continue
                rp=[int(x)%2 for x in (matrix(ZZ,1,rank,residue)*inv).row(0)]
                norm,rep,error=oracle.solve(rp)
                word=list(map(int,(matrix(ZZ,1,rank,rep)*u).row(0)))
                if [x%2 for x in word]!=residue:raise ArithmeticError('CVP parity transport failed')
                extensions.append({'orbit':anchor['orbit'],'shell':anchor['shell'],'extension':ext,'parity':mask,'representative':word,'metric_norm':norm,'lane':'parity'})
        extensions.sort(key=lambda r:(-r['metric_norm'],r['parity']))
        chosen+=extensions[:48]
    return chosen,{'rounded_gram':[list(map(int,r)) for r in g.rows()],'shell_counts':{str(k):len(v) for k,v in buckets.items()},'asymmetry':str(asym)}

def run(index):
    guard();p=protocol();row=read(D/'roster.json')[index];out=D/row['id'];out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists():print('Already terminal',row['id']);return
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from research_runtime.search_state import raw_state
    from research_runtime.pointed_orbit_compression import compress
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl
    model,basis=seed(row);cache=Cache(MemoryFactStore());state=raw_state(model,basis,cache=cache,prime_bound=1000)
    if state.rank!=17 or tuple(tuple(map(F,x)) for x in state.basis)!=basis:raise ArithmeticError('generic seed rank failed')
    checkpoint(out/'seed.json',{'curve':list(map(str,model)),'points':[list(map(str,pt)) for pt in basis],'state':state.record()})
    mapper=load('factor_free_pari_mapping.sage');mapper.pari.allocatemem(256000000,silent=True)
    seen=set();stages=[]
    for wave in range(p['max_waves']):
        wd=out/f'wave-{wave:02d}';wd.mkdir(exist_ok=True)
        if (wd/'result.json').exists():raise FileExistsError('incomplete run requires explicit checkpoint recovery')
        before=state.rank;basis=tuple(tuple(map(F,x)) for x in state.basis);centres,metric=candidates(model,basis,seen)
        checkpoint(wd/'selection.json',{'protocol_sha256':sha(D/'protocol.json'),'rank':before,'basis':[list(map(str,x)) for x in basis],'centres':centres,'metric':metric})
        charts=[];start=time.monotonic()
        for j,c in enumerate(centres):
            mapping=mapper.mapping(model,basis,c)
            search=PointedQuarticSearch(state=state,centre={'coefficients':c['representative']+[0]*(state.rank-before)},coordinate_policy=mapping['coordinate_policy'])
            # The centre point is fixed on the wave basis; padding keeps it fixed after admissions.
            transcript,points=backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
            if backend.replay(search,mapping,transcript)!=points:raise ArithmeticError('exact map replay differs')
            compression=compress(model,state.basis,c['representative']+[0]*(state.rank-before),points)
            charts.append({'index':j,'centre':c,'mapping':mapping,'search':transcript,'admission_compression':compression})
            for k in compression['kept_indices']:state=state.adjoin(points[k],cache=cache)
            seen.add((c['lane'],c['orbit'] if c['lane']=='canonical' else c['parity']))
            checkpoint(wd/f'chart-{j:03d}.json',charts[-1])
            if (j+1)%10==0:print(row['id'],'wave',wave,'chart',j+1,'rank',state.rank,flush=True)
        result={'status':'COMPLETE_DECLARED_WAVE','family':'det1092-cascade','parameter':row['parameter'],'curve':list(map(str,model)),
                'charts':charts,'final_state':state.record(),'rank_lower_bound':state.rank}
        checkpoint(wd/'result.json',result)
        mod2.build(wd/'result.json',wd/'mod2.json',1000,sha(wd/'result.json'));mod2.check(wd/'mod2.json')
        modl.build(wd/'mod2.json',wd/'modl.json');modl.check(wd/'modl.json')
        cloud=read(wd/'mod2.json');basis=tuple(tuple(map(F,x)) for x in cloud['independent_points'])
        state=raw_state(model,basis,cache=cache,prime_bound=1000)
        if state.rank!=cloud['rank_lower_bound']:raise ArithmeticError('recertified seed mismatch')
        stage={'wave':wave,'before':before,'after':state.rank,'charts':len(charts),'completed':sum(c['search']['status']=='bounded_search_complete' for c in charts),'seconds':time.monotonic()-start,'mod2':rel(wd/'mod2.json'),'mod2_sha256':sha(wd/'mod2.json'),'modl':rel(wd/'modl.json'),'modl_sha256':sha(wd/'modl.json')}
        stages.append(stage);checkpoint(out/'stages.json',stages)
        print('CASCADE STAGE',row['id'],stage,flush=True)
        if state.rank==before or state.rank>=32:break
    checkpoint(out/'terminal.json',{'status':'BOUNDED_CASCADE_TERMINAL','id':row['id'],'protocol_sha256':sha(D/'protocol.json'),'stages':stages,'final_rank_lower_bound':state.rank,'read_paths':sorted(READS)})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','run']);ap.add_argument('--index',type=int,default=0);a=ap.parse_args()
    freeze() if a.stage=='freeze' else run(a.index)
