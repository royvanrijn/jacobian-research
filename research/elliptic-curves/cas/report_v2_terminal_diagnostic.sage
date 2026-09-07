#!/usr/bin/env sage-python
"""Independent exact witness and common-metric replay; package terminal diagnosis."""
import hashlib
import json
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import matrix, vector, ZZ, QQ, pari
from research_runtime.store import checkpoint
from research_runtime.search_state import raw_state
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
from pointed_quartic_search import PointedQuarticSearch
from terminal_affine_cvp import AffineCVP
import pari_pointed_backend as backend

CAS=Path(__file__).resolve().parent
l=SourceFileLoader('replay_terminal_landscape',str(CAS/'v2_terminal_landscape.sage')).load_module()
b=l.b;D=b.D
data=b.read(D/'coordinates.json');plan=b.read(D/'state-plan.json');fixed=b.read(D/'fixed-M30.json')
model=tuple(map(F,fixed['curve']));points=tuple(tuple(map(F,p)) for p in fixed['points']);target=tuple(map(F,data['missing_point']))
ambient=tuple(tuple(map(F,p)) for p in data['ambient_points']);full=matrix(ZZ,plan['B31_ambient_words'])
allpts=tuple(tuple(map(F,p)) for p in plan['B31_points']);h=matrix(ZZ,plan['common_ambient_metric'])
assert abs(full.det())==1
for word,point in zip(full.rows(),allpts):assert l.group.linear_combination(model,ambient,word)==point
from memory_rank_certificate import checked_rank
cert=b.read(D/'B31-rank-certificate.json')['proof']
actual=checked_rank(model,allpts,[r['prime'] for r in cert['signatures']],cert['no_rational_2_torsion_prime'])
assert json.dumps(actual,sort_keys=True)==json.dumps(cert,sort_keys=True)
empty=raw_state(model,(),cache=Cache(MemoryFactStore()),prime_bound=1000)
checked=0


def witness(w,basis,missing):
    global checked
    assert tuple(map(F,w['target']))==missing
    cword=w['centre_word'];qword=w['closest_double_target_word'];translation=w['translation_word']
    assert all(c-q==2*z for c,q,z in zip(cword,qword,translation))
    centre=l.group.linear_combination(model,basis,cword)
    moved=l.group.short_add(model,missing,l.group.linear_combination(model,basis,translation))
    assert tuple(map(F,w['translated_target']))==moved
    n,d=map(int,w['coordinate']);assert d>=0 and gcd(abs(n),d)==1
    assert int(w['height'])==max(abs(n),d)
    mapping=w['mapping'];disc=list(map(F,mapping['discriminant_quartic']))
    square=sum(v*n**i*d**(4-i) for i,v in enumerate(disc))
    assert square==F(w['square_root'])**2
    search=PointedQuarticSearch(state=empty,centre={'point':centre},coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search,mapping)
    value=sum(v*n**i*d**(4-i) for i,v in enumerate(search.coefficients))
    assert value.denominator==1 and value>=0
    root=isqrt(value.numerator);assert root*root==value
    assert moved in {search.map_hit(n,d,root),search.map_hit(n,d,-root)}
    checked+=1


def atlas(path,basis,missing):
    allw=[]
    for trial in sorted(path.glob('trial-*.json')):
        r=b.read(trial)
        for w in r.get('witnesses',[]):witness(w,basis,missing);allw.append(w)
    result=b.read(path/'result.json')
    assert int(result['best']['height'])==min(int(w['height']) for w in allw)
    finite=[w for w in allw if int(w['coordinate'][1])]
    return {'best_projective_height':int(result['best']['height']),
        'best_finite_height':min((int(w['height']) for w in finite),default=None),
        'best':result['best'],'witness_count':len(allw),
        'global_projective_minimum_proved':int(result['best']['height'])==1}


terminal=atlas(D/'terminal-landscape',points,target)
g,u,inv,t,reduced,constant=l.setup(points,target,model)
nearest=AffineCVP(reduced.rows()).nearest(t,256,30000000)
assert nearest==b.read(D/'terminal-landscape/nearest.json')
print('REPLAYED UNRESTRICTED TERMINAL CVP AND CHART ATLAS',checked,flush=True)
for path in sorted((D/'retained-bank').glob('anchor-*.json')):
    row=b.read(path)
    for candidate in row['candidates']:
        for w in candidate['witnesses']:witness(w,points,target)
for path in sorted((D/'terminal-searched-charts').glob('chart-*.json')):
    row=b.read(path)
    assert b.sha(b.ROOT/row['chart'])==row['chart_sha256']
    for w in row['witnesses']:witness(w,points,target)
print('REPLAYED RETAINED BANK AND SEARCHED CHART WITNESSES',checked,flush=True)
selected=set(plan['selected_rank29_ids'])|set(b.read(D/'exact-predecessor-plan.json')['selected_ids'])
state_results=[]
for state in plan['states']:
    if state['id']=='omit-13':continue
    if state['rank']==29 and state['id'] not in selected:continue
    basis=tuple(allpts[i] for i in state['kept'])
    cells={str(i):atlas(D/'states'/state['id']/('target-%d'%i),basis,allpts[i]) for i in state['omitted']}
    state_results.append({'state':state,'cells':cells})
    print('REPLAYED OMISSION STATE',state['id'],flush=True)


def exact_cost(rows,missing):
    gram=rows*h*rows.transpose();cross=rows*h*missing
    change=matrix(ZZ,pari(gram).qflllgram()).transpose();reduced=change*gram*change.transpose()
    t=vector(QQ,gram.solve_right(2*cross))*change.inverse()
    orthogonal=4*(missing*h*missing)-t*reduced*t
    proof=AffineCVP(reduced.rows()).nearest(t,1,30000000)
    return int(QQ(proof['rows'][0]['distance'])+orthogonal)


common=[]
for state in plan['states']:
    old=b.read(D/'common-exact'/(state['id']+'.json'))
    for i in state['omitted']:
        assert exact_cost(full[state['kept'],:],full.row(i))==old['targets'][str(i)]['minimum_numerator']
    common.append({'state':state['id'],'rank':state['rank'],'costs':{i:r['minimum_numerator'] for i,r in old['targets'].items()}})
print('REPLAYED ALL105 EXACT COMMON-METRIC STATES',flush=True)
table=b.read(b.ART/'curve302_exceptional_subgroup_landscape_v1.json')
axis_results=[]
identity=matrix.identity(ZZ,31)
for path in sorted((D/'common-exact').glob('agent2-omit-*.json')):
    r=b.read(path);i=r['target_axis'];value=exact_cost(identity[[j for j in range(31) if j!=i],:],identity.row(i))
    assert value==r['result']['minimum_numerator']
    row=table['subset_states'][r['state_mask']]
    assert value<=row['retained_numerators'][i-17]<=row['stage_local_numerators'][i-17]
    axis_results.append({'axis':i,'state_mask':r['state_mask'],'exact_common_numerator':value,
        'agent2_retained_numerator':row['retained_numerators'][i-17]})
trace=b.read(D/'gate-trace.json')
assert trace['best_vetted_retained_bank_witness']['witness']['height']=='1'
seal=b.read(D/'V2-seal.json');assert all(b.sha(b.ROOT/k)==v for k,v in seal.items())
sources=[CAS/n for n in ('diagnose_v2_terminal.sage','terminal_affine_cvp.py','v2_terminal_landscape.sage',
    'v2_terminal_bank.sage','v2_terminal_states.sage','v2_terminal_common_metric.sage',
    'trace_v2_terminal_candidates.py','report_v2_terminal_diagnostic.sage')]
output=b.ART/'v2_terminal_failure_decomposition_v1.json'
if output.exists():raise FileExistsError('preserve terminal report')
result={'schema':'elliptic-curves.v2-terminal-failure-decomposition.v1','status':'PASS_BOUNDED_RETROSPECTIVE_DECOMPOSITION',
    'V2_unchanged_files':len(seal),'sources':{b.rel(p):b.sha(p) for p in sources},
    'fixed_B31_certificate':b.read(D/'B31-rank-certificate.json'),
    'coordinates':data,'terminal_landscape':terminal,'gate_trace':trace,'state_results':state_results,
    'all105_common_metric_states':common,'agent2_axis_M30_comparison':axis_results,
    'equivalent_lifts':b.read(D/'common-exact/equivalent-lifts.json'),'exact_witnesses_replayed':checked,
    'classification':{
        'state_geometry':'V2 M30 is worse than every coordinate-axis M30 in the exact common half-lattice metric, but has globally optimal projective AND finite chart height1. Therefore the state does not force expensive point-search exposure.',
        'selection':'A retained norm10 anchor has a scored extension with a height1 witness. It is discarded before exact CVP; it also loses the deepest-eight cutoff in the counterfactual refined set. The failure is not merely the removed64-mask cap.',
        'exposure':'The best nearest-translate witness vetted in already searched charts has height446175, beyond125000 and within the standard fourfold tier500000.',
        'not_claimed':'No billion-row parity ledger, no claim of globally minimal chart height for every coset or every alternative finite atlas. Global terminal height1 is proved by an attained witness and the universal height>=1 bound. No V3 point search has been run.'},
    'checkpoint_hashes':{b.rel(p):b.sha(p) for p in sorted(D.rglob('*.json'))}}
checkpoint(output,result)
# The user selected a shortlist-only V3. Do not emit or execute a height-tier
# proposal: the larger-box witness above is diagnostic evidence only.
print('PACKAGED TERMINAL DIAGNOSIS',checked,'witnesses; V2 sealed',len(seal),'report sha',b.sha(output),flush=True)
