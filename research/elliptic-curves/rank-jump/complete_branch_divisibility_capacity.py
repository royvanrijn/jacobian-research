#!/usr/bin/env python3
"""Complete the fixed branch panel using a retained generic-basis calibration."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import branch_blocks
import root_curve_real_components as real

PROTOCOL=Path(__file__).with_name('BRANCH_DIVISIBILITY_CAPACITY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_branch_divisibility_capacity_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_branch_divisibility_capacity_v2.json'
WORK=r.ROOT/'artifacts/local/rank-jump-branch-divisibility-capacity-v2'
MODEL=r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
SECTIONS=r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'
ATLAS=r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    from sage.all import QQ,PolynomialRing
    old=r.read(branch_blocks.INPUT);model=r.read(MODEL);sections=r.read(SECTIONS)['sections']
    R=PolynomialRing(QQ,'t');points=[]
    for section in sections:
        x=R(section['x_coefficients_low_to_high'])
        if 'y_coefficients_low_to_high' in section:y=R(section['y_coefficients_low_to_high'])
        else:
            chord=section['chord'];xx,yy=points[chord['reference_basis_index']]
            y=yy+R(chord['slope_coefficients_low_to_high'])*(x-xx)
        points.append((x,y))
    wanted={row['label'] for row in old['covers']};atlas={x['label']:x for x in r.read(ATLAS)['bisections'] if x['label'] in wanted}
    covers=[]
    for row in old['covers']:
        q=row['residual_chord']['q_coefficients']
        assert q==atlas[row['label']]['residual_chord']['q_coefficients']
        covers.append({'label':row['label'],'q':q})
    assert len(covers)==len(wanted)==37
    r.write_new(INPUT,{'schema':'rank-jump.branch-divisibility-capacity-inputs.v1',
        'A':model['A_coefficients_low_to_high'],'B':model['B_coefficients_low_to_high'],
        'sections':[{'x':list(map(str,x.list())),'y':list(map(str,y.list()))} for x,y in points],
        'covers':covers,'bindings':bindings([Path(__file__),PROTOCOL,MODEL,SECTIONS,ATLAS,branch_blocks.INPUT]),
        'boundary':'Equation/generic-only arithmetic projection. Retrospective cover selection, not a rank predictor. No exceptional coordinates, lift maps, quotient labels or curve ranks.'})
    print('Exported',len(covers),'native branch supports and',len(sections),'generic sections')


def compute():
    from sage.all import QQ,GF,PolynomialRing
    data=r.read(INPUT);R=PolynomialRing(QQ,'t');A=R(data['A']);B=R(data['B']);D=-4*A**3-27*B**2
    sections=[(R(s['x']),R(s['y'])) for s in data['sections']]
    assert len(sections)==17 and all(y*y==x**3+A*x+B for x,y in sections)
    assert [A.degree(),B.degree(),D.degree()]==[8,12,24] and D.is_squarefree() and D.gcd(A)==1
    assert D.leading_coefficient()!=0
    polys=[R(c['q']) for c in data['covers']]
    assert all(q.degree()==2 and not q.discriminant().is_square() and q.gcd(D)==1 for q in polys)
    # Exact real root-curve topology; the independent replay checks isolation
    # and recounts connected components by a separate graph traversal.
    events=[]
    for (a,b),multiplicity in (D*B.squarefree_part()).real_root_intervals():
        assert multiplicity==1 and a<b
        if D(a)*D(b)>0:continue
        assert D(a)*D(b)<0 and B(a)*B(b)>0
        events.append({'interval':[str(a),str(b)],'B_sign':int(B(a).sign()),'type':'I1'})
    samples=[(QQ(e['interval'][1])+QQ(events[i+1]['interval'][0]))/2 if i+1<len(events) else QQ(e['interval'][1])+1 for i,e in enumerate(events)]
    arities=[3 if D(t)>0 else 1 for t in samples];edges,groups=real.graph(events,arities)
    cap=10+len(groups)-1
    topology={'genus':10,'events':events,'samples':list(map(str,samples)),'arities':arities,
        'real_components':len(groups),'global_pool_upper_bound':cap}
    rows=[{'label':c['label'],'irreducibility_witness':None,'blocks':[]} for c in data['covers']]
    generic_blocks=[];generic_signatures=[0]*17
    # All reductions are in the fixed order. There is no rank-dependent stop.
    for p in r.primes(2003):
        F=GF(p);P=PolynomialRing(F,'X');X=P.gen()
        try:
            ap,bp=P(A),P(B);qs=[P(q) for q in polys];sp=[(P(x),P(y)) for x,y in sections]
        except (ValueError,ZeroDivisionError):continue
        def block(t):
            a,b=ap(t),bp(t);f=X**3+a*X+b
            if not f.discriminant():return None,f
            roots=f.roots(multiplicities=False)
            if len(roots)!=3:return None,f
            xs=[x(t) for x,y in sp];ys=[y(t) for x,y in sp]
            if any(not y for y in ys):return None,f
            roots=sorted(map(int,roots));signatures=[]
            for x in xs:
                bits=[int(pow((int(x)-z)%p,(p-1)//2,p)==p-1) for z in roots]
                assert sum(bits)%2==0;signatures.append(r.pack(bits))
            return {'p':p,'base_root':int(t),'cubic_roots':roots,'signatures':signatures},f
        b0,f0=block(F(QQ(-2)/377)) if p not in (13,29) else (None,None)
        if b0:
            k=3*len(generic_blocks);generic_blocks.append(b0)
            for j,v in enumerate(b0['signatures']):generic_signatures[j]|=v<<k
        for q,row in zip(qs,rows):
            if q.degree()!=2 or not q.discriminant():continue
            for t in sorted(q.roots(multiplicities=False),key=int):
                bl,f=block(t)
                if f.discriminant() and f.is_irreducible() and row['irreducibility_witness'] is None:
                    row['irreducibility_witness']={'p':p,'base_root':int(t),'cubic_ascending':list(map(int,f.list()))}
                if bl:row['blocks'].append(bl)
    m=r.rank(generic_signatures);assert m==17
    for row in rows:
        signatures=[sum(bl['signatures'][j]<<(3*i) for i,bl in enumerate(row['blocks'])) for j in range(17)]
        rank=r.rank(signatures);row['branch_character_rank']=rank
        row['branch_divisibility_kernel_upper_bound']=17-rank
        row['status']='PASS' if row['irreducibility_witness'] else 'UNKNOWN'
        row['all_scalar_twists_rank_upper_bound']=cap-rank if row['status']=='PASS' else 'UNKNOWN'
    return {'schema':'rank-jump.branch-divisibility-capacity.v2','status':'PASS','topology':topology,
        'generic_mod_two_dimension':m,'generic_independence_parameter':'-2/377','generic_independence_blocks':generic_blocks,
        'rows':rows,'bindings':bindings([Path(__file__),PROTOCOL,INPUT,Path(real.__file__),Path(r.__file__),Path(__file__).with_name('BRANCH_DIVISIBILITY_COMPLETION_PROTOCOL.json'),Path(__file__).with_name('branch_divisibility_capacity.py')]),
        'boundary':'Arithmetic generic twist-rank upper bounds for all nonzero rational scalar multiples of each branch polynomial. No bound on further jumps of an individual fibre; no new point or class is constructed.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        error=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=60)
                if p.returncode:error='worker failure'
            except subprocess.TimeoutExpired:error='bounded timeout'
        if error:r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error});return
    result=r.read(path);r.write_new(OUTPUT,result)
    print('Global capacity',result['topology']['global_pool_upper_bound'],'generic dimension',result['generic_mod_two_dimension'])
    print([(x['label'],x['branch_character_rank'],x['all_scalar_twists_rank_upper_bound']) for x in result['rows']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS branch-divisibility replay')
    else:globals()[a.mode]()
