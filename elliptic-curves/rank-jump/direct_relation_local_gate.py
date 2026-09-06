#!/usr/bin/env python3
"""Complete projective finite-field test from native equations only."""
import argparse
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'DIRECT_RELATION_LOCAL_GATE_PROTOCOL.json'
BASE=r.OUT/'rank_jump_native_triple_intersection_inputs_v1.json'
FORMS=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
SELECTION=r.OUT/'rank_jump_triple_translate_selection_v1.json'
INPUT=r.OUT/'rank_jump_direct_relation_local_gate_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_direct_relation_local_gate_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-direct-relation-local-gate-v1'


def capture():
    base=r.read(BASE);forms=r.read(FORMS)
    roster={c['label']:c for case in forms['cases'] for c in case['covers']}
    lifts=[]
    for c in base['covers']:
        primitive=roster[c['label']];scale=F(primitive['removed_rational_square_root'])
        assert all(F(x)==scale*scale*y for x,y in zip(c['residual_chord']['q_coefficients'],primitive['form'],strict=True))
        lift={'label':c['label'],'q':primitive['form']}
        for name in ('x0','x1','y0','y1'):
            lift[name]=[str(F(x)*(scale if name.endswith('1') else 1)) for x in c['lifted_section'][name+'_coefficients']]
        lifts.append(lift)
    r.write_new(INPUT,{'schema':'rank-jump.direct-relation-local-gate-inputs.v1',
        'A':base['A'],'B':base['B'],'sections':base['sections'],'lifts':lifts,
        'generic_words':r.read(SELECTION)['selected_words'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,BASE,FORMS,SELECTION)},
        'boundary':'Detector inputs contain no characteristic-zero intersection polynomial, exceptional point, original rational parameter or observed rank.'})


def modq(x,p):
    x=F(x);return x.numerator*pow(x.denominator,-1,p)%p


def evaluate(poly,t,p,weight):
    if t is None:return poly[weight] if len(poly)>weight else 0
    value=0
    for x in reversed(poly):value=(value*t+x)%p
    return value


def add(A,P,Q,p):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    if x==u:m=(3*x*x+A)*pow(2*y,-1,p)%p
    else:m=(v-y)*pow(u-x,-1,p)%p
    xx=(m*m-x-u)%p;return (xx,(m*(x-xx)-y)%p)


def neg(P,p):return None if P is None else (P[0],-P[1]%p)


def compute():
    from sage.all import QQ,GF,PolynomialRing
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    out=[]
    for p in r.read(PROTOCOL)['limits']['fixed_primes']:
        R=PolynomialRing(GF(p),'t')
        try:
            A=[modq(x,p) for x in inp['A']];B=[modq(x,p) for x in inp['B']]
            lifts=[{'label':c['label'],**{k:[modq(x,p) for x in c[k]] for k in ('q','x0','x1','y0','y1')}} for c in inp['lifts']]
            basis=[]
            for c in inp['sections']:
                x=R([modq(z,p) for z in c['x_coefficients_low_to_high']])
                if 'y_coefficients_low_to_high' in c:y=R([modq(z,p) for z in c['y_coefficients_low_to_high']])
                else:
                    ch=c['chord'];ref=basis[ch['reference_basis_index']]
                    y=R(ref[1])+R([modq(z,p) for z in ch['slope_coefficients_low_to_high']])*(x-R(ref[0]))
                assert y*y==x**3+R(A)*x+R(B)
                basis.append(([int(z) for z in x.list()],[int(z) for z in y.list()]))
        except (ValueError,ZeroDivisionError):
            out.append({'prime':p,'status':'UNKNOWN_NONINTEGRAL_COEFFICIENT'});continue
        delta=-16*(4*R(A)**3+27*R(B)**2);qs=[R(c['q']) for c in lifts]
        gate={'discriminant_degree':int(delta.degree()),'discriminant_squarefree':bool(delta.is_squarefree()),
              'infinity_smooth':bool(4*evaluate(A,None,p,8)**3+27*evaluate(B,None,p,12)**2),
              'cover_degrees':[int(q.degree()) for q in qs],'cover_squarefree':[bool(q.is_squarefree()) for q in qs],
              'cover_discriminant_gcd_degrees':[int(q.gcd(delta).degree()) for q in qs],
              'pairwise_branch_gcd_degrees':[int(qs[i].gcd(qs[j]).degree()) for i in range(3) for j in range(i)]}
        gate['infinity_smooth']=bool((4*evaluate(A,None,p,8)**3+27*evaluate(B,None,p,12)**2)%p)
        ok=gate['discriminant_degree']==24 and gate['discriminant_squarefree'] and gate['infinity_smooth'] and gate['cover_degrees']==[2]*3 and all(gate['cover_squarefree']) and not any(gate['cover_discriminant_gcd_degrees']+gate['pairwise_branch_gcd_degrees'])
        if not ok:out.append({'prime':p,'status':'UNKNOWN_BAD_GEOMETRY','geometry_gate':gate});continue
        for c,q in zip(lifts,qs,strict=True):
            x0,x1,y0,y1=[R(c[k]) for k in ('x0','x1','y0','y1')]
            assert y0*y0+q*y1*y1==x0**3+3*q*x0*x1*x1+R(A)*x0+R(B)
            assert 2*y0*y1==3*x0*x0*x1+q*x1**3+R(A)*x1
        squares={z:[u for u in range(p) if u*u%p==z] for z in range(p)}
        rows=[];hits=[[] for _ in inp['generic_words']];carrier_count=0
        for t in list(range(p))+[None]:
            a,b=evaluate(A,t,p,8),evaluate(B,t,p,12)
            smooth=(4*a**3+27*b*b)%p!=0
            def point(x,y):
                x%=p;y%=p;assert (y*y-x*x*x-a*x-b)%p==0
                assert y or (3*x*x+a)%p,'A native image meets the nodal singular point'
                return x,y
            pts=[point(evaluate(x,t,p,4),evaluate(y,t,p,6)) for x,y in basis]
            translates=[]
            for word in inp['generic_words']:
                S=None
                for n,P in zip(word,pts,strict=True):
                    for _ in range(abs(n)):S=add(a,S,P if n>0 else neg(P,p),p)
                translates.append(S)
            roots=[squares[evaluate(c['q'],t,p,2)] for c in lifts]
            choices=list(product(*roots));carrier_count+=len(choices)
            for us in choices:
                ps=[point(evaluate(c['x0'],t,p,4)+u*evaluate(c['x1'],t,p,3),evaluate(c['y0'],t,p,6)+u*evaluate(c['y1'],t,p,5)) for c,u in zip(lifts,us,strict=True)]
                value=add(a,add(a,ps[0],neg(ps[1],p),p),ps[2],p)
                for i,S in enumerate(translates):
                    if value==S:hits[i].append({'base_t':'infinity' if t is None else t,'roots':list(us),'smooth_elliptic_fibre':smooth})
            rows.append({'base_t':'infinity' if t is None else t,'smooth_elliptic_fibre':smooth,'root_choice_count':len(choices)})
        out.append({'prime':p,'status':'PASS','geometry_gate':gate,'projective_base_fibres':p+1,
            'carrier_Fp_point_count':carrier_count,'base_fibre_counts':rows,
            'relations':[{'index':i,'Fp_points':hs,'Fp_point_count':len(hs),
                'local_obstruction_proved':not hs} for i,hs in enumerate(hits)]})
    return {'schema':'rank-jump.direct-relation-local-gate.v1','rows':out,'layer':'solubility',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
        r.write_new(execution,status)
    print(r.read(execution));print(log.read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        data=compute()
        if a.mode=='worker':r.write_new(OUTPUT,data)
        else:assert r.read(OUTPUT)==data
        for row in data['rows']:print(row['prime'],row['status'],[(x['index'],x['Fp_point_count']) for x in row.get('relations',[])],row.get('geometry_gate'),flush=True)
