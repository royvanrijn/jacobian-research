#!/usr/bin/env sage-python
"""Fixed old-control replay of the halving-or-cycle seed classifier.

Each invocation: <=1009 proof primes, <=8 halving steps per candidate;
run under timeout 25s. Every completed candidate is immutable/checkpointed.
No search, new parameter or V3 read. The 302 carrier member remains openly
retrospectively calibrated. Prime selection reads only generic sections.
"""
import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, EllipticCurve, PolynomialRing
from split_seed_descent import Classifier, build_frame, point_record, require

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_split_descent_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
REDUCED=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
CHART=ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
CONIC=ART/'det1092_orbit8044_rank18_base_change_v2.json'
FIRST=ART/'det1092_norm8_seed_cover_v2/equation-only-cover-01.json'
BLIND=ART/'det1092_blind_mw16_v1'
PARAMETERS={'dependent-conic':'-528/3635','rank21-seed':'1926/2699',
            'other-conic':'2953/1671','small-conic':'5193/35630'}
SOURCES=[Path(__file__),Path(__file__).with_name('split_seed_descent.py'),
         Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')]
R=PolynomialRing(QQ,'t')


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(d,t):return R(d['numerator'])(t)/R(d['denominator'])(t)
def save(p,d):
    text=json.dumps(d,indent=2,sort_keys=True)+'\n';p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():require(p.read_text()==text,'immutable checkpoint changed: '+str(p))
    else:
        with p.open('x') as f:f.write(text)


def freeze():
    inputs=[PARENT,REDUCED,CHART,CONIC,FIRST,BLIND/'evaluation.json']
    evaluation=read(BLIND/'evaluation.json')
    arms=[r['arm_id'] for r in evaluation['rows'] if r['success']]
    require(len(arms)==9,'fixed nine successful old generic reconstructions')
    inputs += [BLIND/'arms'/a/'certificate.json' for a in arms]
    save(OUT/'protocol.json',dict(schema='det1092-halving-cycle-v1',
        classification='bounded verified application of the halving-or-cycle theorem',
        parameters=PARAMETERS, curve302_parameter='0', blind_arms=arms,
        maximum_halving_steps=8, prime_cap=1009, per_invocation_wall_seconds=25,
        sources={str(p.relative_to(ROOT)):sha(p) for p in SOURCES},
        inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        candidate_rule='Both roots and their anti-traces for each old conic; both roots and anti-traces of the previously calibrated genus-one302 carrier; first already verified generic recovery from each of nine old arms, tested over its M16 core and full M17.',
        boundary='No new address, exceptional point coordinate input, later discovery, V3 landscape, rank catalogue, point search or production change. Carrier01 selection is retrospective and stays labelled. Generic-only proof footprint is frozen before candidate construction. Caps retain UNKNOWN.'))


def protocol():
    p=read(OUT/'protocol.json')
    for name,h in {**p['sources'],**p['inputs']}.items():require(sha(ROOT/name)==h,'source/input changed: '+name)
    return p


def generic(family):
    if family=='302':
        p=read(PARENT);t=QQ(0)
        old=EllipticCurve(QQ,[val(a,t) for a in p['a_invariants']])
        E=EllipticCurve(QQ,[-old.c4()/48,-old.c6()/864])
        def transport(x,y):return E([x+old.b2()/12,y+(old.a1()*x+old.a3())/2])
        basis=[transport(val(x,t),val(y,t)) for x,y in p['basis_weierstrass_coordinates']]
        return E,basis,old,transport
    p=read(REDUCED);s=QQ(PARAMETERS[family]);b=s.denominator()
    E=EllipticCurve(QQ,[val(a,s)*b**k for a,k in zip(p['a_invariants'],[2,4,6,8,12])])
    basis=[E([val(x,s)*b**4,val(y,s)*b**6]) for x,y in p['basis_weierstrass_coordinates']]
    return E,basis,None,None


def candidates(family,E,basis,old,transport):
    rows=[]
    if family=='302':
        cover=read(FIRST);F=R(cover['quartic_coefficients'])(0)
        require(F and F.is_square(),'calibrated carrier must split at302')
        maps=cover['maps'];points=[]
        for sign in [-1,1]:
            W=sign*F.sqrt()
            x=val(maps['x0'],0)+W*val(maps['x1'],0)
            y=val(maps['y0'],0)+W*val(maps['y1'],0)
            old([x,y]);points.append(transport(x,y))
        trace=points[0]+points[1]
        require(trace==basis[14]-basis[15],'norm8 inherited trace')
        origin='Equation-only split of an openly first-witness-calibrated genus-one member'
    else:
        cover=read(CONIC);s=QQ(PARAMETERS[family]);chart=read(CHART)
        aa,bb,cc,dd=map(QQ,chart['parameter_matrix']);tau=(aa*s+bb)/(cc*s+dd)
        old=EllipticCurve(QQ,[val(a,tau) for a in read(PARENT)['a_invariants']])
        rc,rb,ra=[val(v,tau) for v in cover['lift']['residual_coefficients']]
        f0,f1,f2=[R(v)(tau) for v in cover['lift']['line_coefficients']]
        disc=rb*rb-4*ra*rc
        require(ra and f2 and disc and disc.is_square(),'old conic split and map patch')
        h=cc*s.numerator()+dd*s.denominator();w=QQ(chart['weierstrass_u']);points=[]
        for sign in [-1,1]:
            x=(-rb+sign*disc.sqrt())/(2*ra);y=-(f0+f1*x)/f2
            old([x,y])
            points.append(E([(x+old.b2()/12)*h**4/w**2,
                             (y+(old.a1()*x+old.a3())/2)*h**6/w**3]))
        trace=sum((int(n)*P for n,P in zip(cover['lift']['trace_word'],basis)),E(0))
        require(points[0]+points[1]==trace,'conic trace identity')
        origin='Generic orbit8044 residual quadratic and line; no seed packet read'
    for i,P in enumerate(points):
        rows.append((f'branch-{i}',P,None,origin))
        rows.append((f'anti-trace-{i}',2*P-trace,None,'Derived anti-trace; its mod2 class is inherited by construction'))
    if family=='302':
        evaluation=read(BLIND/'evaluation.json')
        for arm in [r for r in evaluation['rows'] if r['success']]:
            name=arm['arm_id'];cert=read(BLIND/'arms'/name/'certificate.json')
            j=arm['first_success_candidate_id'];raw=cert['candidates'][j]
            x,y=[val(d,0) for d in raw['coordinates']]
            P=transport(x,y)
            omitted=arm['omitted_basis_index_one_based']-1
            rows.append((name+'-core16',P,omitted,'Old generic RR reconstruction; diagnostic over its retained core'))
            rows.append((name+'-full17',P,None,'Same old generic RR point; diagnostic over full generic group'))
    return rows


def run(family):
    p=protocol();folder=OUT/family
    E,basis,old,transport=generic(family)
    frame_path=folder/'generic-frame.json'
    if not frame_path.exists():save(frame_path,build_frame(E,basis,p['prime_cap']))
    frame=read(frame_path)
    require(frame['curve']==list(map(str,E.a_invariants())) and frame['basis']==[point_record(P) for P in basis], 'generic frame specialization')
    engine=Classifier(frame)
    # Candidate construction is strictly after the generic frame checkpoint.
    rows=candidates(family,E,basis,old,transport)
    for label,P,omitted,origin in rows:
        destination=folder/(label+'.json')
        if destination.exists():
            saved=read(destination)
            require(saved['original']==point_record(P) and saved['frame_sha256']==sha(frame_path), 'resume input mismatch')
            print('RETAINED',family,label,saved['status'],flush=True)
            continue
        if omitted is not None:
            core=dict(frame);core['basis']=[v for i,v in enumerate(frame['basis']) if i!=omitted]
            core['inherited_rank']=16
            core['finite_rows']=[[v for i,v in enumerate(row) if i!=omitted] for row in frame['finite_rows']]
            core['records']=[{**r,'codes':[v for i,v in enumerate(r['codes']) if i!=omitted]} for r in frame['records']]
            current=Classifier(core)
        else:current=engine
        result=current.classify(P,p['maximum_halving_steps'])
        result.update(family=family,label=label,origin=origin,omitted_basis_index=omitted,
            frame_sha256=sha(frame_path),protocol_sha256=sha(OUT/'protocol.json'))
        save(destination,result)
        print(family,label,result['status'],result.get('reason'),result['steps'],flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('action',choices=['freeze','run'])
    a.add_argument('--family',choices=['302',*PARAMETERS]);args=a.parse_args()
    if args.action=='freeze':freeze();print('FROZEN',flush=True)
    else:
        require(args.family is not None,'family required');run(args.family)
