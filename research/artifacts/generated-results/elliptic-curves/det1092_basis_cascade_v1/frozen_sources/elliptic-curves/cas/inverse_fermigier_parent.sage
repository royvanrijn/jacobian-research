#!/usr/bin/env sage-python
"""Bounded rational (u,v) inverse Fermigier recognition, unbounded rational T.

No parameter or point from the known 245 parent is used in selection. The
control is evaluated after recognition. This is not a complete search of
Fermigier moduli, a basis recognizer, or an original-provenance certificate.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd, lcm
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from icarm_curve245_mestre import fermigier_roots, PUBLIC_MODEL
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS

PRIMES = [101,103,107,109,127,131,137,139,149,151,157,163]
WORK = ROOT/'artifacts/local/elliptic-curves/inverse-fermigier-parent'
RESULTS = ROOT/'artifacts/generated-results/elliptic-curves'


def save(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    tmp.replace(path)


def digest(path): return sha256(path.read_bytes()).hexdigest()


def configurations(height):
    rationals=sorted({QQ(a)/b for a in range(-height,height+1) for b in range(1,height+1)})
    families={}; degenerate=0
    for u in rationals:
        for v in rationals:
            rr=list(map(QQ,fermigier_roots(Fraction(str(u)),Fraction(str(v)))))
            if len(set(rr))<6:
                degenerate+=1; continue
            den=lcm([r.denominator() for r in rr]); ints=sorted(ZZ(den*r) for r in rr)
            shift=ints[0]; gg=gcd([x-shift for x in ints])
            roots=tuple((x-shift)//gg for x in ints)
            reflected=tuple(roots[-1]-x for x in reversed(roots))
            sign=-1 if reflected<roots else 1; canonical=min(roots,reflected)
            families.setdefault(canonical,[]).append({'u':str(u),'v':str(v),
                'native_to_canonical_T_scale':str(sign*QQ(den)/gg)})
    return families, {'rational_parameter_count':len(rationals),
                      'ordered_pairs':len(rationals)**2,'degenerate_pairs':degenerate,
                      'distinct_affine_root_configurations':len(families)}


def build_invariants(roots):
    Rt=PolynomialRing(QQ,'T'); T=Rt.gen(); RX=PolynomialRing(Rt,'X'); X=RX.gen()
    product=RX(1)
    for a in roots: product *= (X-a-T)*(X-a+T)
    g=X**6
    for k in range(5,-1,-1): g += (product[6+k]-(g*g)[6+k])/2*X**k
    rem=g*g-product
    assert rem.degree()<=4, 'Mestre condition failed'
    coeff=[]
    for i in range(5):
        q,r=rem[i].quo_rem(T*T); assert r==0; coeff.append(q)
    e,d,c,b,a=coeff
    I=12*a*e-3*b*d+c*c
    J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    assert all(I[i]==J[i]==0 for i in range(1,max(I.degree(),J.degree())+1,2))
    Rz=PolynomialRing(QQ,'Z')
    iz=Rz([I[2*i] for i in range(I.degree()//2+1)])
    jz=Rz([J[2*i] for i in range(J.degree()//2+1)])
    if 4*iz**3-jz**2==0: return None
    numerator=6912*iz**3; denominator=4*iz**3-jz**2
    common=numerator.gcd(denominator)
    return iz,jz,numerator//common,denominator//common


def primitive(f):
    den=lcm([c.denominator() for c in f])
    coeff=[ZZ(den*c) for c in f]
    gg=gcd(coeff)
    return f.parent()([c/gg for c in coeff])


def recognize(roots, origins, target):
    inv=build_invariants(roots)
    row={'roots':list(map(str,roots)),'origins':origins}
    if inv is None:
        row['outcome']='IDENTICALLY_SINGULAR'; return row
    I,J,num,den=inv; targetj=target.j_invariant()
    comp=num*targetj.denominator()-den*targetj.numerator()
    if not comp:
        row['outcome']='UNRESOLVED_CONSTANT_J_MATCH'; return row
    comp=primitive(comp); maximum=max(num.degree(),den.degree())
    infinity_match=comp.degree()<maximum
    # A rational T reduces to an Fp point of P1, even if its denominator
    # is divisible by p. Finite Z=T^2 values are precisely the squares.
    # Never exclude on degree drop or an identically-zero reduction.
    for p in PRIMES:
        if infinity_match or not ZZ(comp[maximum])%p: continue
        coeff=[int(ZZ(c)%p) for c in comp]
        squares=sorted({i*i%p for i in range(p)})
        def evaluate(z):
            ans=0
            for c in reversed(coeff): ans=(ans*z+c)%p
            return ans
        if all(evaluate(z) for z in squares):
            row.update(outcome='EXCLUDED_MODULAR_PROJECTIVE_T',prime=p,
                       comparison_Z_coefficients_mod_p=coeff,
                       comparison_degree_Z=int(comp.degree()))
            return row
    factorization=comp.factor(); matches=[]; candidates=[]; unresolved=False
    for f,m in factorization:
        if f.degree()!=1: continue
        z=-f[0]/f[1]
        square=z>=0 and z.is_square()
        item={'Z':str(z),'multiplicity':int(m),'rational_T':bool(square)}
        candidates.append(item)
        if not square: continue
        for t in sorted({z.sqrt(),-z.sqrt()}):
            aa,bb=-27*I(z),-27*J(z)
            if 4*aa**3+27*bb**2==0:
                # A cancelled j-map factor can hide a removable nonminimal
                # fibre. Do not reject without constructing its local model.
                item['singular_raw_jacobian']=True; unresolved=True; continue
            E=EllipticCurve(QQ,[aa,bb]); assert E.j_invariant()==targetj
            isos=E.isomorphisms(target)
            candidate={'T_canonical':str(t),'Q_isomorphic':bool(isos),
                       'raw_short_model':list(map(str,E.a_invariants()))}
            if isos:
                candidate['to_target_u_r_s_t']=list(map(str,isos[0].tuple()))
                candidate['native_parameters']=[{'u':o['u'],'v':o['v'],
                    'T':str(t/QQ(o['native_to_canonical_T_scale']))} for o in origins]
                matches.append(candidate)
            item.setdefault('fibre_candidates',[]).append(candidate)
    row.update(outcome=('UNRESOLVED_INFINITY' if infinity_match else
                        'UNRESOLVED_LOCAL_MINIMALIZATION' if unresolved else
                        'Q_ISOMORPHISM_MATCH' if matches else 'EXCLUDED_EXACT'),
               factor_degrees_multiplicities=[[int(f.degree()),int(m)] for f,m in factorization],
               rational_Z_candidates=candidates,matches=matches,
               comparison_Z_coefficients=list(map(str,comp.list())))
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',choices=['245','302'],required=True)
    parser.add_argument('--height',type=int,default=8)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not 1<=args.height<=8: parser.error('frozen maximum height is eight')
    label=f'{args.target}_h{args.height}'
    output=RESULTS/f'curve302_inverse_fermigier_{label}_v1.json'
    checkpoint=WORK/f'{label}_checkpoint.json'
    dependencies=[Path(__file__),Path(__file__).with_name('icarm_curve245_mestre.py'),
                  Path(__file__).with_name('mestre_root_tuples.py'),
                  Path(__file__).with_name('icarm_curve302.py')]
    bindings={str(p.relative_to(ROOT)):digest(p) for p in dependencies}
    # Search control must pass before the target is admitted. The evaluator
    # (below) reads the known control parameters only after recognition.
    if args.target=='302':
        control=RESULTS/f'curve302_inverse_fermigier_245_h{args.height}_v1.json'
        cc=json.loads(control.read_text())
        assert cc['status']=='COMPLETE' and cc['control_recovered'] is True
        assert cc['input_sha256']==bindings
    protocol={'height_uv':args.height,'height_definition':'max(abs(numerator),denominator) in lowest terms',
              'T_search':'all rational T by exact factorization in Z=T^2',
              'prime_budget':PRIMES,'maximum_wall_seconds':900,'workers':1,
              'failure':'timeout or unresolved models stay UNKNOWN; preserve checkpoint',
              'target':args.target,'input_sha256':bindings}
    WORK.mkdir(parents=True,exist_ok=True)
    def timeout(*unused): raise TimeoutError('frozen 900-second budget exhausted')
    signal.signal(signal.SIGALRM,timeout); signal.alarm(900)
    before=time.monotonic()
    families,counts=configurations(args.height)
    print('CONFIGURATIONS',counts,flush=True)
    target=EllipticCurve(QQ,list(map(QQ,PUBLIC_MODEL if args.target=='245' else GENERAL_WEIERSTRASS_COEFFICIENTS)))
    rows=[]; previous=None
    if args.check: previous=json.loads(output.read_text())
    elif args.resume and checkpoint.exists():
        prior=json.loads(checkpoint.read_text()); assert prior['protocol']==protocol
        rows=prior['rows']
    elif checkpoint.exists():
        raise FileExistsError('use --resume for the existing checkpoint')
    complete=True
    try:
        for index,(roots,origins) in enumerate(sorted(families.items())):
            if index<len(rows):
                assert rows[index]['roots']==list(map(str,roots)) and rows[index]['origins']==origins
                continue
            row=recognize(roots,origins,target); rows.append(row)
            if not args.check:
                save(checkpoint,{'protocol':protocol,'counts':counts,'rows':rows})
            if row.get('matches') or (index+1)%25==0:
                print('PROGRESS',index+1,len(families),row['outcome'],round(time.monotonic()-before,2),flush=True)
    except TimeoutError:
        complete=False
    signal.alarm(0)
    hits=[{'roots':r['roots'],'matches':r['matches']} for r in rows if r.get('matches')]
    # Retrospective control evaluation, isolated from candidate generation.
    control_recovered=None
    if args.target=='245':
        control_recovered=any(QQ(p['u'])==QQ(3)/2 and QQ(p['v'])==2 and
            abs(QQ(p['T']))==QQ(5801)/160 for h in hits for a in h['matches'] for p in a['native_parameters'])
    result={'schema':'curve302.inverse-fermigier-recognition.v1',
            'status':'COMPLETE' if complete else 'INCOMPLETE_TIMEOUT',
            'input_sha256':bindings,'protocol':protocol,'software':{'sage':SAGE_VERSION},
            'counts':counts,'rows_completed':len(rows),'rows':rows,'hits':hits,
            'outcome_counts':dict(Counter(r['outcome'] for r in rows)),
            'control_recovered':control_recovered,
            'boundary':'Bounded rational u,v inverse construction search with exact elimination of rational T. Affine-root duplicates removed only; no claim of complete fibration deduplication. A Q-isomorphic fibre is not a full generic MW basis or original provenance.'}
    if args.target=='302': result['calibration_sha256']=digest(control)
    if args.check:
        assert complete and result==previous,'replay mismatch'
    else: save(output,result)
    print('RESULT',result['status'],result['outcome_counts'],'control',control_recovered,flush=True)


if __name__=='__main__': main()
