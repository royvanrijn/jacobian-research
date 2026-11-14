"""Uniform seven-choice normalization interface and literal unit replays."""
from pathlib import Path
import sys,json,gzip,hashlib,time,importlib.util
from fractions import Fraction
from collections import defaultdict
import sympy as S
from sympy.polys.rings import ring

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G
spec=importlib.util.spec_from_file_location('quintic_implicit',ROOT/'archive/non-elliptic/scripts/verify_f2_affine_target_k1_implicit_conductor.py')
I=importlib.util.module_from_spec(spec);spec.loader.exec_module(I)
A,M,s=S.symbols('A M s');CHECKS=[]
def check(name,condition):
    assert condition,name
    CHECKS.append(name)
def zero(name,value):check(name,S.expand(value)==0)
def terms(poly):return [[list(m),str(co)] for m,co in poly.terms()]
def stored(rows):return {tuple(m):Fraction(co) for m,co in rows if Fraction(co)}
def add_product(result,left,right):
    for a,c in left.items():
        for b,d in right.items():result[tuple(x+y for x,y in zip(a,b))]+=c*d
def read_unit(label):return json.loads(gzip.decompress((HERE/f'unit-{label}.json.gz').read_bytes()))
def replay_unit(label,generators):
    data=read_unit(label)
    check(label+' exact source generators',data['generators']==generators)
    result=defaultdict(Fraction)
    for gen,co in zip(data['generators'],data['cofactors']):add_product(result,stored(gen),stored(co))
    result={m:co for m,co in result.items() if co}
    rhs={(0,0,m[0]):Fraction(co) for m,co in data['rhs'] if Fraction(co)}
    check(label+' literal Fraction cofactor identity',result==rhs)
    return stored(data['rhs'])

def reconstruct(family,p,q,param):
    equation,_,_=G.approximate_quintic(p,q,G.F.t)
    b,c,e,f,n=S.symbols('b c e f n')
    pb=s**3+A*s+M;qb=s**5+b*s**4+c*s**3+e*s**2+f*s+n
    value=S.Poly(S.expand(equation.subs({G.U:pb,G.V:qb})),s)
    recovered={}
    for degree,var in [(14,b),(13,c),(12,e),(11,f),(10,n)]:
        pol=S.Poly(S.expand(value.nth(degree).subs(recovered)),var)
        check(family+' triangular coefficient '+str(degree),pol.degree()==1 and pol.nth(1)==3)
        recovered[var]=S.factor(-pol.nth(0)/3)
    qb=S.expand(qb.subs(recovered))
    residual={k:S.expand(value.nth(k).subs(recovered)) for k in range(10)}
    # The known resultant formula directly constructs the irreducible image.
    implicit=-I.expected_implicit_quintic().subs({I.a:A,I.b:recovered[b],I.c:recovered[e],
        I.d:recovered[f]-recovered[c]*A,I.P:G.U-M,
        I.Q:G.V-recovered[c]*(G.U-M)-recovered[n]},simultaneous=True)
    implicit=S.Poly(S.expand(implicit),G.U,G.V)
    difference=S.Poly(S.expand(equation-implicit.as_expr()),G.U,G.V)
    support={(3,0),(1,1),(2,0),(0,1),(1,0),(0,0)}
    check(family+' exact implicit difference support',set(difference.monoms())<=support)
    zero(family+' U3 residual coefficient',difference.coeff_monomial(G.U**3)-residual[9])
    zero(family+' UV residual coefficient',difference.coeff_monomial(G.U*G.V)-residual[8])
    zero(family+' redundant seventh coefficient',residual[7]-3*A*residual[9]-recovered[b]*residual[8])
    # Monic weighted Groebner basis over QQ[parameter]; no parameter inversion.
    R,aa,mm=ring('A,M',S.QQ.poly_ring(param),order=lambda mon:(2*mon[0]+3*mon[1],mon[0]))
    ff=R.from_expr(S.Rational(27,5)*residual[9]);gg=R.from_expr(S.Rational(9,5)*residual[8])
    hh=-(mm*ff-aa*gg)/9;basis=[ff,gg,hh]
    check(family+' three constant leading units',[(v.LM,v.LC) for v in basis]==[((3,0),R.domain.one),((2,1),R.domain.one),((0,3),R.domain.one)])
    for i in range(3):
        for j in range(i+1,3):
            lm=tuple(max(x,y) for x,y in zip(basis[i].LM,basis[j].LM))
            di=tuple(x-y for x,y in zip(lm,basis[i].LM));dj=tuple(x-y for x,y in zip(lm,basis[j].LM))
            quotients,remainder=(basis[i].mul_monom(di)-basis[j].mul_monom(dj)).div(basis)
            check(family+f' Buchberger pair {i}{j}',not remainder)
    standard=[(i,j) for i in range(3) for j in range(3) if not(i>=2 and j>=1)]
    check(family+' seven free standard monomials',len(standard)==7)
    # One reduced exact fibre establishes that seven distinct curves occur on
    # a nonempty open; every fibre still has scheme length seven.
    witness=S.groebner([v.as_expr().subs(param,1) for v in [ff,gg]],M,A,order='lex')
    univariate=next(S.Poly(v,A) for v in witness.polys if not v.as_expr().has(M))
    linear=next(v.as_expr() for v in witness.polys if v.as_expr().has(M))
    check(family+' witness elimination degree7',univariate.degree()==7)
    check(family+' witness seven simple choices',S.gcd(univariate,univariate.diff()).degree()==0)
    check(family+' witness M uniquely recovered',S.Poly(linear,M).degree()==1 and S.Poly(linear,M).LC()==1)
    generators=[terms(S.Poly(residual[k],A,M,param).clear_denoms()[1].primitive()[1]) for k in range(9,3,-1)]
    rhs=replay_unit(family,generators)
    if family=='K':rhs2=replay_unit('K2',generators)
    else:check('L nonzero constant unit RHS',set(rhs)=={(0,)} and rhs[(0,)]!=0)
    row=dict(family=family,parameter=str(param),normalization_p=str(pb),normalization_q=str(qb),
        necessary_sufficient_equations=[str(residual[9]),str(residual[8])],
        groebner_basis=[str(v.as_expr()) for v in basis],standard_monomials=[list(v) for v in standard],
        high_contact_implicit='G_C minus the U2,V,U,1 coefficients of the displayed exact implicit difference',
        implicit_difference={str(mon):str(difference.coeff_monomial(G.U**mon[0]*G.V**mon[1])) for mon in sorted(support)},
        witness_univariate=str(univariate.as_expr()),witness_linear=str(linear),
        high_contact_count_with_multiplicity=7,possible_high_contacts=[38,40])
    return row,(rhs,rhs2) if family=='K' else None

def main():
    start=time.monotonic()
    line,_=reconstruct('L',G.F.p.subs(G.F.line),G.F.q.subs(G.F.line),G.F.d)
    laurent,rhs=reconstruct('K',G.K.P,G.K.Q,G.K.z)
    bezout=json.loads((HERE/'K-rhs-bezout.json').read_text());result=defaultdict(Fraction)
    add_product(result,stored(bezout['a']),rhs[0]);add_product(result,stored(bezout['b']),rhs[1])
    check('K literal denominator Bezout equals z7',{m:c for m,c in result.items() if c}=={(7,):Fraction(1)})
    result=dict(schema='r7-high-contact-polynomial-quintics-v1',families=[line,laurent],
        maximal_contact41='IMPOSSIBLE for every coefficient point in both charts',
        finite_intersection_minimum=10,
        conditional_scope='Additional image B must have polynomial normalization degrees(3,5) in the same coordinates. Seven choices count only infinity contact at least38, with multiplicity; no high contact is forced for actual B.',
        assurance='Authored normalization/ideal interpretation with exact algebra; no formalization or independent external review; no Keller realization.')
    (HERE/'normalizations.json').write_text(json.dumps(result,indent=2)+'\n')
    names=['verify.py','PROOF.md','README.md','REVIEW.md','unit-L.json.gz','unit-K.json.gz','unit-K2.json.gz','K-rhs-bezout.json','discovery.json']
    inputs=[HERE/name for name in names]
    inputs += [G.HERE/'geometry.py',HERE.parent/'jc2-r7-image-family/family.py',HERE.parent/'jc2-r7-laurent-collisions/laurent.py',ROOT/'archive/non-elliptic/scripts/verify_f2_affine_target_k1_implicit_conductor.py',ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(passed=True,check_count=len(CHECKS),checks=CHECKS,elapsed_seconds=round(time.monotonic()-start,3),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        output_sha256={'normalizations.json':hashlib.sha256((HERE/'normalizations.json').read_bytes()).hexdigest()},assurance=result['assurance'])
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} QUINTIC NORMALIZATION CHECKS ({time.monotonic()-start:.3f}s)',flush=True)

if __name__=='__main__':main()
