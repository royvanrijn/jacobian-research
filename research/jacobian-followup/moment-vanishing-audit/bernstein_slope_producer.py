"""Finite postprocessing of VALID actual Bernstein certificate data.

This module does not construct or verify the underlying D-module certificates.
See EFFECTIVE_BERNSTEIN_SLOPE_SELECTION_AND_CUTOFFS.md for the exact contract.
All algebra in the postprocessor is rational and exact. Its optional original-
power evaluator accepts rational phase coefficients and a complete finite jet.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import math
import sympy as S

X,Y,Z=S.symbols('X Y Z')

def polynomial(terms, variables):
    out=0
    for exponents, coefficient in terms:
        if len(exponents)!=len(variables) or any(int(e)!=e or e<0 for e in exponents):
            raise ValueError('Invalid nonnegative exponent vector')
        out+=S.Rational(coefficient)*S.prod(v**int(e) for v,e in zip(variables,exponents))
    return S.Poly(out,*variables,domain=S.QQ)

def terms(poly):
    return [[list(m),str(c)] for m,c in poly.terms() if c]

def root_bound(poly):
    """All complex roots have modulus <= this rational Cauchy bound."""
    if poly.is_zero:raise ValueError('Zero polynomial has no finite root bound')
    if poly.degree()==0:return S.Rational(0)
    co=poly.all_coeffs()
    return 1+max(abs(c/co[0]) for c in co[1:])

def integer_roots(poly):
    if poly.is_zero:raise ValueError('Specialized polynomial is identically zero')
    return sorted(int(z) for z in poly.ground_roots() if z.q==1 and z>=0)

def make_plan(certificate_summary):
    """Assume C1/C2 in the proof; construct a finite exact query plan."""
    p=polynomial(certificate_summary['polynomial_terms'],(X,Y))
    if p.is_zero:raise ValueError('The order-constraint polynomial must be nonzero')
    k0=certificate_summary['validity_start'];J=certificate_summary['step_bound']
    if type(k0)!=int or k0<1 or type(J)!=int or J<0:
        raise ValueError('Positive integer start and nonnegative integer step bound required')
    normalized=dict(polynomial_terms=terms(p),validity_start=k0,step_bound=J)
    slopes={}
    for f,_ in S.factor_list(p.as_expr(),X,Y)[1]:
        f=S.Poly(f,X,Y,domain=S.QQ)
        if f.total_degree()!=1:continue
        ay=f.coeff_monomial(Y)
        if ay==0:continue
        rho=-f.coeff_monomial(X)/ay
        if not (0<rho<=J):continue
        a,b=int(rho.p),int(rho.q)
        offset=-b*f.coeff_monomial(1)/ay
        slopes.setdefault((a,b),[]).append(offset)
    strips=[];K=k0
    for (a,b),offsets in sorted(slopes.items(),key=lambda z:S.Rational(*z[0])):
        q=S.Poly(p.as_expr().subs(Y,(a*X+Z)/b),X,Z,domain=S.QQ)
        degree_x=q.degree(X)
        f=S.Poly(q.as_expr().coeff(X,degree_x),Z,domain=S.QQ)
        R=root_bound(f);B=int(S.ceiling(R))+1;L=a+b*J;A=B+L+1
        mu=abs(f.LC())
        error=sum(abs(c)*A**j for (i,j),c in q.terms() if i<degree_x)
        start=max(1,int(S.floor(error/mu))+1)
        K=max(K,start)
        strips.append(dict(a=a,b=b,rho=str(S.Rational(a,b)),normal_radius=B,
            normal_step_bound=L,collar_outer_radius=A,collar_validity_start=start,
            transformed_polynomial_terms=terms(q),leading_x_coefficient_terms=terms(f),
            leading_coefficient_root_bound=str(R),collar_lower_bound=str(mu),error_numerator_bound=str(error),
            integral_line_offsets=sorted(set(int(c) for c in offsets if c.q==1))))
    for i,s in enumerate(strips):
        for t in strips[i+1:]:
            gap=abs(S.Rational(s['a'],s['b'])-S.Rational(t['a'],t['b']))
            widths=S.Rational(s['normal_radius'],s['b'])+S.Rational(t['normal_radius'],t['b'])
            K=max(K,int(S.floor(widths/gap))+1)
    # This one nonzero coefficient prevents any identically zero Y specialization.
    coefficient=S.Poly(p.as_expr(),Y).LC()
    K=max(K,int(S.floor(root_bound(S.Poly(coefficient,X,domain=S.QQ))))+1)
    specialized=S.Poly(p.as_expr().subs(X,K),Y,domain=S.QQ)
    queries=integer_roots(specialized)
    return dict(schema_version=1,kind='finite postprocessing of supplied actual Bernstein certificates',
        certificate_summary=normalized,
        assumptions_not_verified_by_this_program=[
            'For every k>=validity_start,each finite actual tau_k satisfies p(k,tau_k)=0.',
            'For every k>=validity_start,nonzero v_k forces nonzero v_(k+1),and finite orders have increments in[0,J].',
            'The input is an actual positive-degree polynomial phase and a nonzero zero-constant formal symbol over a characteristic-zero field.'],
        query_operator_power=K,phase_power_candidates=queries,
        symbol_jet_order_required_per_phase_degree=max(queries,default=0),strips=strips,
        specialized_polynomial_terms=terms(specialized),
        original_power_queries='For each listed j,determine exactly whether T^K(P^j) is nonzero.',
        universal_certificate_generator_implemented=False,formal_verification=False)

def resolve(plan,evaluations,*,all_positive_pure=False):
    """Resolve from exact original-power nonzero predicates, not guessed roots."""
    if type(all_positive_pure)!=bool:raise ValueError('Purity premise must be an explicit boolean')
    if plan!=make_plan(plan['certificate_summary']):raise ValueError('Plan failed exact recomputation')
    if len(evaluations)!=len(plan['phase_power_candidates']):raise ValueError('Missing or extra evaluations')
    outcomes={}
    for e in evaluations:
        j=e['phase_power']
        if j in outcomes or type(e['nonzero'])!=bool:raise ValueError('Invalid evaluation record')
        outcomes[j]=e['nonzero']
    if set(outcomes)!=set(plan['phase_power_candidates']):raise ValueError('Evaluation powers differ from plan')
    K=plan['query_operator_power'];nonzero=[j for j in sorted(outcomes) if outcomes[j]]
    common=dict(conditional_on_supplied_actual_certificates=True,
                all_positive_purity_assumed=bool(all_positive_pure),query_operator_power=K,
                independent_replay=False,formal_verification=False,external_review=False)
    if not nonzero:
        return dict(common,alternative='eventually infinite first-survival orders',nilpotence_index_bound=K,
                    all_marker_eventual_safety=True,affine_cutoff=dict(C=1,E=K-1),
                    pure_linear_cutoff=K if all_positive_pure else None)
    tau=min(nonzero)
    matching=[s for s in plan['strips'] if abs(s['b']*tau-s['a']*K)<s['normal_radius']]
    if len(matching)!=1:raise ValueError('Original-power evaluations contradict the supplied certificate contract')
    s=matching[0];a,b=s['a'],s['b'];rho=S.Rational(a,b)
    out=dict(common,alternative='all first-survival orders finite',first_survival_at_query=tau,
             selected_rational_slope=str(rho),normal_radius=s['normal_radius'],
             all_marker_eventual_safety=bool(a>b))
    if a>b:
        H=max(K-1,math.ceil(Fraction(s['normal_radius'],a)))
        C=math.ceil(Fraction(a,a-b));E=math.ceil(Fraction(a*H,a-b))
        out.update(global_height_intercept=H,affine_cutoff=dict(C=C,E=E),
                   pure_linear_cutoff=C+E if all_positive_pure else None)
    elif all_positive_pure:
        if a<b:raise ValueError('The supplied all-positive purity premise contradicts slope<1')
        candidates=[c for c in s['integral_line_offsets'] if c>0]
        if not candidates:raise ValueError('No positive critical offset is compatible with supplied purity')
        out.update(critical_bad_marker_candidates_P_powers=candidates,
                   marker_guarantee='At least ONE listed fixed marker P^d has infinitely many nonzero outputs;the program does not select which one.')
    return out

def _multiply(f,g):
    out={}
    for a,c in f.items():
        for b,d in g.items():
            e=tuple(x+y for x,y in zip(a,b));out[e]=out.get(e,Fraction(0))+c*d
    return {e:c for e,c in out.items() if c}

def _power(f,j,n):
    out={(0,)*n:Fraction(1)}
    while j:
        if j&1:out=_multiply(out,f)
        j//=2
        if j:f=_multiply(f,f)
    return out

def _operator(f,symbol):
    out={}
    for exponents,coefficient in f.items():
        for gamma,c in symbol.items():
            if any(a<b for a,b in zip(exponents,gamma)):continue
            v=coefficient*c
            for a,b in zip(exponents,gamma):
                v*=math.prod(range(a-b+1,a+1))
            e=tuple(a-b for a,b in zip(exponents,gamma))
            out[e]=out.get(e,Fraction(0))+v
    return {e:c for e,c in out.items() if c}

def evaluate_rational(plan,phase_terms,symbol_jet_terms,complete_through):
    """Evaluate the finite plan on Q coefficients and an asserted complete jet.

    Omitted jet coefficients through complete_through are asserted to be zero.
    The complete infinite formal symbol need not be supplied or convergent.
    """
    if plan!=make_plan(plan['certificate_summary']):raise ValueError('Plan failed exact recomputation')
    if not phase_terms:raise ValueError('Nonzero phase required')
    n=len(phase_terms[0][0])
    def sparse(ts):
        out={}
        for e,c in ts:
            if len(e)!=n or any(type(a)!=int or a<0 for a in e):raise ValueError('Invalid exponent')
            e=tuple(e);out[e]=out.get(e,Fraction(0))+Fraction(c)
        return {e:c for e,c in out.items() if c}
    phase=sparse(phase_terms);symbol=sparse(symbol_jet_terms)
    D=max((sum(e) for e in phase),default=0)
    if D<=0:raise ValueError('Positive phase degree required')
    needed=D*plan['symbol_jet_order_required_per_phase_degree']
    if type(complete_through)!=int or complete_through<needed:raise ValueError('Insufficient complete symbol jet')
    if symbol.get((0,)*n,0):raise ValueError('The formal symbol must have zero constant term')
    if any(sum(e)>complete_through for e in symbol):raise ValueError('Jet contains terms beyond its stated order')
    K=plan['query_operator_power'];outputs=[]
    for j in plan['phase_power_candidates']:
        f=_power(phase,j,n)
        for _ in range(K):
            f=_operator(f,symbol)
            if not f:break
        serialized=[[list(e),str(c)] for e,c in sorted(f.items())]
        digest=hashlib.sha256(json.dumps(serialized,separators=(',',':')).encode()).hexdigest()
        outputs.append(dict(phase_power=j,nonzero=bool(f),output_terms=len(f),output_sha256=digest))
    return dict(phase_degree=D,complete_jet_required=needed,complete_jet_supplied=complete_through,evaluations=outputs)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['plan','evaluate','resolve'])
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();data=json.loads(args.input.read_text())
    if args.mode=='plan':result=make_plan(data)
    elif args.mode=='evaluate':result=evaluate_rational(data['plan'],data['phase_terms'],data['symbol_jet_terms'],data['complete_through'])
    else:result=resolve(data['plan'],data['evaluations'],all_positive_pure=data.get('all_positive_pure',False))
    args.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
