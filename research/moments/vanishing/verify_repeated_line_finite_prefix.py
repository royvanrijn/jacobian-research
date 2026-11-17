#!/usr/bin/env python3
"""Exact degree bounds for the written repeated-line99-moment corollary."""
from pathlib import Path
import hashlib,json
import sympy as S


def connection(U,D,u,t):
    d=S.Poly(U,u).degree();a=d*S.expand(U).coeff(u,d)
    H=S.expand(D*D*a**d)
    source=S.Poly(S.expand(S.diff(U,t)*D-U*S.diff(D,t)-S.diff(D,t)*D/2),u)
    rows=[];counts=[]
    for j in range(d-1):
        # Work with the numerator, keeping the common D^2 outside.
        poly=S.Poly(source.as_expr()*u**j,u)
        denominator=S.Integer(1);count=0
        while poly.degree()>=d-1:
            n=poly.degree();lead=poly.LC()
            replacement=-sum(S.expand(S.diff(U,u)).coeff(u,k)*u**(n-d+1+k) for k in range(d-1))
            if n-d+1:replacement-=(n-d+1)*D*u**(n-d)
            poly=S.Poly(S.expand(a*(poly.as_expr()-lead*u**n)+lead*replacement),u)
            denominator*=a;count+=1
        assert count<=d
        padded=S.Poly(S.expand(poly.as_expr()*a**(d-count)),u)
        rows.append([padded.coeff_monomial(u**k) for k in range(d-1)])
        counts.append(count)
    N=S.Matrix(d-1,d-1,lambda i,j:rows[j][i])
    return N,H,counts


def main():
    p=Path(__file__).resolve().parent
    u,t=S.symbols('u t')
    degree_rows=[]
    for d in [3,4]:
        # Every numerator coefficient has the allowed maximal t-degree;
        # the leading coefficient has full cubic degree and no t=0 root.
        U=sum((1+(j+2)*t+(2*j+1)*t*t+(j*j+1)*t**3)*u**j for j in range(d+1))
        for D in [1+2*t+3*t*t,t*(2*t-1)]:
            N,H,counts=connection(U,D,u,t);B=4+3*d
            hd=S.Poly(H,t).degree()
            nd=max(S.Poly(e,t).degree() for e in N)
            assert hd<=B and nd<=B
            h=S.eye(d-1)[:,0]
            iterated=[]
            for j in range(1,d):
                h=S.expand(H)*h.diff(t)-(j-1)*S.diff(H,t)*h+N*h
                h=h.applyfunc(S.expand)
                degrees=[S.Poly(e,t).degree() for e in h]
                maximum=max(degrees);assert maximum<=j*B
                iterated.append({'iterate':j,'max_numerator_degree':int(maximum),'bound':j*B})
            C=B*(d-1)*d//2
            degree_rows.append({'phase_degree':d,'denominator':str(D),'denominator_at_zero':str(D.subs(t,0)),
                                'B':B,'H_degree':int(hd),'N_max_degree':int(nd),
                                'reduction_counts':counts,'iterated_columns':iterated,
                                'scalar_cofactor_degree_bound':C,'pure_exclusion_prefix':C+d-1})
    assert [r['pure_exclusion_prefix'] for r in degree_rows]==[41,41,99,99]
    # Both eliminations have the universal t-degree3/2 numerator/denominator
    # bound. Verify the nontrivial isotropic completion with arbitrary
    # polynomial coefficient placeholders; no conditional t=0 average.
    Z,T=S.symbols('Z T');lam,A=S.symbols('lam A',nonzero=True)
    H0=S.Symbol('h20')*Z**2+S.Symbol('h10')*Z+S.Symbol('h01')*T+S.Symbol('h00')
    G=S.Symbol('g30')*Z**3+S.Symbol('g21')*Z**2*T+S.Symbol('g20')*Z**2+S.Symbol('g11')*Z*T+S.Symbol('g02')*T*T+S.Symbol('g10')*Z+S.Symbol('g01')*T+S.Symbol('g00')
    psi=S.expand(t*(G-H0**2/(4*lam))+Z*H0/(2*lam)-T*T/2-Z*Z/(4*lam*t))
    cT=S.expand(psi).coeff(T,2);lT=S.expand(psi).coeff(T,1);kT=S.expand(psi).coeff(T,0)
    assert S.Poly(cT,t).degree()<=1 and S.Poly(lT,t).degree()<=1
    numerator=S.expand(t*cT*kT-t*lT*lT/4)
    assert S.Poly(numerator,t).degree()<=3
    assert S.Poly(t*cT,t).degree()<=2
    # Nonisotropic resolvent: D quadratic, adjugate linear.
    aa,ab,bb=S.symbols('aa ab bb');matrix=S.Matrix([[aa,ab],[ab,bb]])
    resolvent=S.eye(2)-2*t*matrix;det=S.expand(resolvent.det());adj=resolvent.adjugate()
    assert S.Poly(det,t).degree()==2
    assert all(S.Poly(e,t).degree()<=1 for e in adj)
    assert S.Poly(t*det,t).degree()==3
    assert all(S.Poly(t*t*e,t).degree()<=3 for e in adj)
    # The actual surviving frames have q2=0 or d*Z^2, hence no WZT/WT^2.
    W=S.symbols('W');alpha,beta,delta=S.symbols('alpha beta delta')
    assert S.diff(T*T*(alpha*T+beta*Z),W)==0
    assert S.diff(Z*Z*(alpha*Z+beta*T+delta*W),W)==delta*Z*Z
    # The only c*d interior has G3=alpha*Z^3 after the null rotation.
    a,b,c,d=S.symbols('a b c d',nonzero=True)
    g20,g11,g02,g10,g01,g00=S.symbols('g20 g11 g02 g10 g01 g00')
    h=(Z-t*(a+b*Z+d*Z*Z))/(c*t)
    mixed=S.expand(t*(alpha*Z**3+g20*Z*Z+g11*Z*h+g02*h*h+g10*Z+g01*h+g00)-h*h/2)
    assert S.Poly(mixed,Z).degree()==4
    assert S.factor(S.expand(mixed).coeff(Z,4)-d*d*(t*g02-S.Rational(1,2))/(c*c))==0
    assert S.Poly(S.expand(t*t*mixed),t).degree()<=3
    assert (3*4+3)*3*4//2+3==93
    prior=json.loads((p/'cubic_gmc_closure_verification.json').read_text())
    frame=json.loads((p/'cubic_frame_prefix_verification.json').read_text())
    assert prior['completed'] and frame['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_repeated_line_checker_sha256':prior['checker_sha256'],
            'prior_frame_prefix_checker_sha256':frame['checker_sha256'],
            'elimination_t_degree_bounds':{'numerator':3,'denominator':2},
            'connection_models':degree_rows,'pure_exclusion_bound':99,'all_repeated_line_pure_gate':99,
            'adapted_frame_family':'W*(a+bZ+cT+dZ^2)+G, with G3=alpha*Z^3 whenever d!=0',
            'adapted_frame_prefixes':{'c_zero':35,'d_zero':92,'c_d_nonzero_restricted_G3':93},
            'global_arbitrary_cubic_pure_gate':None,
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Universal degree algebra for both quadratic eliminations; exact quotient reductions, denominators and iterated-column degree bounds on four rational degree-complete models including singular t=0 denominators; adapted leading-form algebra and restricted mixed quartic/numerator-degree identities;35/92/93/99 arithmetic. The q0!=0 analytic exclusions, universal cofactor proof, frame choices and all-order99 corollary remain written mathematics. No moment99 expansion, optimal pure gate, global arbitrary-cubic gate, independent review or formal verification is asserted.'}
    (p/'repeated_line_finite_prefix_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','pure_exclusion_bound','all_repeated_line_pure_gate','global_arbitrary_cubic_pure_gate','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
