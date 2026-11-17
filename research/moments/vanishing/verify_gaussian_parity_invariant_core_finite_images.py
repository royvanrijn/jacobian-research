#!/usr/bin/env python3
"""Bounded exact kernel and original Gaussian controls for parity invariant cores.

The checker records finite identities and primary-mirror pins.  It does
not certify the all-m continuation, monodromy argument or EIT theorem.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

def kernel_controls():
    eta, rho, v = sp.symbols("eta rho v")
    R = sp.symbols("R", positive=True)


    def zero(value):
        assert sp.simplify(sp.expand(value)) == 0, value


    order = 8
    em = []
    for m in range(11):
        s = sp.Rational(3, 2)-m
        em.append(-sum((-eta)**k*R**(s+k)/(sp.factorial(k)*(s+k))
                       for k in range(order+1)))
    for m in range(1, 11):
        s = sp.Rational(3, 2)-m
        exponential = sum((-eta*R)**k/sp.factorial(k) for k in range(order+1))
        expression = s*em[m]-eta*em[m-1]+R**s*exponential
        zero(sp.series(expression, eta, 0, order+1).removeO())
        zero(s*sp.gamma(s)-sp.gamma(s+1))

    # I0 is an elementary exact primitive. The Im integration recurrence is
    # checked before using either a special function or a continued branch.
    q = 1-rho+rho*v**2
    zero(sp.diff(v/sp.sqrt(q), v)-(1-rho)*q**sp.Rational(-3, 2))
    for m in range(2, 11):
        power = sp.Rational(2*m-3, 2)
        zero(sp.diff(v*q**power, v)
             -(2*m-2)*q**power+(2*m-3)*(1-rho)*q**(power-1))
        cm = sp.prod(sp.Rational(2*k-3, 2*k-2) for k in range(2, m+1))
        zero(cm-sp.gamma(m-sp.Rational(1, 2))/(sp.sqrt(sp.pi)*sp.gamma(m)))
        zero(-sp.gamma(sp.Rational(3, 2)-m)*cm/sp.sqrt(sp.pi)
             -(-1)**m/sp.factorial(m-1))
    zero(-sp.gamma(sp.Rational(1, 2))/sp.sqrt(sp.pi)+1)

    # Root-filter trace without floating roots of unity.
    for degree in range(1, 9):
        companion = sp.zeros(degree)
        for j in range(degree):
            companion[(j+1) % degree, j] = 1
        zero(sp.trace((sp.eye(degree)-rho*companion).inv())
             -degree/(1-rho**degree))

    # The original beta-normalization is rational, and the angular substitution
    # has the correct (normalized) Jacobian and radial exponent.
    for n in range(7):
        for degree in range(1, 6):
            zero(sp.gamma(sp.Rational(3, 2))/sp.gamma(degree*n+sp.Rational(3, 2))
                 -1/sp.rf(sp.Rational(3, 2), degree*n))
    r, x = sp.symbols("r x", positive=True)
    s, u = r*(1-x**2), r*x**2
    jacobian = sp.det(sp.Matrix([[sp.diff(s, r), sp.diff(s, x)],
                                 [sp.diff(u, r), sp.diff(u, x)]]))
    zero(jacobian/sp.sqrt(u)-2*sp.sqrt(r))
    return {
        "status":"passed","coefficient_field":"QQ",
        "incomplete_Laplace_m_indices":[0,10],
        "entire_remainder_coefficient_orders":[0,8],
        "angular_recurrence_m_indices":[2,10],
        "root_filter_trace_degrees":[1,8],
        "normalization_n_indices":[0,6],
        "normalization_radial_degrees":[1,5],
        "scope":"Finite exact identities only; global continuation is the written all-m convergence argument.",
    }


def literal_controls():
    S,T,Z,W,Y,x,r,t=sp.symbols('S T Z W Y x r t')
    beta=sp.Rational(3,2)
    def gamma(value):
        return sp.expand(sum(c*sp.factorial(i)*sp.rf(sp.Rational(1,2),j)
                             for (i,j),c in sp.Poly(sp.expand(value),S,T).terms()))
    def original(value):
        return sp.expand(sum(c*sp.factorial(i)*sp.factorial2(k-1)
                             for (i,j,k),c in sp.Poly(sp.expand(value),Z,W,Y).terms()
                             if i==j and k%2==0))
    def angular(value):
        expanded=value.subs({S:r*(1-x*x),T:r*x*x},simultaneous=True)
        return sp.expand(sum(c*sp.rf(beta,i)/(j+1)
                             for (i,j),c in sp.Poly(sp.expand(expanded),r,x).terms()))
    even=[]
    fixtures=[S+2*T+1,S*S+S*T+T*T+S-2*T+1,(S+T)**2+3*S+T+1,S**3+T**2+1]
    for f in fixtures:
        phase=f.subs({S:Z*W,T:Y*Y/2},simultaneous=True)
        D=sp.Poly(f,S,T).total_degree()
        values=[]
        for n in range(5):
            raw=gamma(f**n)
            assert sp.expand(raw-original(phase**n))==0
            assert sp.expand(raw-angular(f**n))==0
            values.append(str(raw))
        top=sum(co*S**i*T**j for (i,j),co in sp.Poly(f,S,T).terms() if i+j==D)
        A=sp.expand(top.subs({S:1-x*x,T:x*x},simultaneous=True))
        assert A!=0 and sp.expand(A.subs(x,-x)-A)==0
        assert sp.degree(A,x)<=2*D
        even.append({'f':str(f),'radial_degree':D,'angular_top':str(A),
                     'angular_degree':int(sp.degree(A,x)),'raw_moments_0_to_4':values})
    odd=[]
    for b in [S+T+1,S-T+T*T]:
        phase=Y*b.subs({S:Z*W,T:Y*Y/2},simultaneous=True)
        Q=2*T*b*b
        m=sp.Poly(b,S,T).total_degree()
        assert sp.Poly(Q,S,T).total_degree()==2*m+1
        values=[]
        for n in range(4):
            raw=gamma(Q**n)
            assert sp.expand(raw-original(phase**(2*n)))==0
            assert original(phase**(2*n+1))==0
            assert sp.rf(beta,(2*m+1)*n)>=sp.factorial(2*n)
            values.append(str(raw))
        top=sum(co*S**i*T**j for (i,j),co in sp.Poly(b,S,T).terms() if i+j==m)
        A=sp.expand(2*x*x*top.subs({S:1-x*x,T:x*x},simultaneous=True)**2)
        assert A!=0 and A.subs(x,0)==0 and sp.degree(A,x)>0
        odd.append({'b':str(b),'Q':str(sp.expand(Q)),'Q_degree':2*m+1,
                    'nonflat_angular_top':str(A),'Q_raw_moments_0_to_3':values})
    def truncate(value,order):
        return sp.expand(sum(co*t**int(power[0])
            for power,co in sp.Poly(sp.expand(value),t).terms() if int(power[0])<=order))
    def exp_series(value,order):
        answer=1
        power=1
        for j in range(1,order+1):
            power=truncate(power*value,order)
            answer+=power/sp.factorial(j)
        return truncate(answer,order)
    remainder=Z**5*W**3*Y**4+Z**2
    finite=[]
    order=6
    for c,b in [(0,0),(3,0),(3,2)]:
        phase=c+b*Y+remainder
        powers=[sp.Integer(1)]
        for _ in range(order):
            powers.append(sp.expand(powers[-1]*phase))
        factor=exp_series(c*t+b*b*t*t/2,order)
        controls=[]
        for marker in [sp.Integer(1),W**2,Y*W**2,W**4]:
            D=sp.Poly(marker,Z,W,Y).total_degree()
            normalized=sp.expand(sum(t**n*original(marker.subs(Y,Y+b*t)
                             *remainder.subs(Y,Y+b*t)**n)/sp.factorial(n)
                                    for n in range(D//2+1)))
            actual=sum(original(marker*powers[j])*t**j/sp.factorial(j) for j in range(order+1))
            assert sp.expand(actual-truncate(factor*normalized,order))==0
            bound=D//2 if b==0 else D+5*(D//2)
            assert normalized==0 or sp.degree(normalized,t)<=bound
            if marker==1:
                assert normalized==1
            if marker==W**2:
                assert sp.expand(normalized-(362*t+720*b*b*t**3+120*b**4*t**5))==0
            controls.append({'marker':str(marker),'normalized_exact_formula':str(normalized),
                             'sufficient_degree_bound':bound})
        finite.append({'constant_core':c,'linear_Y_core':b,'R':'Z^5 W^3 Y^4+Z^2',
                       'gamma':2,'original_moment_orders':[0,order],'markers':controls})
    # Nonfinite core plus one-sided remainder: retain actual scalar identity.
    f=fixtures[1]
    phase=f.subs({S:Z*W,T:Y*Y/2},simultaneous=True)+remainder
    for n in range(4):
        assert sp.expand(original(phase**n)-gamma(f**n))==0
    return {'even_Gamma_original_angular_fixtures':even,'even_fixture_orders':[0,4],
            'odd_squared_original_fixtures':odd,'odd_even_moment_orders':[0,6],'odd_odd_moment_orders':[1,7],
            'finite_original_branches':finite,'nonfinite_core_vacuum_invariance_orders':[0,3],
            'scope':'Bounded exact identities only; the all-degree continuation and classification are written proofs.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,help="Write the exact JSON receipt to this path")
    args=parser.parse_args()
    checks={"kernel":kernel_controls(),"literal_original":literal_controls()}
    here=Path(__file__).resolve()
    repository=here.parents[3]
    canonical=[here]+[here.with_name(name) for name in [
        "GAUSSIAN_PARITY_INVARIANT_CORE_FINITE_IMAGES.md",
        "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
        "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
        "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
        "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
    ]]
    base=repository/"research/factorial/image-rigidity"
    primary=["sources/long-fc2.tex","sources/long-eit-pinned.tex",
             "sources/pakovich-muzychuk-07104085v2.html","sources/beukers-math0405549.html"]
    references=[base/name for name in ["LONG_EIT_AND_FC2_AUDIT.tex","source_ledger.json",
                                      "long_source_audit.json"]+primary]
    ledger=json.loads((base/"source_ledger.json").read_text(encoding="utf-8"))
    for name in primary:
        entry=next(item for item in ledger["cached_primary_sources"] if item["file"]==name)
        assert hashlib.sha256((base/name).read_bytes()).hexdigest()==entry["sha256"],name
    checks["primary_mirror_ledger_consistency"]={"status":"passed","files":primary}
    paths=canonical+references
    receipt={
        "status":"passed","checked_on":"2026-11-17",
        "environment":{"python":platform.python_version(),"sympy":sp.__version__},
        "sha256":{p.relative_to(repository).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in paths},
        "checks":checks,
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "unrestricted_result_from_finite_prefix":False,
        "scope":"Bounded exact kernel, Gamma/Gaussian and finite-branch controls. The unrestricted result is the written proof with pinned unrefereed primary manuscripts, retained audits and published transitive inputs; no external or formal assurance upgrade.",
    }
    encoded=json.dumps(receipt,indent=2)+"\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8")
    print(encoded,end="")

if __name__=="__main__":
    main()
