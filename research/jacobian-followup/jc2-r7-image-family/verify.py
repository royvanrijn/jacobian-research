"""Replay the complete two-family coefficient classification, without Groebner calls."""
import hashlib
import json
import time
from fractions import Fraction as QF
import sympy as S
import family as F

CHECKS=[]
def check(name, value):
    assert value, name
    CHECKS.append(name)

def zero(name, value):
    check(name, S.cancel(value)==0)

def decode(rows):
    result={}
    for m,c0 in rows:
        assert len(m)==6 and all(type(k) is int and k>=0 for k in m)
        m=tuple(m)
        assert m not in result
        result[m]=QF(c0)
    return result

def sparse(expr):
    return {m:QF(str(c0)) for m,c0 in S.Poly(expr,*F.VARIABLES).terms() if c0}

def add_product(out, left, right):
    for m,c0 in left.items():
        for n,c1 in right.items():
            k=tuple(a+b for a,b in zip(m,n))
            out[k]=out.get(k,QF(0))+c0*c1

def main():
    start=time.monotonic()
    a,b,c,d,B,D,E,t,x,X,Y,u=F.a,F.b,F.c,F.d,F.B,F.D,F.E,F.t,F.x,F.X,F.Y,F.u
    data=json.loads((F.HERE/'certificates.json').read_text())
    check('explicit coefficient symbols',data['variables']==list(map(str,F.VARIABLES)))
    check('membership dimensions',len(data['membership'])==5 and all(len(r)==3 for r in data['membership']))
    for i,f in enumerate(F.generators):check(f'original generator{i}',decode(data['generators'][i])==sparse(f))
    for j,g in enumerate(F.targets):
        check(f'graph target{j}',decode(data['targets'][j])==sparse(g))
        result={}
        for i,f in enumerate(data['generators']):add_product(result,decode(f),decode(data['membership'][i][j]))
        result={m:c0 for m,c0 in result.items() if c0}
        check(f'independent Fraction graph membership{j}',result==sparse(g))
        result[(0,)*6]=result.get((0,)*6,QF(0))+1
        check(f'mutated membership rejected{j}',result!=sparse(g))
    zero('graph c equation',F.g0.subs(F.graph))
    zero('graph D equation',F.g2.subs(F.graph))
    zero('graph E equation',F.rho[1].as_expr().subs(b,1).subs(F.graph))
    zero('complete split relation',F.g1.subs(F.graph)-9*(a-2*d)*F.Q)
    zero('rho5 is redundant on the graph',
         (F.rho[5].as_expr()+2*a*F.rho[3].as_expr()/3).subs(b,1).subs(F.graph))

    # Exceptional a=2d chart; no division by d or u before it is justified.
    f=16*d**3+d*u**2+8
    g=416*d**6-16*d**5*u+26*d**4*u**2-d**3*u**3+216*d**3-7*d**2*u+4
    ell=d**2*u+8*d**3+4
    sub={a:2*d,B:(u-7*d)/6}
    zero('linear chart rho3',F.reduced[0].subs(sub)-d*u*f)
    zero('linear chart rho7',F.reduced[1].subs(sub)-u*g)
    zero('linear chart polynomial reduction',g-(-d**2*u+26*d**3)*f-ell)
    zero('linear chart exceptional parameter identity',
         d**3*f-(d**2*u-8*d**3-4)*ell-8*(2*d**3+1)*(5*d**3+2))
    zero('d0 forces u0',g.subs(d,0)-4)
    zero('first exceptional root forces u0',S.rem(ell,2*d**3+1,d)-d**2*u)
    zero('second exceptional root forces u2d',S.rem(ell-d**2*(u-2*d),5*d**3+2,d))
    zero('off-conic exceptional root impossible',F.Q.subs({a:2*d,B:-S.Rational(5,6)*d}))

    # The Q-conic and its exceptional homogeneous direction.
    linear={a:2*d+(X-Y)/6,B:(X+7*(2*d+(X-Y)/6)-21*d)/6}
    zero('conic homogeneous relation',F.Q.subs(linear)-(F.N-6*F.Z*d)/6)
    zero('Z0 on conic forces X0',F.N.subs(Y,5*X/3)-10*X**2)
    zero('conic parameterization',F.Q.subs(F.conic))
    K=5*X**2-6*X*Y+5*Y**2
    zero('conic rho3 factor',F.reduced[0].subs(F.conic)-X*K*F.R/F.Z**4)
    P7=S.cancel(2*F.Z**6*F.reduced[1].subs(F.conic)/F.R)
    check('conic rho7 quotient is polynomial',S.fraction(P7)[1]==1)
    zero('X0 cofactor is nonzero off Z0',P7.subs(X,0)+54*Y**4)
    r=S.Symbol('r');kr=5*r**2-6*r+5
    cofactor=S.Poly(S.expand(P7.subs(Y,r*X)/X**4),X)
    zero('quadratic slope removes high term',S.rem(cofactor.nth(3),kr,r))
    zero('quadratic slope residual linear term',
         S.rem(cofactor.nth(0),kr,r)+S.Rational(256,125)*(139*r-95))
    check('quadratic slope cannot equal forced value',kr.subs(r,S.Rational(95,139))==S.Rational(62500,19321))

    # Inverse and forward Laurent chart, plus the exceptional main component.
    XY={X:-(3*x**3+2)/x**2,Y:S.Rational(5,2)*x*(3*x**3+2)}
    zero('Laurent relation',F.R.subs(XY))
    zero('Laurent inverse',((X*Y)/F.Z).subs(XY)-x)
    zero('Laurent a reconstruction',F.conic[a].subs(XY)-F.laurent[a])
    zero('Laurent d reconstruction',F.conic[d].subs(XY)-F.laurent[d])
    zero('Laurent B reconstruction',F.conic[B].subs(XY)-F.laurent[B])
    zero('line a2d',F.line[a]-2*d)
    zero('line B relation',6*F.line[B]+7*d)
    intersection=S.Poly(S.together((6*B+7*d).subs(F.laurent)).as_numer_denom()[0],x)
    zero('family intersection support',intersection.as_expr()+(3*x**3+2)*(35*x**3+26))
    # The second factor of this B equation is removed by the simultaneous a equation.
    test=S.together((a-2*d).subs(F.laurent)).as_numer_denom()[0]
    gcd=S.Poly(S.gcd(intersection.as_expr(),test),x).monic().as_expr()
    zero('three common coefficient points',gcd-x**3-S.Rational(2,3))

    for name,sub,parameter in [('line',F.line,d),('laurent',F.laurent,x)]:
        for k in [1,3,5,7,9]:zero(f'{name} gap{k}',F.rho[k].as_expr().subs(sub))
        gap=F.rho[11].as_expr().subs(sub)
        zero(f'{name} gap exactly21',gap-(S.Rational(1,19683) if name=='line' else -S.Rational(2,59049)/x**3))
        pp=S.factor(F.p.subs(sub));qq=S.factor(F.q.subs(sub))
        dp=S.together(S.diff(pp,t)).as_numer_denom()[0]
        dq=S.together(S.diff(qq,t)).as_numer_denom()[0]
        resultant=S.factor(S.resultant(dp,dq,t))
        expected=(-1011316948992*(3*d**3+2) if name=='line' else
                  5263136739342894229381323813814272*x**67*(3*x**3+2)*(15*x**3-4))
        zero(f'{name} exact nonimmersion resultant',resultant-expected)
        check(f'{name} exact degrees',(S.degree(pp,t),S.degree(qq,t))==(6,10))
        zero(f'{name} monic p',S.Poly(pp,t).LC()-1)
        zero(f'{name} monic q',S.Poly(qq,t).LC()-1)
    # Finite normalization of the Laurent coefficient curve.
    zero('Laurent parameter is integral over B',5*x**3+3*F.laurent[B]*x**2+1)
    zero('Laurent inverse parameter is integral over B',x**-3+3*F.laurent[B]/x+5)
    zero('sharp witness retained',F.p.subs(F.line).subs(d,0)-t**6-t**3)
    beta=S.Symbol('beta')
    zero('quadratic involution centering',S.Poly(F.p.subs(t,-t+beta),t).nth(5)+6*beta)
    zero('centered involution contradicts b1',
         S.Poly(F.p.subs(b,1).subs(t,-t)-F.p.subs(b,1),t).nth(3)+2)
    families=dict(schema='r7-image-family-v1',scope='Complete complex coefficient points with b=1; no coefficient-scheme equality or global Keller realization.',
                  families={name:dict(parameter=str(parameter),domain=domain,
                                      coefficients={str(v):str(S.factor(sub[v])) for v in [a,b,c,d,B,D,E]},
                                      rho11=str(S.factor(F.rho[11].as_expr().subs(sub))),
                                      nonimmersion_equation=bad)
                            for name,sub,parameter,domain,bad in [
                                ('line',F.line,d,'A1','3*d^3+2=0'),
                                ('laurent',F.laurent,x,'Gm: x!=0','(3*x^3+2)*(15*x^3-4)=0')]},
                  unknowns=['second ramified image equation','optional unramified image equation',
                            'actual image intersections','global complement relations','finite normal cover and Keller source realization'])
    (F.HERE/'families.json').write_text(json.dumps(families,indent=2)+'\n')
    inputs=[F.HERE/n for n in ['PROOF.md','README.md','REVIEW.md','family.py','verify.py','certificates.json']]
    inputs.append(F.ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py')
    receipt=dict(schema='r7-image-family-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),
                 input_sha256={str(p.relative_to(F.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                 output_sha256={'families.json':hashlib.sha256((F.HERE/'families.json').read_bytes()).hexdigest()},
                 coefficient_set='Union of a polynomial A1 family and a Laurent Gm family; set-theoretic classification, not equality of coefficient schemes.',
                 assurance='Exact graph memberships, rational identities and resultants; authored normalization and completeness proof, no formalization or independent external review.')
    (F.HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} TWO-FAMILY CHECKS ({receipt["elapsed_seconds"]}s)')

if __name__=='__main__':main()
