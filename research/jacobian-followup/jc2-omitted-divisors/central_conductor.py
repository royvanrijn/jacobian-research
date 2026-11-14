"""Exact conductor and degree-at-most-seven target representatives in both charts."""
from pathlib import Path
import json,sys
import sympy as s
from flint import fmpq,fmpq_poly,fmpq_mpoly_ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F

U,V,W,z=s.symbols('U V W z')
PT=fmpq_mpoly_ctx.get(('t','a'))
UV=fmpq_mpoly_ctx.get(('U','V','a'))
tt,aa=PT.gens()
uu,vv,ab=UV.gens()


def convert(expr,variables,ctx):
    return ctx.from_dict({mon:fmpq(str(co)) for mon,co in s.Poly(expr,*variables).terms()})


def evaluate(poly,p,q):
    rows=poly.to_dict()
    pp=[p**i for i in range(max((i for i,j,k in rows),default=0)+1)]
    qq=[q**i for i in range(max((j for i,j,k in rows),default=0)+1)]
    return sum((co*pp[i]*qq[j]*aa**k for (i,j,k),co in rows.items()),PT.constant(0))


def leading_t(poly):
    rows=poly.to_dict()
    degree=max(i for i,j in rows)
    last=max(j for (i,j),co in rows.items() if i==degree)
    return int(degree),fmpq_poly([rows.get((degree,j),fmpq(0)) for j in range(last+1)])


def embed(poly,ctx):
    size=ctx.nvars()
    return ctx.from_dict({(0,)*(size-1)+(i,):co for i,co in enumerate(poly.coeffs()) if co})


def encode(poly):
    return [[[int(i) for i in mon],str(co)] for mon,co in sorted(poly.to_dict().items())]


def decode(rows,ctx):
    assert len(rows)==len({tuple(mon) for mon,co in rows})
    return ctx.from_dict({tuple(mon):fmpq(co) for mon,co in rows})


def specialize(poly,value):
    rows=poly.to_dict()
    degree=max(i for i,j in rows)
    return fmpq_poly([sum((co*fmpq(value)**j for (i,j),co in rows.items() if i==n),fmpq(0))
                      for n in range(degree+1)])


def build_models(check):
    source=json.loads((HERE.parent/'jc2-r7-elliptic-contact/models.json').read_text())
    result=[]
    for row in source['families']:
        name=row['family'];parameter=F.d if name=='L' else z
        if name=='L':
            p0=s.expand(F.p.subs(F.line));q0=s.expand(F.q.subs(F.line))
        else:
            def scaled(poly,degree):
                raw=s.Poly(s.expand(poly.subs(F.laurent).subs(F.t,F.t/F.x)*F.x**degree),F.x,F.t)
                assert all(mon[0]%3==0 for mon,co in raw.terms())
                return s.expand(sum(co*z**(mon[0]//3)*F.t**mon[1] for mon,co in raw.terms()))
            p0,q0=scaled(F.p,6),scaled(F.q,10)
        p,q=convert(p0,(F.t,parameter),PT),convert(q0,(F.t,parameter),PT)
        Gexpr=s.sympify(row['key_G'],locals={'U':U,'V':V,str(parameter):parameter})
        key=s.sympify(row['key_relation'],locals={'U':U,'V':V,'W':W,str(parameter):parameter})
        ell=convert(s.Poly(key,W).nth(1),(U,V,parameter),UV)
        L=convert(s.Poly(key,V).nth(1),(U,V,parameter),UV)
        H=convert(key.subs({W:0,V:0}),(U,V,parameter),UV)
        G=convert(Gexpr,(U,V,parameter),UV)
        g=evaluate(G,p,q)
        equation=G*G+ell*G+H+L*vv
        check(name+' exact key/implicit restriction',not(g*g+evaluate(ell,p,q)*g+evaluate(H,p,q)+evaluate(L,p,q)*q))
        check(name+' implicit degree10 monic V6',max(i+j for i,j,k in equation.to_dict())==10
              and {(i,j,k):co for (i,j,k),co in equation.to_dict().items() if j==6}=={(0,6,0):fmpq(1)})
        derivative_V=(2*g+evaluate(ell,p,q))*evaluate(G.derivative(1),p,q)+evaluate(L,p,q)
        conductor=derivative_V/p.derivative(0)
        check(name+' exact gradient quotient and degree24',conductor*p.derivative(0)==derivative_V and leading_t(conductor)[0]==24)
        derivative_U=((2*g+evaluate(ell,p,q))*evaluate(G.derivative(0),p,q)
            +evaluate(ell.derivative(0),p,q)*g+evaluate(H.derivative(0),p,q)+evaluate(L.derivative(0),p,q)*q)
        check(name+' second gradient/conductor identity',derivative_U==-q.derivative(0)*conductor)
        expected=fmpq_poly([fmpq(-1,6561)]) if name=='L' else fmpq_poly([0]*6+[fmpq(2,19683)])
        check(name+' conductor leading coefficient unit on entire domain',leading_t(conductor)[1]==expected)
        residual=conductor
        adjoint=UV.constant(0)
        basis=sorted([(6*i+10*j+9*k,i,j,k) for i in range(5) for j in range(3) for k in range(2) if 6*i+10*j+9*k<=24],reverse=True)
        check(name+' thirteen distinct filtered basis poles',len(basis)==13 and len({row[0] for row in basis})==13)
        for degree,i,j,k in basis:
            image=p**i*q**j*g**k
            if not residual or leading_t(residual)[0]<degree:
                continue
            assert leading_t(residual)[0]==degree
            coefficient,remainder=divmod(leading_t(residual)[1],leading_t(image)[1])
            assert not remainder
            residual-=embed(coefficient,PT)*image
            adjoint+=embed(coefficient,UV)*uu**i*vv**j*G**k
        check(name+' exact conductor adjoint restriction',not residual and evaluate(adjoint,p,q)==conductor)
        check(name+' adjoint total target degree at most7',max(i+j for i,j,k in adjoint.to_dict())<=7)
        result.append(dict(family=name,parameter=str(parameter),domain='all complex d' if name=='L' else 'z!=0',
            normalization_p=encode(p),normalization_q=encode(q),conductor=encode(conductor),
            adjoint_representative=encode(adjoint),conductor_degree=24,
            adjoint_degree=int(max(i+j for i,j,k in adjoint.to_dict())),
            conductor_leading_parameter_polynomial=[str(co) for co in expected.coeffs()],
            parameter_denominators_constant=True,extra_parameter_loci_discarded=False))
    return result
