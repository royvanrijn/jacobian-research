"""Controls for the coupled three-cell comparison and scalar radial obstruction.

The full-rho simplex inequality is not asserted. Run with the existing
Python environment containing python-flint; all retained interval inputs
are fixed below. The analytic proofs are DG5--DG8 in the direct note.
"""
from fractions import Fraction
import json
import math
from flint import arb, ctx


def enclose(value):
    return [math.nextafter(float(value.lower()),-math.inf),
            math.nextafter(float(value.upper()),math.inf)]


def quadrant(prod,c):
    delta=prod-c*c
    return 1/delta-c*(c/prod.sqrt()).acos()/(delta*delta.sqrt())


def first_mode(prod,rho,d,j1):
    delta=prod-rho*rho
    angle=(rho/prod.sqrt()).acos()
    return j1*rho*d/2*((prod+2*rho*rho)*angle/(delta*delta*delta.sqrt())
                      -3*rho/(delta*delta))


def certificate(panels=8192,precision=100,rho_value=Fraction(1,2)):
    ctx.prec=precision
    assert rho_value in (Fraction(1,2),Fraction(3,5),Fraction(7,10))
    rho=arb(rho_value.numerator)/rho_value.denominator
    d_exact=1-rho_value*rho_value
    d=arb(d_exact.numerator)/d_exact.denominator;pi=arb.pi()
    j1=(9-pi*arb(3).sqrt())/(4*pi*pi)
    exponents=(0,1,16)
    coefficients=(Fraction(-21,34),Fraction(1),Fraction(2))
    assert sum(c/Fraction(1+a) for c,a in zip(coefficients,exponents))==0
    coeff=[arb(c.numerator)/c.denominator for c in coefficients]
    nterms={};dterms={}
    for a,c in zip(exponents,coeff):
        for b,z in zip(exponents,coeff):
            # Products are rational here; identify equal ones without using
            # an approximate floating-point key.
            kn=(Fraction(1)+a*d_exact)*(Fraction(1)+b*d_exact)
            kd=Fraction(1)+(a+b)*d_exact
            nterms[kn]=nterms.get(kn,arb(0))+c*z
            dterms[kd]=dterms.get(kd,arb(0))+c*z
    ndata=[(arb(p.numerator)/p.denominator,c) for p,c in nterms.items()]
    ddata=[(arb(p.numerator)/p.denominator,c) for p,c in dterms.items()]
    nvalue=dvalue=arb(0)
    width=arb(2)/(3*panels)
    for i in range(panels):
        lo=arb(2*i)/(3*panels);hi=arb(2*(i+1))/(3*panels)
        v=lo.union(hi)
        c=-rho*(pi*v).cos()
        fn=sum((weight*d*quadrant(p,c) for p,weight in ndata),arb(0))
        fd=sum((weight*d*quadrant(p,c) for p,weight in ddata),arb(0))
        nvalue+=width*(2+3*v)*fn/2
        dvalue+=width*(1-3*v)*fd/2
    nb=sum((weight*first_mode(p,rho,d,j1) for p,weight in ndata),arb(0))
    db=sum((weight*first_mode(p,rho,d,j1) for p,weight in ddata),arb(0))
    nvalue-=4*nb
    dvalue+=db
    residual=nvalue-2*dvalue
    correction=(db-nb)/(4*j1)
    corrected=residual-correction
    sharp_corrected=residual-2*correction
    assert dvalue>0 and residual>0,(str(nvalue),str(dvalue),str(residual))
    assert correction>0
    if rho_value==Fraction(1,2):
        assert corrected<0,str(corrected)
    else:
        assert corrected>0,str(corrected)
    if rho_value in (Fraction(1,2),Fraction(3,5)):
        assert sharp_corrected<0,str(sharp_corrected)
    else:
        assert sharp_corrected>0,str(sharp_corrected)
    return {'status':'pass','rho':str(rho_value),'exponents':list(exponents),
            'coefficients':[str(c) for c in coefficients],
            'exact_mean':'0','panels':panels,'arb_precision_bits':precision,
            'angular_rule':'interval Riemann enclosures over every subinterval of [0,2/3]',
            'numerator_enclosure':enclose(nvalue),'denominator_enclosure':enclose(dvalue),
            'numerator_minus_twice_denominator':enclose(residual),
            'rayleigh_quotient_enclosure':enclose(nvalue/dvalue),
            'closure_correction_enclosure':enclose(correction),
            'corrected_quadratic_coefficient':enclose(corrected),
            'sharp_first_moment_coefficient':enclose(sharp_corrected),
            'exact_gram_controls':gram_controls(),
            'scope':'The scalar radial bound fails at this selected correlation. The elementary correction repairs the profile at 1/2; the sharp local first-moment correction repairs it also at 3/5 but fails at 7/10. Positive upper-bound coefficients do not prove improving Gaussian partitions. Finite Gram controls support the written algebraic lemma, not an unrestricted proof by sampling.'}


def dot(x,y):
    return sum((a*b for a,b in zip(x,y)),Fraction(0))


def gram_bound(a,b,u,v):
    assert all(sum(row[j] for row in u)==0 for j in (0,1))
    assert all(sum(row[j] for row in v)==0 for j in (0,1))
    assert all(dot(row,row)<=length*length for row,length in zip(u,a))
    assert all(dot(row,row)<=length*length for row,length in zip(v,b))
    r=[y/x for x,y in zip(a,b)]
    eta=[(r[i]-r[(i+1)%3])*(r[i]-r[(i+2)%3]) for i in range(3)]
    numerator=sum(x*x*y for x,y in zip(a,eta))
    defect=[x*x-dot(row,row) for x,row in zip(a,u)]
    actual=tuple(sum(r[i]*u[i][j] for i in range(3)) for j in (0,1))
    assert numerator==dot(actual,actual)+sum(x*y for x,y in zip(defect,eta))
    limit=max([Fraction(0)]+[x/y for x,y in zip(eta,r)])
    bound=numerator/(2*(sum(r)+limit))
    delta=sum(x*y for x,y in zip(a,b))-sum(dot(x,y) for x,y in zip(u,v))
    assert 0<=bound<=delta
    return bound,delta


def gram_controls():
    a=tuple(map(Fraction,(3,4,5)))
    b=tuple(map(Fraction,(5,5,6)))
    u=tuple(tuple(map(Fraction,row)) for row in ((3,0),(0,4),(-3,-4)))
    v=((Fraction(5),Fraction(0)),(Fraction(-7,5),Fraction(24,5)),
       (Fraction(-18,5),Fraction(-24,5)))
    cases=[gram_bound(a,b,u,v),gram_bound(b,a,v,u),gram_bound(a,a,u,u)]
    contracted=tuple((x*Fraction(3,4),y*Fraction(3,4)) for x,y in u)
    rotated=tuple(((3*x-4*y)/5,(4*x+3*y)/5) for x,y in v)
    cases.append(gram_bound(a,b,contracted,rotated))
    bad=(tuple(x+Fraction(1,10) for x in u[0]),*u[1:])
    try:
        gram_bound(a,b,bad,v)
    except AssertionError:
        pass
    else:
        raise AssertionError('nonclosing triples must be rejected')
    return {'closing_rational_cases':len(cases),'nonclosing_control_rejected':True,
            'bound_and_defect_pairs':[[str(x),str(y)] for x,y in cases]}


if __name__=='__main__':
    result=certificate()
    result['three_fifths_limit_test']=certificate(rho_value=Fraction(3,5))
    result['seven_tenths_limit_test']=certificate(rho_value=Fraction(7,10))
    print(json.dumps(result,indent=2))
