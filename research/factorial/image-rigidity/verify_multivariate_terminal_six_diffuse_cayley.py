#!/usr/bin/env python3
"""Exact shared terminal-six angular separation on .60 <= u <= .72, x5 <= 1/8."""

from fractions import Fraction as Q
from functools import lru_cache
from sys import argv

from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, KNOTS as EARLIER_KNOTS, ONE, bisect, imaginary_coordinate, interval,
    linear_form, log, one_fifth_power_complement, projections, rational,
    require, scalar, third_tail, verify_scalar_constants,
)


def skew_bounds(strength, reciprocal, cap, log_cap, sextic):
    # Sum x <= min((g^-2-1)/(j-1), (cap/s_cap)*strength).
    n = scalar(5 if sextic else 4)
    total_cap = strength * cap / log_cap
    reciprocal_bound = (ONE - reciprocal.square()) / (n*reciprocal.square())
    total_cap = I(0, min(total_cap.hi, reciprocal_bound.hi))
    outer = cap * total_cap
    link = cap**3 * total_cap
    if sextic:
        outer = scalar(5, 3)*outer
        link = scalar(6)*link
        inner = (scalar(4)*outer+link) * scalar(1, 5)
    else:
        link = scalar(3, 5)*link
        inner = (outer+link)*scalar(1, 2)
    require((outer+link).hi < scalar(3,50).lo and ((scalar(5,4) if sextic else scalar(2))*inner).hi < scalar(3,50).lo,
            "small hyperbolic arguments")
    log_defect = log_cap.square()*strength*scalar(1, 5 if sextic else 6)
    return total_cap, outer, inner, link, log_defect


def contributions_five(p, lam, cap, log_cap):
    u = ONE-p
    U = -log(u)
    _, alpha, chi, link, V = skew_bounds(U,u,cap,log_cap,False)
    R = scalar(3,8) + p*third_tail(scalar(3,2),p)
    v0 = ONE-scalar(3,2)*p+p.square()*R
    error = (v0 * V / (ONE-V)).positive_error()
    v = v0 + error
    d_over_p = scalar(1,2)-p*R-error/p
    a, b = scalar(501,1000)*u*alpha.square(), error+scalar(501,1000)*v*chi.square()
    f = (scalar(1001,1000)*chi).symmetric_error()
    epsilon = (u*(scalar(501,500)*link+scalar(101,100)*chi**3)).symmetric_error()
    nominal = ((1,scalar(5,8)*R),
               (3,scalar(5,4)*(ONE-p*R)),
               (5,ONE-scalar(5,4)*p+scalar(5,8)*p.square()*R))
    families = (
        (a, ((1,-scalar(15,16)/p.square()),(3,scalar(5,8)/p),(5,scalar(5,16)))),
        (b, ((1,scalar(5,8)/p.square()),(3,-scalar(5,4)/p),(5,scalar(5,8)))),
        (f, ((0,-scalar(5,8)*d_over_p/(lam**3)),
             (2,-scalar(5,4)*(u+v)/(lam**3)),
             (4,scalar(5,8)*(scalar(3)*u+v)/lam))),
        (epsilon, ((0,-scalar(5,16)/(lam**5)),
                   (2,-scalar(5,8)/(lam**3)),
                   (4,scalar(15,16)/lam))),
    )
    return nominal, families


def contributions_six(p, lam, rho, cap, log_cap):
    r = rho-p*(rho*(rho-ONE)*scalar(1,2)+p*third_tail(rho,p))
    q = r*p
    w = ONE-q
    require(w.lo > scalar(57,100).hi, "w lower bound for central defect")
    W = rho*(-log(ONE-p))
    total, beta, eta, link, Y = skew_bounds(W,w,cap,log_cap,True)
    J = scalar(6,5)*cap**4*total
    R_y, R_z = third_tail(scalar(8,5),q), third_tail(scalar(9,5),q)
    y0 = ONE-scalar(8,5)*q+scalar(12,25)*q.square()+R_y*q**3
    D = (y0*Y/(ONE-Y)).positive_error()
    y = y0 + D
    d_over_q = scalar(3,5)-scalar(12,25)*q-R_y*q.square()-D/q
    a, b = scalar(501,1000)*w*beta.square(), scalar(501,1000)*y*eta.square()
    E = I(-(scalar(71,50)*J).hi,(scalar(17,20)*D.square()).hi)
    f = (scalar(1001,1000)*eta).symmetric_error()
    epsilon = (w*(scalar(501,2000)*link+eta**3)).symmetric_error()
    nominal = ((0,r**3*(scalar(15)*R_y-scalar(10)*R_z)*scalar(1,32)),
               (2,r.square()*(scalar(9,20)+scalar(15,32)*(scalar(2)*R_z-R_y)*q)),
               (4,r*(scalar(3,2)-scalar(9,10)*q-scalar(15,32)*(R_y+scalar(2)*R_z)*q.square())),
               (6,ONE-scalar(3,2)*q+scalar(9,20)*q.square()+
                  (scalar(15)*R_y+scalar(10)*R_z)*q**3*scalar(1,32)))
    a_terms = (-scalar(3,16)/(p**3), -scalar(15,16)/p.square(), scalar(15,16)/p, scalar(3,16))
    b_terms = (scalar(15,32)/(p**3), -scalar(15,32)/p.square(), -scalar(15,32)/p, scalar(15,32))
    z_terms = (-scalar(5,16)/(p**3), scalar(15,16)/p.square(), -scalar(15,16)/p, scalar(5,16))
    proportional = scalar(3,2)*one_fifth_power_complement(q)
    even = (
        (a,tuple(zip((0,2,4,6),a_terms))),
        (b,tuple(zip((0,2,4,6),b_terms))),
        (D,tuple(zip((0,2,4,6),(x+proportional*y for x,y in zip(b_terms,z_terms))))),
        (E,tuple(zip((0,2,4,6),z_terms))),
    )
    odd = (
        (f, ((1,-scalar(15,16)*r*d_over_q/(lam**3)),
             (3,-scalar(15,8)*y/(lam**3)),
             (5,scalar(15,16)*(w+y)/lam))),
        (epsilon, ((1,-scalar(3,4)/(lam**5)),(5,scalar(3,4)/lam))),
    )
    return nominal, even+odd


KNOTS = {
    280: (1224,1010,568,422),
    300: (1240,1017,572,423),
    320: (1259,1024,576,424),
    340: (1278,1030,579,424),
    350: (1287,1036,581,425),
    360: (1298,1043,583,426),
    380: (1318,1051,587,427),
    400: (1340,1062,591,429),
}
BANDS = ((280,340,Q(1,8)),(340,400,Q(1,8)))


def ray_value(left,right,index,p):
    a,b=KNOTS[left][index],KNOTS[right][index]
    return scalar(a,1000)+scalar(b-a,right-left)*(p-scalar(left,1000))


def verify_ray(order,index,left,right,cap):
    cap_i=rational(cap)
    log_ratio=log(ONE+scalar(25,6)*cap_i)/log(ONE+scalar(4)*cap_i)
    rho=I(log_ratio.lo,(ONE-log(ONE-cap_i*scalar(1,6))/log(ONE+scalar(4)*cap_i)).hi) if order==6 else ONE
    cap_i=scalar(5)*cap_i/(scalar(6)-cap_i) if order==6 else cap_i
    log_cap=scalar(1,2)*log(ONE+scalar(order-1)*cap_i)
    pmax=Q(right,1000)
    bmax=Q(max(KNOTS[left][index],KNOTS[right][index]),1000)
    bound=pmax*bmax**2
    require(bound<1,"left physical ray")
    c=(bound-1)/(bound+1)
    limit=63/(65-16*c)
    initial=(interval(Q(left,1000),Q(right,1000)),rho,interval(Q(0),limit))
    widths=[max(v.hi-v.lo,1) for v in initial]

    @lru_cache(maxsize=8192)
    def get_coefficients(plo,phi,mlo,mhi):
        p=I(plo,phi)
        lam=p.sqrt()
        if order==5:
            return lam,contributions_five(p,lam,cap_i,log_cap)
        return lam,contributions_six(p,lam,I(mlo,mhi),cap_i,log_cap)

    stack=[(initial,0)]
    count=0
    while stack:
        (p,mu,xi),depth=stack.pop()
        count+=1
        lam, (nominal, families)=get_coefficients(p.lo,p.hi,mu.lo,mu.hi)
        b=ray_value(left,right,index,p)
        x,Y=xi/lam,imaginary_coordinate(p,b,xi)
        for directions in projections(order,x,Y):
            value=linear_form(nominal,directions)
            for parameter,terms in families:
                value=value+parameter*linear_form(terms,directions)
            if value.lo>0 or value.hi<0:
                break
        else:
            require(depth<41,("unresolved ray",order,index,left,right,[v.shown() for v in (p,mu,xi)]))
            box=(p,mu,xi)
            weights=(160,1 if depth>=16 else 0,32)
            k=max(range(3),key=lambda i: (box[i].hi-box[i].lo)*weights[i]*GRID//widths[i])
            old=box[k]
            for half in bisect(old):
                stack.append((box[:k]+(half,)+box[k+1:],depth+1))
    return count


def verify_diffuse_constants():
    verify_scalar_constants()
    require(KNOTS[280] == EARLIER_KNOTS[280], "matching sector anchor")
    require(all(all(a>b>0 for a,b in zip(rays,rays[1:])) for rays in KNOTS.values()),
            "strictly ordered ray knots")
    require(Q(3,8)**5 < Q(17,20)**5*Q(57,100)**7
            and 1 < Q(71,50)**5*Q(57,100)**3, "central derivative caps")
    z = Q(3,50)
    require(1/(2*(1-z*z/12)) < Q(501,1000)
            and 1/(6*(1-z*z/20)) < Q(1,5)
            and z*z/(6*(1-z*z/20)) < Q(1,1000)
            and Q(501,1000)*z*z < Q(1,500), "hyperbolic caps")
    require(Q(18,25)/(1-z)<Q(19,20), "interior multiplier cap")


def main():
    verify_diffuse_constants()
    for low,high,cap in BANDS:
        if len(argv)>1 and int(argv[1])!=low:
            continue
        knots=[k for k in KNOTS if low<=k<=high]
        for left,right in zip(knots,knots[1:]):
            out=[verify_ray(order,index,left,right,cap) for order in (5,6) for index in range(4)]
            print(f"Outer deficit [{left}, {right}]/1000: {out}",flush=True)


if __name__ == "__main__":
    main()
