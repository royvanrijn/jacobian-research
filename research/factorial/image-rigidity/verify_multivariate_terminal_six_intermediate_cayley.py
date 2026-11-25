#!/usr/bin/env python3
"""Exact angular separation on 18/25 <= u <= 9/10 at terminal six."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial, isqrt

from verify_multivariate_terminal_six_private_angles import require


GRID = 2**48


class I:
    """Closed intervals with endpoints on the dyadic lattice 1/GRID."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def __add__(self, other):
        return I(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        a, b, c, d = self.lo, self.hi, other.lo, other.hi
        terms = (a*c, a*d, b*c, b*d)
        return I(min(terms)//GRID, -(-max(terms)//GRID))

    def __truediv__(self, other):
        require(other.lo > 0, "positive denominator")
        return self * I(GRID**2//other.hi, -(-GRID**2//other.lo))

    def __pow__(self, exponent):
        if exponent == 0:
            return ONE
        result = self
        for _ in range(1, exponent):
            result = result * self
        return result

    def square(self):
        lower = min(abs(self.lo), abs(self.hi)) if self.lo * self.hi > 0 else 0
        upper = max(abs(self.lo), abs(self.hi))
        return I(lower*lower//GRID, -(-upper*upper//GRID))

    def sqrt(self):
        require(self.lo >= 0, "nonnegative radical")
        return I(isqrt(self.lo * GRID), isqrt(self.hi * GRID) + 1)

    def positive_error(self):
        require(self.hi >= 0, "nonnegative bound")
        return I(0, self.hi)

    def symmetric_error(self):
        require(self.hi >= 0, "nonnegative bound")
        return I(-self.hi, self.hi)

    def shown(self):
        return self.lo/GRID, self.hi/GRID


@lru_cache(maxsize=None)
def scalar(numerator, denominator=1):
    return I(numerator*GRID//denominator, -(-numerator*GRID//denominator))


def rational(value):
    return scalar(value.numerator, value.denominator)


def interval(lower, upper):
    return I(rational(lower).lo, rational(upper).hi)


ZERO, ONE = scalar(0), scalar(1)


def log(value):
    """Log of a positive interval via the atanh power sum."""
    z = (value - ONE) / (value + ONE)
    square = z.square()
    term = z
    total = z
    for k in range(1, 12):
        term = term * square
        total = total + term * scalar(1, 2*k+1)
    tail = (term * square * scalar(2, 25) / (ONE - square))
    return scalar(2) * total + I(-max(abs(tail.lo), abs(tail.hi)), max(abs(tail.lo), abs(tail.hi)))


def third_tail(a, deficit):
    """Binomial tail after degree two divided by deficit cubed, 1<a<2."""
    coefficient = a * (a-ONE) * (scalar(2)-a) * scalar(1, 6)
    term = coefficient
    total = term
    for k in range(3, 10):
        term = term * (scalar(k)-a) * scalar(1, k+1) * deficit
        total = total + term
    term = term * (scalar(10)-a) * scalar(1, 11) * deficit
    return total + (term / (ONE-deficit)).positive_error()


def one_fifth_power_complement(q):
    """Enclose (1-q)^(1/5) by the positive series for its deficit."""
    a = scalar(1, 5)
    term = a * q
    total = term
    for k in range(1, 9):
        term = term * (scalar(k)-a) * scalar(1, k+1) * q
        total = total + term
    next_term = term * (scalar(9)-a) * scalar(1, 10) * q
    return ONE - total - (next_term/(ONE-q)).positive_error()


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
    require((outer+link).hi < scalar(1,25).lo and ((scalar(5,4) if sextic else scalar(2))*inner).hi < scalar(1,25).lo,
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
    epsilon = (u*(scalar(1001,1000)*link+scalar(101,100)*chi**3)).symmetric_error()
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
    require(w.lo > scalar(7,10).hi, "w lower bound for central defect")
    W = rho*(-log(ONE-p))
    total, beta, eta, link, Y = skew_bounds(W,w,cap,log_cap,True)
    J = scalar(6,5)*cap**4*total
    R_y, R_z = third_tail(scalar(8,5),q), third_tail(scalar(9,5),q)
    y0 = ONE-scalar(8,5)*q+scalar(12,25)*q.square()+R_y*q**3
    D = (y0*Y/(ONE-Y)).positive_error()
    y = y0 + D
    d_over_q = scalar(3,5)-scalar(12,25)*q-R_y*q.square()-D/q
    a, b = scalar(501,1000)*w*beta.square(), scalar(501,1000)*y*eta.square()
    E = I(-(scalar(5,4)*J).hi,(scalar(13,20)*D.square()).hi)
    f = (scalar(1001,1000)*eta).symmetric_error()
    epsilon = (w*(scalar(1001,4000)*link+eta**3)).symmetric_error()
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


def imaginary_coordinate(p, b, xi):
    t = xi.square()
    numerator = b*t
    factor = p*b
    Y = b*(ONE-t)
    Y = I(Y.lo,b.hi)
    for _ in range(3):
        update = b - numerator/(ONE+factor*Y)
        Y = I(max(Y.lo,update.lo),min(Y.hi,update.hi))
    return Y


def projections(order, x, Y):
    modulus = x.square()+Y.square()
    re, ims = ONE,[ZERO]
    for k in range(order):
        re, im = re*x-ims[-1]*Y, re*Y+ims[-1]*x
        ims.append(im)
    moduli = [ONE]
    for _ in range(order):
        moduli.append(moduli[-1]*modulus)
    for h in range(order+1):
        yield [moduli[min(k,h)]*(ims[k-h] if k>=h else -ims[h-k]) for k in range(order+1)]


def linear_form(terms,directions):
    projected = ZERO
    for degree, coefficient in terms:
        projected = projected + coefficient*directions[degree]
    return projected


# p=1-u followed by scaled ray parameters. Linear interpolation is exact.
KNOTS = {
    100: (1092,940,545,404),
    120: (1106,950,548,406),
    140: (1119,959,550,408),
    143: (1120,960,551,408),
    160: (1130,969,551,410),
    180: (1144,975,554,412),
    200: (1158,982,555,414),
    220: (1175,989,559,416),
    240: (1189,994,562,418),
    250: (1195,999,564,419),
    260: (1208,1004,565,420),
    280: (1224,1010,568,422),
}
BANDS = ((100,143,Q(1,16)),(143,250,Q(1,11)),(250,280,Q(27,250)))


def ray_value(left,right,index,p):
    a,b=KNOTS[left][index],KNOTS[right][index]
    return scalar(a,1000)+scalar(b-a,right-left)*(p-scalar(left,1000))


def bisect(item):
    mid=(item.lo+item.hi)//2
    return I(item.lo,mid),I(mid,item.hi)


def verify_ray(order,index,left,right,cap):
    cap_i=rational(cap)
    log_ratio=log(ONE+scalar(25,6)*cap_i)/log(ONE+scalar(4)*cap_i)
    rho=I(log_ratio.lo,(ONE-log(ONE-cap_i*scalar(1,6))/log(ONE+scalar(4)*cap_i)).hi) if order==6 else ONE
    cap_i=scalar(5)*cap_i/(scalar(6)-cap_i) if order==6 else cap_i
    log_cap=scalar(1,2)*log(ONE+scalar(order-1)*cap_i)
    initial=(interval(Q(left,1000),Q(right,1000)),rho,interval(Q(0),Q(15,17)))
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





# Deficit of the residual order-five outer mean and physical ray cosines.
EXCEPTION_RAYS = {
    0:(-517,-682,-853,-934),
    15:(-493,-656,-853,-923),
    30:(-471,-630,-852,-912),
    45:(-448,-605,-851,-901),
    60:(-423,-583,-849,-890),
    71:(-404,-567,-848,-884),
}


def nth_root(value, order):
    """Grid enclosure of the positive order-th root of an exact rational."""
    numerator,denominator=value.numerator,value.denominator
    target=numerator*GRID**order
    left,right=0,GRID
    require(right**order*denominator>=target,"normalized row multiplier")
    while left+1<right:
        middle=(left+right)//2
        if middle**order*denominator<=target:
            left=middle
        else:
            right=middle
    return I(left,right)


@lru_cache(maxsize=None)
def marked_row(order):
    def coefficient(k):
        multiplier=Q(factorial(order+k)**order,
                     factorial(order)**(order-k)*factorial(2*order)**k)
        return scalar(comb(order,k))*nth_root(multiplier,order)
    return tuple(coefficient(k) for k in range(order+1))


@lru_cache(maxsize=None)
def cayley_basis(order,k):
    return tuple(scalar(sum((-1)**(k-j)*comb(k,j)*comb(order-k,h-j)
                            for j in range(max(0,h-order+k),min(k,h)+1)),2**order)
                 for h in range(order+1))


def cayley_transform(order,coefficients):
    entries=[ZERO]*(order+1)
    for k,a in enumerate(coefficients):
        entries=[value+a*b for value,b in zip(entries,cayley_basis(order,k),strict=True)]
    return tuple(enumerate(entries))


@lru_cache(maxsize=None)
def residual_anchor_constants(order):
    exponents = (0, 10, 15, 15, 10, 0) if order==5 else (0, 10, 16, 18, 16, 10, 0)
    marked=marked_row(order)
    return (
        tuple(scalar(g,10) for g in exponents),
        cayley_transform(order,marked),
        cayley_transform(order,[a*scalar(g,10) for a,g in zip(marked,exponents,strict=True)]),
    )


def exception_contributions(order,p,rho,cap,s_cap):
    exponents,leading,linear=residual_anchor_constants(order)
    marked=marked_row(order)
    exponents=tuple(a*rho for a in exponents)
    m=[]
    quad=[]
    for e in exponents:
        quad.append(e*(e-ONE)*scalar(1,2) if e.hi else ZERO)
        m.append(ONE-e*p+quad[-1]*p.square()+p**3*third_tail(e,p) if e.hi else ONE)
    secondary=cayley_transform(order,[a*(d+p*third_tail(e,p)) if e.hi else ZERO
                                         for a,e,d in zip(marked,exponents,quad,strict=True)])
    nominal=tuple((k,b-p*rho*c+p.square()*d)
                  for (k,b),(_,c),(_,d) in zip(leading,linear,secondary,strict=True))
    anchors=tuple(a*b for a,b in zip(marked,m,strict=True))
    strength=rho*(-log(ONE-p))
    strength=I(0,strength.hi)
    sem=cap*strength/s_cap
    outer=cap*sem
    link=cap**3*sem
    if order==5:
        link=scalar(3,5)*link
        inner=(outer+link)*scalar(1,2)
        defect=s_cap.square()*strength*scalar(1,6)
        f=(0,2,1,-1,-2,0)
        e=(0,1,0,0,-1,0)
        g=(0,0,1,1,0,0)
        h=(0,)*6
        J=ZERO
    else:
        outer=scalar(5,3)*outer
        link=scalar(6)*link*scalar(1,4)
        inner=(outer+link)*scalar(4,5)
        defect=s_cap.square()*strength*scalar(1,5)
        f=(0,5,4,0,-4,-5,0)
        e=(0,1,0,0,0,-1,0)
        g=(0,0,2,3,2,0,0)
        h=(0,0,0,-1,0,0,0)
        J=scalar(6,5)*cap**4*sem
    fden=1 if order==5 else 4
    gden=1 if order==5 else 2
    families=[]
    for data,denominator,scale in ((f,fden,inner.symmetric_error()),(e,1,link.symmetric_error()),
                                   (g,gden,defect),(h,1,J)):
        if scale.hi:
            families.append((scale,cayley_transform(order,[a*scalar(d,denominator)
                                 for a,d in zip(anchors,data,strict=True)])))
    for k,(fk,ek,gk,hk) in enumerate(zip(f,e,g,h,strict=True)):
        if k==0 or k==order:
            continue
        sigma=scalar(abs(fk),fden)*inner+scalar(abs(ek))*link+scalar(abs(gk),gden)*defect+scalar(abs(hk))*J
        # exp(t)-1-t is nonnegative and at most t^2/[2(1-|t|)].
        error=(scalar(1,2)*sigma.square()/(ONE-sigma)).positive_error()
        families.append((error,tuple((h,anchors[k]*b) for h,b in enumerate(cayley_basis(order,k)))))
    return nominal,families


def exception_ray(left,right,index,p):
    a,b=EXCEPTION_RAYS[left][index],EXCEPTION_RAYS[right][index]
    cos=scalar(a,1000)+scalar(b-a,right-left)*(p-scalar(left,1000))
    return ((ONE+cos)/(ONE-cos)).sqrt()


def verify_exception_ray(order,index,left,right):
    t=scalar(1,26)
    log_base=log(ONE+scalar(4)*t)
    rho=I((log(ONE+scalar(25,6)*t)/log_base).lo,
          (ONE-log(ONE-t*scalar(1,6))/log_base).hi) if order==6 else ONE
    cap=scalar(5)*t/(scalar(6)-t) if order==6 else t
    s_cap=scalar(1,2)*log(ONE+scalar(order-1)*cap)
    initial=interval(Q(left,1000),Q(right,1000)),rho,interval(Q(-15,17),Q(15,17))
    widths=[max(v.hi-v.lo,1) for v in initial]

    @lru_cache(maxsize=8192)
    def coefficients(plo,phi,rlo,rhi):
        return exception_contributions(order,I(plo,phi),I(rlo,rhi),cap,s_cap)

    stack=[(initial,0)]
    count=0
    while stack:
        (p,r,x),depth=stack.pop()
        count+=1
        nominal,families=coefficients(p.lo,p.hi,r.lo,r.hi)
        b=exception_ray(left,right,index,p)
        y=imaginary_coordinate(ONE,b,x)
        for directions in projections(order,x,y):
            total=linear_form(nominal,directions)
            for parameter,terms in families:
                total=total+parameter*linear_form(terms,directions)
            if total.lo>0 or total.hi<0:
                break
        else:
            require(depth<37,("exceptional unresolved ray",order,index,left,right,[a.shown() for a in (p,r,x)]))
            box=(p,r,x)
            weights=(120,1 if depth>=16 else 0,32)
            k=max(range(3),key=lambda i: (box[i].hi-box[i].lo)*weights[i]*GRID//widths[i])
            for half in bisect(box[k]):
                stack.append((box[:k]+(half,)+box[k+1:],depth+1))
    return count


def complex_mul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def modulus(z):
    return (z[0].square()+z[1].square()).sqrt()


def exception_anchor():
    points={5:((-504,774),(-924,481),(-1079,0)),6:((-443,792),(-837,588),(-1056,216))}
    sectors={5:(1,3,4),6:(0,2,4)}
    rays=EXCEPTION_RAYS[0]
    radius=scalar(1,500)
    for order in (5,6):
        marked=marked_row(order)
        for (px,py),sector in zip(points[order],sectors[order],strict=True):
            x,y=scalar(px,1000),scalar(py,1000)
            center=x,y
            powers=[(ONE,ZERO)]
            for _ in range(order):
                powers.append(complex_mul(powers[-1],center))
            local=[]
            for k in range(order+1):
                f=g=ZERO
                for h in range(k,order+1):
                    factor=scalar(comb(h,k))*marked[h]
                    a,b=powers[h-k]
                    f,g=f+factor*a,g+factor*b
                local.append(modulus((f,g)))
            remainder=local[0]
            for k in range(2,order+1):
                remainder=remainder+local[k]*radius**k
            require((radius*local[1]).lo>remainder.hi,("marked disk",order,px,py))
            for index,cos in enumerate(rays):
                c=scalar(cos,1000)
                s=(ONE-c.square()).sqrt()
                direction=c*y-s*x
                # A point is above an upper ray iff its argument is larger.
                require((direction.lo>radius.hi if index<sector else direction.hi<(-radius).lo),
                        ("marked sector",order,px,py,index,sector))
            if py:
                require(y.lo>radius.hi,"upper marked disk")
        full=points[order]+tuple((x,-y) for x,y in points[order] if y)
        for i,(x,y) in enumerate(full):
            for X,Y in full[:i]:
                require((x-X)**2+(y-Y)**2>16,"disjoint marked disks")


def exception_main():
    exception_anchor()
    knots=tuple(EXCEPTION_RAYS)
    for left,right in zip(knots,knots[1:]):
        out=[verify_exception_ray(order,index,left,right) for order in (5,6) for index in range(4)]
        print("exception",left,right,out,flush=True)


def verify_gaussian_anchor():
    p = scalar(1, 10)
    rho = scalar(25, 24)
    ordinary_cap = scalar(1, 16)
    for order, signs in ((5, (1, -1, -1, 1)), (6, (1, 1, -1, -1))):
        cap = ordinary_cap if order == 5 else scalar(5)*ordinary_cap/(scalar(6)-ordinary_cap)
        s = scalar(1, 2)*log(ONE+scalar(order-1)*cap)
        if order == 5:
            terms, _ = contributions_five(p, p.sqrt(), cap, s)
        else:
            terms, _ = contributions_six(p, p.sqrt(), rho, cap, s)
            require(terms[0][1].lo > 0, "Gaussian sextic at the origin")
        require(terms[-1][1].lo > 0, "Gaussian leading term")
        for value, sign in zip(KNOTS[100], signs, strict=True):
            b = scalar(value, 1000)
            evaluation = ZERO
            for degree, coefficient in terms:
                power = degree - (order % 2)
                evaluation = evaluation + scalar((-1)**(power//2))*coefficient*b**power
            require((scalar(sign)*evaluation).lo > 0, ("Gaussian anchor", order, value))


def verify_scalar_constants():
    # Outward radial bounds and the upper bound on the cosine of every ray.
    require(Q(19, 20)*(Q(9, 8)**6-1)+Q(1, 8)**5 < 1, "root annulus")
    require(Q(7, 25)*max(Q(r[0], 1000)**2 for r in KNOTS.values()) < Q(3, 7),
            "scaled ray cosine")
    for table in (KNOTS, EXCEPTION_RAYS):
        for rays in table.values():
            require(all(a>b for a,b in zip(rays, rays[1:])), "strict ray order")
    require(max(r[0] for r in EXCEPTION_RAYS.values()) < -400, "marked ray cosine")
    require(Q(4, 5) < Q(19, 20) and all(
        marked_row(j)[k].hi*5 < scalar(comb(j,k)).lo*4
        for j in (5, 6) for k in range(1,j)
    ), "marked interior coefficient cap")
    require(Q(3, 8)**5 < Q(13, 20)**5*Q(7, 10)**7
            and 1 < Q(5, 4)**5*Q(7, 10)**3, "central derivative caps")
    z = Q(1, 25)
    require(1/(2*(1-z*z/12)) < Q(501,1000)
            and 1/(6*(1-z*z/20)) < Q(1,5)
            and z*z/(6*(1-z*z/20)) < Q(1,1000)
            and Q(501,1000)*z*z < Q(1,1000), "hyperbolic caps")


def main():
    verify_scalar_constants()
    verify_gaussian_anchor()
    for low, high, cap in BANDS:
        knots = [k for k in KNOTS if low <= k <= high]
        for left, right in zip(knots, knots[1:]):
            out = [verify_ray(order, index, left, right, cap)
                   for order in (5, 6) for index in range(4)]
            print(f"Ordinary p in [{left}, {right}]/1000: {out}", flush=True)
    exception_main()


if __name__ == "__main__":
    main()
