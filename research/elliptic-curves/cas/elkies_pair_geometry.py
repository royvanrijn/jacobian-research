"""Exact conics and pointed genus-one fibre products for the X1092 pair search.

Sage-only arithmetic. Norm-ten construction is the existing Euclidean formula;
its scalar is retained. New code supplies the odd-divisor reduction on the
quartic, pointed maps, and a finite-reduction nontorsion certificate. It never
calls rank(), a class-group algorithm, or a known-fibre point oracle.
"""
from __future__ import annotations

from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, gcd
from elkies_pair_policy import require


def polynomial_record(p):
    return list(map(str,p.list())) or ['0']


def function_record(f):
    return {'numerator':polynomial_record(f.numerator()),
            'denominator':polynomial_record(f.denominator())}


def function(field, record):
    ring=field.ring()
    return field(ring(record['numerator']))/ring(record['denominator'])


def load_parent(parent):
    ring=PolynomialRing(QQ,'s'); field=ring.fraction_field()
    a,b=(ring(parent[k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    curve=EllipticCurve(field,[a,b])
    def value(r):
        return field(ring(r['numerator_coefficients_low_to_high']))/ring(r['denominator_coefficients_low_to_high'])
    points=[curve([value(row['X']),value(row['Y'])]) for row in parent['sections']]
    gram=matrix(ZZ,parent['generic_height_gram'])
    require(len(points)==17 and gram.det()==1092 and gram.is_positive_definite(),'wrong X1092 parent')
    require(a.degree()<=8 and b.degree()<=12,'not the degree-8/12 parent')
    delta=4*a**3+27*b*b
    require(delta.degree()==24 and delta.gcd(delta.derivative()).degree()==0, 'not the squarefree 24-I1 chart')
    return curve,points,gram


def exact_division(a,b):
    q,r=a.quo_rem(b)
    require(not r,'nonexact polynomial division')
    return q


def odd_reduce(f,g,v):
    """An actual odd effective divisor, not just odd intersection parity, is input."""
    ring=f.parent(); g=g.monic();v%=g;steps=[]
    require(g.degree()>0 and g.degree()%2==1 and not (v*v-f)%g,'invalid odd divisor')
    while g.degree()>1:
        nxt=exact_division(v*v-f,g)
        require(0<nxt.degree()<g.degree() and nxt.degree()%2==1,
                'odd reduction needs another infinity chart')
        steps.append({'g':polynomial_record(g),'v':polynomial_record(v),'quotient':polynomial_record(nxt)})
        g=nxt.monic();v%=g
    t=-g[0]/g[1];y=v(t)
    require(y*y==f(t),'odd reduction lost the rational point')
    return t,y,steps


def finite_divisor(errors, ordinate, f, extra_denominators=()):
    ring=f.parent(); field=ring.fraction_field();ordinate=field(ordinate)
    errors=[field(e) for e in errors if e]
    require(errors,'intersection is not zero-dimensional')
    g=gcd([ring(e.numerator()) for e in errors]).monic()
    forbidden=ring(ordinate.denominator())
    for e in errors: forbidden*=ring(e.denominator())
    for d in extra_denominators: forbidden*=ring(d)
    while g.degree()>0:
        d=gcd(g,forbidden)
        if d.degree()==0: break
        g=exact_division(g,d).monic()
    require(g.degree()>0 and g.degree()%2==1,'no odd finite intersection divisor in this chart')
    v=(ring(ordinate.numerator())*ring(ordinate.denominator()).inverse_mod(g))%g
    require(not (v*v-f)%g,'intersection does not lift to the carrier')
    return g,v


def lift(c,t,w):
    ring=PolynomialRing(QQ,'s')
    h,m,b,k=(ring(c[key])(t) for key in ('h','m','b','k'))
    return (b+h*w)/2, -(h*k+m*w)/2


def inverse_ordinate(c,t,x,y):
    ring=PolynomialRing(QQ,'s')
    h,m,b,k=(ring(c[key]) for key in ('h','m','b','k'))
    d,u,v=h.xgcd(m);require(d.degree()==0,'nonprimitive Euclidean line')
    u/=d;v/=d
    return u(t)*(2*x-b(t))-v(t)*(2*y+h(t)*k(t))


def verify_conic(c,curve,points,gram):
    ring=curve.base_ring().ring();h,m,b,k,q=(ring(c[key]) for key in ('h','m','b','k','q'))
    word=vector(ZZ,c['word']);require(word*gram*word==10,'lost norm-ten marking')
    negative=-sum((n*p for n,p in zip(word,points)),curve(0))
    require(negative[0]==m*m/(h*h)-b, 'trace x-word does not match the conic')
    require(negative[1]==-m**3/h**3+3*m*b/(2*h)-h*k/2, 'trace y-word does not match the conic')
    x0,x1=b/2,h/2;y0,y1=-h*k/2,-m/2;a,beta=curve.a4(),curve.a6()
    require(y0*y0+y1*y1*q-x0**3-3*x0*x1*x1*q-a*x0-beta==0,'bad constant lift identity')
    require(2*y0*y1-3*x0*x0*x1-x1**3*q-a*x1==0,'bad linear lift identity')
    require(q.degree()==2 and q.discriminant()!=0 and gcd(h,m).degree()==0,'degenerate conic')
    require(gcd(q,ring(4*a**3+27*beta*beta)).degree()==0, 'branch meets bad fibre; needs correction-aware proof')
    field=curve.base_ring();T=function(field,c['T']);W=function(field,c['W'])
    require(W*W==q(T) and max(T.numerator().degree(),T.denominator().degree())==2,'bad conic parametrization')
    x,y=lift(c,T,W); x=field(T.denominator()**4*x)
    height=max(x.numerator().degree(),x.denominator().degree()+8)
    require(height==8, 'unexpected base-change height: cannot claim the Schur value 3')


def conic(curve,points,gram,row):
    ring=curve.base_ring().ring();field=curve.base_ring();s=ring.gen();word=vector(ZZ,row['word'])
    require(word*gram*word==10,'selected word is not norm ten')
    p=-sum((n*point for n,point in zip(word,points)),curve(0))
    den=ring(p[0].denominator());require(den.is_square(),'trace denominator is not square')
    h=ring(den.sqrt()).monic()
    require(h.degree()==3,'trace needs another base chart; retained as UNKNOWN')
    nx,ny=ring(p[0]*h*h),ring(p[1]*h**3)
    require(gcd(nx,h).degree()==0,'trace has uncancelled poles')
    m=(-ny*nx.inverse_mod(h*h))%(h*h)
    g=exact_division(m*nx+ny,h*h)
    b=exact_division(m*m-nx,h*h)
    k=exact_division(m*b-2*g,h*h)
    q=exact_division(4*m*k-3*b*b-4*ring(curve.a4()),h*h)
    require(q.degree()==2 and q.discriminant()!=0,'nonquadratic or singular residual')
    c={'mask':row['mask'],'word':row['word'],**{n:polynomial_record(v) for n,v in [('h',h),('m',m),('b',b),('k',k),('q',q)]}}
    proof=None
    if QQ(q[2]).is_square():
        root=QQ(q[2]).sqrt();T=field(q[0]-s*s)/(2*root*s-q[1]);W=root*T+s
        proof={'kind':'rational-infinity','root':str(root)}
    else:
        for i,point in enumerate(points):
            intersection=int(gram[i,i]-(gram*word)[i])
            if intersection<=0 or intersection%2==0:continue
            ordinate=inverse_ordinate(c,s,point[0],point[1])
            x,y=lift(c,s,ordinate)
            try:
                div,v=finite_divisor([point[0]-x,point[1]-y,ordinate*ordinate-q],ordinate,q,
                                    [point[0].denominator(),point[1].denominator()])
                t0,w0,steps=odd_reduce(q,div,v)
            except ValueError:
                continue
            T=field(t0)+(q.derivative()(t0)-2*w0*s)/(s*s-q[2]);W=w0+s*(T-t0)
            proof={'kind':'odd-section-divisor','section':i,'intersection':intersection,
                   'g':polynomial_record(div),'v':polynomial_record(v),'steps':steps,'point':[str(t0),str(w0)]}
            break
    require(proof is not None,'odd rationality is known, but this bounded chart did not construct its point')
    c.update(T=function_record(T),W=function_record(W),rationality=proof,
             coefficient_bits=sum(max(abs(v.numerator()).nbits(),v.denominator().nbits()) for v in q))
    verify_conic(c,curve,points,gram)
    return c


def pointed_model(f,t0,y0):
    ring=f.parent();s=ring.gen();shift=ring(f(s+t0));c=[QQ(shift[i]) for i in range(5)]
    require(c[0]==y0*y0 and f.gcd(f.derivative()).degree()==0,'singular or wrongly pointed quartic')
    if y0:
        a4=c[1]*c[3]-4*y0*y0*c[4]
        a6=y0*y0*c[3]**2+c[1]**2*c[4]-4*y0*y0*c[2]*c[4]
    else:
        require(c[1]!=0,'singular branch origin');a4=c[1]*c[3];a6=c[1]**2*c[4]
    J=EllipticCurve(QQ,[0,c[2],0,a4,a6])
    return J,c


def to_quartic(J,c,t0,y0,point):
    if point.is_zero():return t0,y0
    x,y=point.xy()
    if y0:
        special=J([c[1]**2/(4*y0*y0)-c[2],
                   -(c[1]*(c[1]**2/(4*y0*y0)-c[2])+2*y0*y0*c[3])/(2*y0)])
        if point==special:return t0,-y0
        den=x*x-4*y0*y0*c[4];require(den!=0,'exceptional inverse chart')
        u=(2*y0*y+c[1]*x+2*y0*y0*c[3])/den
        z=x*u*u/(2*y0)-y0-c[1]*u/(2*y0)
    else:
        require(x!=0,'exceptional branch inverse chart');u=c[1]/x;z=y*u*u/c[1]
    return t0+u,z


def from_quartic(J,c,t0,y0,t,z):
    u=t-t0
    if u==0:
        if z==y0:return J(0)
        require(y0 and z==-y0,'wrong origin ordinate')
        x=c[1]**2/(4*y0*y0)-c[2]
        return J([x,-(c[1]*x+2*y0*y0*c[3])/(2*y0)])
    if y0:
        x=(2*y0*(z+y0)+c[1]*u)/(u*u)
        y=((x*x-4*y0*y0*c[4])*u-c[1]*x-2*y0*y0*c[3])/(2*y0)
    else:x=c[1]/u;y=c[1]*z/(u*u)
    return J([x,y])


def verify_pointed_maps(f,t0,y0,J,c):
    """Check both rational maps in function fields, not at finitely many points."""
    ring=f.parent();field=ring.fraction_field();s=ring.gen()
    ar=PolynomialRing(field,'z');z=ar.gen();K=field.extension(z*z-field(f),'z');z=K.gen();u=K(s-t0)
    if y0:
        x=(2*y0*(z+y0)+c[1]*u)/(u*u)
        y=((x*x-4*y0*y0*c[4])*u-c[1]*x-2*y0*y0*c[3])/(2*y0)
        ui=(2*y0*y+c[1]*x+2*y0*y0*c[3])/(x*x-4*y0*y0*c[4])
        zi=x*ui*ui/(2*y0)-y0-c[1]*ui/(2*y0)
    else:
        x=c[1]/u;y=c[1]*z/(u*u);ui=c[1]/x;zi=y*ui*ui/c[1]
    require(y*y==x**3+J.a2()*x*x+J.a4()*x+J.a6(),'pointed Weierstrass identity failed')
    require(ui==u and zi==z,'maps are not birational inverses')


def nontorsion(curve,point):
    """Gcd of good reduction orders bounds rational torsion; no BSD/rank oracle."""
    orders=[];bound=ZZ(0)
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]:
        try:
            E=EllipticCurve(GF(p),[GF(p)(a) for a in curve.a_invariants()])
            n=ZZ(E.cardinality())
        except (ValueError,ZeroDivisionError,ArithmeticError):continue
        orders.append([p,int(n)]);bound=gcd(bound,n)
        if len(orders)>=3:break
    require(len(orders)>=2,'insufficient good reductions for a torsion bound')
    multiple=bound*point
    if multiple.is_zero():return None
    return {'point':list(map(str,point.xy())),'good_orders':orders,'torsion_order_divides':int(bound),
            'multiple':list(map(str,multiple.xy())),'rank_lower_bound':1}


def carrier(first,second):
    ring=PolynomialRing(QQ,'u');field=ring.fraction_field();u=ring.gen()
    q1,q2=ring(first['q']),ring(second['q'])
    require(q1.degree()==q2.degree()==2 and q1.resultant(q2)!=0,'overlapping branch divisors')
    T=function(field,first['T']);W=function(field,first['W']);D=field(T.denominator())
    f=ring(D*D*q2(T))
    require(f.degree()==4 and f.gcd(f.derivative()).degree()==0,'carrier needs another projective chart')
    origin_proof={}
    if QQ(f[4]).is_square():
        root=QQ(f[4]).sqrt();f=ring(list(reversed([f[i] for i in range(5)])))
        T=field(T(1/u));W=field(W(1/u));D=field(u*u*D(1/u))
        t0,y0=QQ(0),root
        origin_proof={'kind':'rational-infinity-chart','root':str(root)}
    else:
        x,y=lift(first,T,W);v=inverse_ordinate(second,T,x,y);x2,y2=lift(second,T,v);z=D*v
        div,ordinate=finite_divisor([x-x2,y-y2,z*z-f],z,f,
                                    [T.denominator(),W.denominator(),x.denominator(),y.denominator(),D.denominator()])
        t0,y0,steps=odd_reduce(f,div,ordinate)
        origin_proof={'kind':'odd-intersection-divisor','g':polynomial_record(div),
                      'v':polynomial_record(ordinate),'steps':steps}
    require(W*W==q1(T) and field(f)==D*D*q2(T),'fibre-product maps changed squareclasses')
    J,c=pointed_model(f,t0,y0);verify_pointed_maps(f,t0,y0,J,c)
    witnesses=[]
    if y0:witnesses.append(from_quartic(J,c,t0,y0,t0,-y0))
    # The other deck involution can supply a useful point when the first is torsion.
    try:
        s0,w10,w20=T(t0),W(t0),y0/D(t0)
        roots=ring((T-s0).numerator()).roots(QQ)
        for r,mult in roots:
            try:
                if W(r)==-w10:
                    z=D(r)*w20
                    if z*z==f(r):witnesses.append(from_quartic(J,c,t0,y0,r,z))
            except (ZeroDivisionError,ValueError):pass
    except (ZeroDivisionError,ValueError):pass
    proof=None
    for p in witnesses:
        if not p.is_zero():
            proof=nontorsion(J,p)
            if proof:break
    require(proof is not None,'deck-orbit points do not certify positive rank in this bounded attempt')
    return {'kind':'genus-one','masks':[first['mask'],second['mask']],
            'f':polynomial_record(f),'T':function_record(T),'W1':function_record(W),'D':function_record(D),
            'origin':[str(t0),str(y0)],'origin_proof':origin_proof,
            'curve':list(map(str,J.a_invariants())),'nontorsion':proof,
            'function_field_rank_lower_bound':19,
            'independence':'Two distinct quadratic characters, each with a nonzero anti-invariant section.',
            'boundary':'Positive base rank and generic subgroup rank are not a specialized M19 certificate.'}


def carrier_value(data,point):
    ring=PolynomialRing(QQ,'u');field=ring.fraction_field()
    f=ring(data['f']);t0,y0=map(QQ,data['origin']);J,c=pointed_model(f,t0,y0)
    r,z=to_quartic(J,c,t0,y0,point)
    require(z*z==f(r),'inverse carrier point failed')
    t=function(field,data['T'])(r);w1=function(field,data['W1'])(r);d=function(field,data['D'])(r)
    require(d!=0,'undefined second-sheet chart')
    return t,w1,z/d
