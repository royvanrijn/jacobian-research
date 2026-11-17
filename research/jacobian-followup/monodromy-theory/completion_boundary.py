"""Exact boundary controls for the all-rank weighted completion."""
import sympy as S

def cycle_lengths(perm):
    remaining=set(range(len(perm)));out=[]
    while remaining:
        a=min(remaining);seen={a};b=perm[a]
        while b!=a:seen.add(b);b=perm[b]
        remaining-=seen;out.append(len(seen))
    return sorted(out)

def permutation_power(perm,n):
    out=tuple(range(len(perm)))
    for unused in range(n):out=tuple(perm[out[i]] for i in range(len(perm)))
    return out

def truncated_substitution(expression,old_variables,values,ring,degree):
    def cut(poly):
        return ring.from_dict({m:c for m,c in poly.items() if sum(m)<=degree})
    powers=[]
    source=S.Poly(expression,*old_variables,domain=S.QQ)
    maxima=[max((m[i] for m,c in source.terms()),default=0) for i in range(len(values))]
    for val,top in zip(values,maxima):
        cur=[ring.one]
        for j in range(top):cur.append(cut(cur[-1]*val))
        powers.append(cur)
    answer=ring.zero
    for mon,co in source.terms():
        term=ring.ground_new(ring.domain.convert(co))
        for axis,power in enumerate(mon):term=cut(term*powers[axis][power])
        answer+=term
    return cut(answer)

def binary_high_weight_data():
    from sympy.polys.rings import ring
    from branch_centered import binary_sweeps,binary_constants
    from higher_completion import weighted_stage,rectify
    w,z,g,u,v=S.symbols('w z gamma u v')
    rr=S.sqrt(-3);K=S.QQ.algebraic_field(rr)
    canonical=lambda x:K.to_sympy(K.from_sympy(x))
    X=binary_sweeps(w,z,g)['centered_X']
    stage=weighted_stage(X,[w,z],g,[u,v]);E=stage['E'];theta=[g,u,v]
    point=S.Matrix([4*rr/3,3*(1-rr)/16,S.Rational(3,8)])
    gold=4*rr/3
    Mi=S.Matrix([[1,0,0],[-point[1]/gold,1/gold,0],[-point[2]/gold,0,1/gold]])
    L=(Mi*binary_constants(rr)['L']*S.diag(1,gold,gold)).applyfunc(canonical)
    R,t,a,b=ring('t,a,b',K);new=[t,a,b]
    values=[R.ground_new(K.from_sympy(point[i]))+
            sum((R.ground_new(K.from_sympy(L[i,j]))*new[j] for j in range(3)),R.zero)
            for i in range(3)]
    jet=S.Matrix([truncated_substitution(expr,theta,values,R,4).as_expr() for expr in E])
    variables=[S.Symbol('t'),S.Symbol('a'),S.Symbol('b')]
    t0=variables[0];us=variables[1:];D=4
    def uni(expr):
        p=S.Poly(S.expand(expr),t0,domain=K)
        return S.Poly.from_dict({m:c for m,c in p.terms() if m[0]<D},(t0,),domain=K).as_expr()
    h=S.zeros(2,1)
    for iteration in range(D):
        value=jet[1:,:].subs(dict(zip(us,h)),simultaneous=True)
        h=(h-value).applyfunc(uni)
    M=jet[1:,:].jacobian(us).subs(dict(zip(us,h)),simultaneous=True).applyfunc(uni)
    pivot=S.Poly(M[1,1],t0,domain=K)
    coefficients=[K.from_sympy(pivot.nth(i)) for i in range(D)]
    inv=[K.one/coefficients[0]]
    for n in range(1,D):
        inv.append(-sum((coefficients[j]*inv[n-j] for j in range(1,n+1)),K.zero)/coefficients[0])
    numer=S.Poly(M[1,0],t0,domain=K)
    nc=[K.from_sympy(numer.nth(i)) for i in range(D)]
    quotient=[sum((nc[j]*inv[n-j] for j in range(n+1)),K.zero) for n in range(D)]
    entry=-sum((K.to_sympy(c)*t0**i for i,c in enumerate(quotient)),S.Integer(0))
    A=S.Matrix([[1,0],[entry,1]])
    phi=S.Matrix([t0,*list(h+A*S.Matrix(us))])
    inverse=S.Matrix([t0,*list(A.inv()*(S.Matrix(us)-h))])
    rect=dict(h=h,A=A,M=M,phi=phi,inverse=inverse,D=4,weights=[1,3,4])
    return dict(X=X,E=E,theta=theta,point=point,L=L,jet=jet,variables=variables,rectifier=rect,
                coefficient_field=K)
