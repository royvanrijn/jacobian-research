"""Actual fixed germs and an inductive sharp dessin family; no Keller realizations."""
import sympy as S


def cycles(p):
    seen=set();out=[]
    for i in range(len(p)):
        if i in seen:continue
        row=[];j=i
        while j not in seen:
            seen.add(j);row.append(j);j=p[j]
        out.append(row)
    return out


def dessin_rotations(m):
    if type(m) is not int or m<1:raise ValueError('positive integer parameter required')
    sigma=[1,0,3,2];tau=[1,2,0,3];leaf=(0,1,2)
    for _ in range(m-1):
        a,b,c=leaf;aw,bz,w0,w1,z0,z1=range(len(sigma),len(sigma)+6)
        sigma += [a,b,w1,w0,z1,z0];sigma[a]=aw;sigma[b]=bz
        tau += [w0,z0,w1,aw,z1,bz]
        leaf=(w0,w1,aw)
    face=[sigma[tau[i]] for i in range(len(sigma))]
    return tuple(sigma),tuple(tau),tuple(face)


def sharp_profile(m):
    if type(m) is not int or m<1:raise ValueError('positive integer parameter required')
    n=3*m-1;e=5*m-2;ell=6*m-2
    return dict(parameter=m,quotient_order=n,marked_tangent_character=m,index=e,
                local_degree=ell,residue=1,excess=m,regular_fibre_order=1,
                remaining_cone_order=3,residual_Milnor=n*(m-1),
                etale_open='C* times C',whole_plane_Keller_realization=False,
                local_bound_is_sharp=True,global_Wright_bound_optimality_claimed=False)


def explicit_sharp_pair(m):
    t,X,Y=S.symbols('t X Y')
    if m==1:p=1+4*t;q=1+6*t+6*t*t
    elif m==2:
        p=1+12*t+60*t*t+96*t**3
        q=1+18*t+144*t*t+576*t**3+1080*t**4+432*t**5
    else:raise ValueError('only the two pinned explicit examples are expanded')
    n=3*m-1;P=S.expand(Y**(2*n)*p.subs(t,X/Y**3));Q=S.expand(Y**(3*n)*q.subs(t,X/Y**3))
    return dict(t=t,X=X,Y=Y,p=p,q=q,P=P,Q=Q,profile=sharp_profile(m))


def higher_characteristic_packet():
    x,y,u,v=S.symbols('x y u v')
    P=21*x*x+42*x*y**3-7*y**6+42*x*x*y+28*x*y**4-6*y**7;Q=x*y
    F=u**4-6272*u*u*v**3-917504*u*v**5-16777216*v**7+9834496*v**6
    A=x-y**3;H=S.cancel(F.subs({u:P,v:Q})/A**3)
    a1=-458752*u*v**3-733824*u*v**2+12500992*v**4
    b1=13*u*u-512*u*v*v-262144*v**4-529984*v**3
    c1=1835008*v**3+3179904*v*v
    a2=-91*u*u-3584*u*v*v+21952*v**3
    b2=-56*u*v-2048*v**3;c2=364*u+14336*v*v
    M=S.Matrix([[a1,a2],[b1,b2]])
    J=S.Matrix([[S.diff(P,x),S.diff(P,y)],[S.diff(Q,x),S.diff(Q,y)]])
    L=(J.inv()*M.subs({u:P,v:Q})).applyfunc(S.cancel)
    return dict(x=x,y=y,u=u,v=v,P=P,Q=Q,F=F,A=A,H=H,
                target_fields=M,target_log_factors=[-c1,-c2],Jacobian_matrix=J,lifted_fields=L)


def truncated_values(generators,t,bound):
    # This exact linear-algebra helper is used only after the written generator
    # order bounds show that EVERY omitted series term vanishes modulo t^bound.
    rows=[[S.expand(g).coeff(t,k) for k in range(bound)] for g in generators]
    matrix=S.Matrix(rows)
    _,pivots=matrix.rref()
    return list(pivots)
