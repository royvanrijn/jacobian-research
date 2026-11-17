"""Factorized homogeneous sweeps and an exact degree11 nonfinite completion."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
import sympy as S
def degree11_input():
    a,z,A=S.symbols('a z A');ps=list(S.symbols('p1:5'))
    R=S.Poly(330*a**4-165*a**3+55*a*a-11*a+1,a,domain=S.QQ)
    def red(e):
        num,den=S.fraction(S.cancel(e))
        return (S.Poly(num,a)*S.invert(S.Poly(den,a),R)).rem(R).as_expr()
    def redz(e):
        return S.Add(*[red(co)*z**mon[0] for mon,co in S.Poly(S.expand(e),z).terms()])
    c1=2-4*a;c2=12*a*a-6*a+2;c3=-30*a**3+15*a*a-5*a+1
    D3=z**3-c1*z*z+c2*z-c3
    C=S.Matrix([[0,0,c3],[1,0,-c2],[0,1,c1]])
    M=S.diag(a,C)
    Qp=ps[1]*S.eye(3)+ps[2]*C+ps[3]*C*C
    F=S.expand(ps[0]**7*Qp.det(method='domain-ge'))
    F=S.Add(*[red(co)*S.prod(p**j for p,j in zip(ps,mon)) for mon,co in S.Poly(F,*ps).terms()])
    f=-F;g=(F*(M-S.eye(4))*S.Matrix(ps)).applyfunc(S.expand)
    g=g.applyfunc(lambda e:S.Add(*[red(co)*S.prod(p**j for p,j in zip(ps,mon)) for mon,co in S.Poly(e,*ps).terms()]))
    Q=S.expand(z*(A+a*z)**7*(A**3+c1*A*A*z+c2*A*z*z+c3*z**3))
    Q=redz(Q)
    Qm=redz(Q.subs(A,-1));b=red(Qm.subs(z,1))
    tail,rem=S.div(Qm-b,(z-1)**5,z);tail=redz(tail);rem=redz(rem)
    resultant=S.Poly(S.resultant(R.as_expr(),tail,a),z).clear_denoms()[1].primitive()[1]
    lin=next(v for v in reversed(S.subresultants(R.as_expr(),tail,a)) if S.degree(v,a)==1)
    lc=S.Poly(lin,a).coeff_monomial(a);ct=S.Poly(lin,a).coeff_monomial(1)
    arec=(-S.Poly(ct,z)*S.invert(S.Poly(lc,z),resultant)).rem(resultant).as_expr()
    return dict(a=a,z=z,A=A,parameters=ps,R=R,red=red,redz=redz,C=C,M=M,F=F,
      X=S.Matrix([f,*g]),Q=Q,critical_value=b,tail=tail,tail_remainder=rem,
      resultant=resultant,a_recovery=arec,cubic=D3)

def completion_data(rectify=True,progress=None):
    from sympy.polys.rings import ring
    report=progress or (lambda message:None)
    inp=degree11_input();a=inp['a'];zsym=inp['z'];ps=inp['parameters'];X=inp['X']
    beta=S.CRootOf(inp['resultant'].as_expr(),0);K=S.QQ.algebraic_field(beta)
    bb=K.from_sympy(beta);one=K.one;zero=K.zero
    def qeval(expr,var,value):
        answer=zero
        for co in S.Poly(expr,var,domain=S.QQ).all_coeffs():answer=answer*value+K.convert(co)
        return answer
    aa=qeval(inp['a_recovery'],zsym,bb)
    def convert(expr):return qeval(inp['red'](expr),a,aa)
    M=[[convert(inp['M'][i,j]) for j in range(4)] for i in range(4)]
    def sparse(expr):
        return [(mon,convert(co)) for mon,co in S.Poly(expr,*ps).terms() if co!=0]
    polynomials=[sparse(e) for e in X]
    derivatives=[[sparse(S.diff(e,p)) for p in ps] for e in X]
    def evalpoly(poly,values):
        out=zero
        for mon,co in poly:
            term=co
            for val,n in zip(values,mon):term*=val**n
            out+=term
        return out
    def mv(A,x):return [sum((co*val for co,val in zip(row,x)),zero) for row in A]
    def mm(A,B):
        return [[sum((A[i][k]*B[k][j] for k in range(len(B))),zero) for j in range(len(B[0]))] for i in range(len(A))]
    def inv(A):
        n=len(A);aug=[row.copy()+[one if i==j else zero for j in range(n)] for i,row in enumerate(A)]
        det=one
        for j in range(n):
            k=next(i for i in range(j,n) if aug[i][j]!=zero)
            if k!=j:aug[k],aug[j]=aug[j],aug[k];det=-det
            co=aug[j][j];det*=co;aug[j]=[v/co for v in aug[j]]
            for i in range(n):
                if i==j:continue
                co=aug[i][j]
                aug[i]=[u-co*v for u,v in zip(aug[i],aug[j])]
        return [row[n:] for row in aug],det
    pstar=[one,one,zero,zero];star=[evalpoly(e,pstar) for e in polynomials]
    DXstar=[[evalpoly(e,pstar) for e in row] for row in derivatives]
    P0=[[DXstar[i+1][j]-pstar[i]*DXstar[0][j] for j in range(4)] for i in range(4)]
    power=P0
    for n in range(2,4):power=mm(power,P0)
    cyclic=next(j for j in range(4) if any(power[i][j]!=zero for i in range(4)))
    chain=[[one if i==cyclic else zero for i in range(4)]]
    for n in range(3):chain.append(mv(P0,chain[-1]))
    T=[[chain[3-j][i] for j in range(4)] for i in range(4)];Ti,detT=inv(T)
    # A0=-1, Z0=beta, and gamma0=beta-1.
    recover=[[bb*M[i][j]-(one if i==j else zero) for j in range(4)] for i in range(4)]
    recoveri,_=inv(recover);oldp=mv(recoveri,star[1:]);gamma0=bb-one
    w0=mv(Ti,[oldp[i]-pstar[i] for i in range(4)])
    weights=[1,5,6,7,8];powers=[4,5,6,7]
    theta0=[gamma0]+[w0[i]/gamma0**powers[i] for i in range(4)]
    def centered_numeric(values):
        old=[evalpoly(e,values) for e in polynomials]
        vec=[old[i+1]-pstar[i]*old[0]-(star[i+1]-pstar[i]*star[0]) for i in range(4)]
        return [old[0]-star[0],*mv(Ti,vec)]
    point_center=centered_numeric(oldp)
    marked_stage=[(point_center[0]+gamma0)/gamma0]+[
       (point_center[i+1]+gamma0*w0[i])/gamma0**weights[i+1] for i in range(4)]
    DXold=[[evalpoly(e,oldp) for e in row] for row in derivatives]
    targetDX=[DXold[0]]+mm(Ti,[[DXold[i+1][j]-pstar[i]*DXold[0][j] for j in range(4)] for i in range(4)])
    DXcenter=mm(targetDX,T)
    Sw=[[DXcenter[i][j]+(gamma0 if i==j+1 else zero) for j in range(4)] for i in range(5)]
    wg=[K.convert(powers[i])*w0[i]/gamma0 for i in range(4)]
    delta=[one,*w0]
    raw=[[delta[i]+sum((Sw[i][j]*wg[j] for j in range(4)),zero)]+[
        Sw[i][j]*gamma0**powers[j] for j in range(4)] for i in range(5)]
    J=[[co/gamma0**weights[i] for co in row] for i,row in enumerate(raw)]
    L,detJ=inv(J)
    report('degree24 point, Jordan frame and inverse stage jet computed')
    out=dict(input=inp,beta=beta,K=K,aa=aa,bb=bb,M=M,X_polynomials=polynomials,
      X_derivatives=derivatives,evalpoly=evalpoly,mv=mv,mm=mm,inv=inv,
      pstar=pstar,Ystar=star,DXstar=DXstar,P0=P0,T=T,Ti=Ti,detT=detT,
      old_marked_point=oldp,gamma0=gamma0,theta0=theta0,marked_stage=marked_stage,
      J=J,L=L,det_J=detJ,source_weights=powers,output_weights=weights)
    if not rectify:return out

    R,z=ring('z',K);D=8
    def ground(c):return R.ground_new(c)
    def cut(v):return R.from_dict({m:c for m,c in v.items() if m[0]<D})
    def mul(a,b):return cut(a*b)
    def inverse_unit(v):
        co=[v.get((i,),zero) for i in range(D)];answer=[one/co[0]]
        for n in range(1,D):
            answer.append(-sum((co[j]*answer[n-j] for j in range(1,n+1)),zero)/co[0])
        return R.from_dict({(i,):co for i,co in enumerate(answer) if co!=zero})
    def lvec(A,x):
        return [sum((ground(co)*val for co,val in zip(row,x)),R.zero) for row in A]
    def cleft(A,B):
        return [[sum((ground(A[i][k])*B[k][j] for k in range(len(B))),R.zero)
            for j in range(len(B[0]))] for i in range(len(A))]
    def cright(B,A):
        return [[sum((B[i][k]*ground(A[k][j]) for k in range(len(A))),R.zero)
            for j in range(len(A[0]))] for i in range(len(B))]
    allpolys=polynomials+[v for row in derivatives for v in row]
    maxima=[max(mon[i] for poly in allpolys for mon,co in poly) for i in range(4)]
    def series_stage(theta,with_jac=False):
        gam=theta[0];us=theta[1:];gp=[R.one];gi=[R.one];gaminv=inverse_unit(gam)
        for j in range(1,9):gp.append(mul(gp[-1],gam));gi.append(mul(gi[-1],gaminv))
        w=[mul(gp[powers[i]],us[i]) for i in range(4)]
        oldp=[ground(pstar[i])+sum((ground(T[i][j])*w[j] for j in range(4)),R.zero) for i in range(4)]
        pp=[]
        for val,top in zip(oldp,maxima):
            cur=[R.one]
            for j in range(top):cur.append(mul(cur[-1],val))
            pp.append(cur)
        def ev(poly):
            answer=R.zero
            for mon,co in poly:
                term=ground(co)
                for i,pow in enumerate(mon):term=mul(term,pp[i][pow])
                answer+=term
            return cut(answer)
        oldX=[ev(poly) for poly in polynomials]
        vx=[oldX[i+1]-ground(pstar[i])*oldX[0]-ground(star[i+1]-pstar[i]*star[0]) for i in range(4)]
        center=[oldX[0]-ground(star[0]),*lvec(Ti,vx)]
        swept=[center[0]+gam]+[center[i+1]+mul(gam,w[i]) for i in range(4)]
        E=[mul(gi[weights[i]],v) for i,v in enumerate(swept)]
        if not with_jac:return E
        oldDX=[[ev(poly) for poly in row] for row in derivatives]
        vx=[[oldDX[i+1][j]-ground(pstar[i])*oldDX[0][j] for j in range(4)] for i in range(4)]
        tarDX=[oldDX[0],*cleft(Ti,vx)];cenDX=cright(tarDX,T)
        Sw=[[cenDX[i][j]+(gam if i==j+1 else R.zero) for j in range(4)] for i in range(5)]
        wg=[ground(K.convert(powers[i]))*mul(gp[powers[i]-1],us[i]) for i in range(4)]
        delta=[R.one,*w]
        raw=[[delta[i]+sum((mul(Sw[i][j],wg[j]) for j in range(4)),R.zero)]+[
            mul(Sw[i][j],gp[powers[j]]) for j in range(4)] for i in range(5)]
        jac=[]
        for i in range(5):
            row=[mul(gi[weights[i]],v) for v in raw[i]]
            row[0]-=ground(K.convert(weights[i]))*mul(gaminv,E[i])
            jac.append([cut(v) for v in row])
        return E,jac
    def values(h):
        return [ground(theta0[i])+ground(L[i][0])*z+
          sum((ground(L[i][j+1])*h[j] for j in range(4)),R.zero) for i in range(5)]
    h=[R.zero]*4
    for iteration in range(D):
        E=series_stage(values(h));h=[cut(h[i]-E[i+1]) for i in range(4)]
        if iteration in [1,3,5,7]:report('finite curve iteration '+str(iteration+1)+'/8')
    vals=values(h);E,jac=series_stage(vals,with_jac=True)
    matrix=[[cut(sum((jac[i+1][k]*ground(L[k][j+1]) for k in range(5)),R.zero))
             for j in range(4)] for i in range(4)]
    U=[row.copy() for row in matrix]
    correction=[[R.one if i==j else R.zero for j in range(4)] for i in range(4)]
    for i in range(3,0,-1):
        for j in range(i):
            co=mul(U[i][j],inverse_unit(U[i][i]))
            for row in range(4):
                U[row][j]=cut(U[row][j]-mul(co,U[row][i]))
                correction[row][j]=cut(correction[row][j]-mul(co,correction[row][i]))
    report('finite lower-unitriangular rectifier completed')
    out.update(series_ring=R,z=z,cut=cut,mul=mul,inverse_unit=inverse_unit,
       series_stage=series_stage,curve_h=h,curve_values=vals,curve_residual=E[1:],
       M_high=matrix,U=U,correction_A=correction)
    return out
