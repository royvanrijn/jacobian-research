"""Polynomial triangular normalized sweeps from an integrated power."""
import sympy as S

def moment_target(rank,Y):
    T=S.Dummy('T');u=sum(Y[i]*T**(i+1) for i in range(rank))
    F=-S.expand(sum((-1)**(j+1)*u**j/S.Integer(j) for j in range(1,rank+2))).coeff(T,rank+1)
    return S.expand(F-F.subs({Y[j]:0 for j in range(1,rank)},simultaneous=True))

def family(rank,h,t=None):
    t=S.Symbol('t') if t is None else t;g=S.Symbol('gamma')
    ps=[t,*S.symbols('p2:'+str(rank+1))]
    unknown=list(S.symbols('q2:'+str(rank+1)))
    target_variables=[S.Symbol('A'),*S.symbols('Y2:'+str(rank+1))]
    H=moment_target(rank,target_variables)
    A=g+h
    Ys=[A]+[g*ps[i-1]+unknown[i-1] for i in range(1,rank)]
    primitives=[t]+[S.integrate(h**j,t) for j in range(1,rank+1)]
    integral=sum((-1)**j*S.binomial(rank,j)*target_variables[0]**(rank-j)*primitives[j] for j in range(rank+1))
    expr=S.expand(H.subs(dict(zip(target_variables,Ys)),simultaneous=True)+
                  (-1)**(rank-1)*integral.subs(target_variables[0],A))
    solved={}
    for j in range(2,rank+1):
        power=rank+1-j;want=ps[-1] if j==rank else 0
        equation=S.expand(expr.subs(solved)).coeff(g,power)-want
        unknown_j=unknown[j-2]
        coefficient=S.diff(equation,unknown_j)
        if coefficient not in [1,-1]:raise AssertionError('non-unit triangular coefficient')
        solved[unknown_j]=S.expand(unknown_j-equation/coefficient)
    X=S.Matrix([h,*[solved[q] for q in unknown],S.expand(expr.subs(solved).subs(g,0))])
    # For rank1 the unknown list is empty and the coefficient gamma is supplied by I directly.
    if rank==1:X=S.Matrix([h,S.integrate(t*S.diff(h,t),t)])
    sweep=X+g*S.Matrix([1,*ps])
    return dict(rank=rank,h=h,t=t,g=g,parameters=ps,X=X,S=sweep,
       target_variables=target_variables,target_polynomial=H,integral=integral,
       degree=rank*int(S.degree(h,t))+1)

def degree9_data(rectify=True):
    """Full-rank A9 input and a deterministic number-field completion recipe."""
    from sympy.polys.rings import ring
    from higher_completion import weighted_stage
    t=S.Symbol('t');f=family(4,t*t,t);ps=f['parameters'];g=f['g'];X=f['X']
    pstar=S.Matrix([1,0,0,0]);sstar=dict(zip(ps,pstar));Ystar=X.subs(sstar)
    N=S.Matrix.hstack(-S.Matrix(ps),S.eye(4));P=N*X.jacobian(ps);P0=P.subs(sstar)
    cyclic=S.eye(4)[:,3];T=S.Matrix.hstack(P0**3*cyclic,P0**2*cyclic,P0*cyclic,cyclic)
    Ti=T.inv();newps=list(S.symbols('w1:5'))
    oldpoint=T*S.Matrix(newps)+pstar
    oldX=X.subs(dict(zip(ps,oldpoint)),simultaneous=True)
    centered=S.Matrix([oldX[0]-Ystar[0],*list(Ti*(oldX[1:,:]-pstar*oldX[0]-(Ystar[1:,:]-pstar*Ystar[0])))])
    centered=centered.applyfunc(S.expand);us=list(S.symbols('u1:5'))
    st=weighted_stage(centered,newps,g,us);E=st['E'];theta=[g,*us]
    rootpoly=35*t**4+175*t**3+345*t*t+325*t+128
    alpha=S.CRootOf(rootpoly,0);K=S.QQ.algebraic_field(alpha);aa=K.from_sympy(alpha)
    def evaluate(expr,variables,values):
        poly=S.Poly(expr,*variables,domain=S.QQ);out=K.zero
        for mon,co in poly.terms():
            term=K.convert(co)
            for value,power in zip(values,mon):term*=value**power
            out+=term
        return out
    gamma0=K.one-aa*aa;oldp=[aa]
    for i in range(1,4):
        remainder=X[i]-ps[i]
        value=evaluate(remainder,ps[:i],oldp)
        oldp.append(K.convert(Ystar[i])-gamma0*oldp[i-1]-value)
    pp=[]
    for i in range(4):
        pp.append(sum((K.convert(Ti[i,j])*(oldp[j]-K.convert(pstar[j])) for j in range(4)),K.zero))
    theta0=[gamma0]+[pp[i]/gamma0**st['source_weights'][i] for i in range(4)]
    J=[[evaluate(S.diff(E[i],theta[j]),theta,theta0) for j in range(5)] for i in range(5)]
    aug=[row.copy()+[K.one if i==j else K.zero for j in range(5)] for i,row in enumerate(J)]
    det=K.one
    for col in range(5):
        pivot=next(i for i in range(col,5) if aug[i][col]!=K.zero)
        if pivot!=col:aug[pivot],aug[col]=aug[col],aug[pivot];det=-det
        value=aug[col][col];det*=value;aug[col]=[v/value for v in aug[col]]
        for i in range(5):
            if i==col:continue
            co=aug[i][col]
            aug[i]=[a-co*b for a,b in zip(aug[i],aug[col])]
    L=[row[5:] for row in aug]
    out=dict(family=f,T=T,P0=P0,pstar=pstar,Ystar=Ystar,centered_X=centered,centered_parameters=newps,
       E=E,theta=theta,rootpoly=rootpoly,alpha=alpha,K=K,old_marked_point=oldp,theta0=theta0,
       J=J,L=L,det_J=det,evaluate=evaluate,stage=st)
    if not rectify:return out
    R,z=ring('z',K);D=8
    def cut(p):return R.from_dict({m:c for m,c in p.items() if m[0]<D})
    def evaluate_series(expr,values):
        poly=S.Poly(expr,*theta,domain=S.QQ)
        maximum=[max((m[i] for m,c in poly.terms()),default=0) for i in range(5)]
        powers=[]
        for value,top in zip(values,maximum):
            pp=[R.one]
            for j in range(top):pp.append(cut(pp[-1]*value))
            powers.append(pp)
        result=R.zero
        for mon,co in poly.terms():
            term=R.ground_new(K.convert(co))
            for i,power in enumerate(mon):term=cut(term*powers[i][power])
            result+=term
        return cut(result)
    h=[R.zero for i in range(4)]
    def values(h):
        return [R.ground_new(theta0[i])+R.ground_new(L[i][0])*z+
          sum((R.ground_new(L[i][j+1])*h[j] for j in range(4)),R.zero) for i in range(5)]
    for iteration in range(D):
        vals=values(h);err=[evaluate_series(E[i+1],vals) for i in range(4)]
        h=[cut(a-b) for a,b in zip(h,err)]
    vals=values(h)
    residual=[evaluate_series(E[i+1],vals) for i in range(4)]
    M=[]
    for i in range(4):
        M.append([sum((evaluate_series(S.diff(E[i+1],theta[k]),vals)*L[k][j+1]
                       for k in range(5)),R.zero) for j in range(4)])
    M=[[cut(e) for e in row] for row in M];U=[row.copy() for row in M]
    A=[[R.one if i==j else R.zero for j in range(4)] for i in range(4)]
    def inverse_unit(p):
        co=[p.get((i,),K.zero) for i in range(D)];inv=[K.one/co[0]]
        for n in range(1,D):
            inv.append(-sum((co[j]*inv[n-j] for j in range(1,n+1)),K.zero)/co[0])
        return R.from_dict({(i,):v for i,v in enumerate(inv) if v!=K.zero})
    for i in range(3,0,-1):
        for j in range(i):
            co=cut(U[i][j]*inverse_unit(U[i][i]))
            for row in range(4):
                U[row][j]=cut(U[row][j]-co*U[row][i])
                A[row][j]=cut(A[row][j]-co*A[row][i])
    out.update(series_ring=R,z=z,h=h,A=A,M=M,U=U,curve_residual=residual,
               values=vals,series_evaluate=evaluate_series,cut=cut)
    return out
