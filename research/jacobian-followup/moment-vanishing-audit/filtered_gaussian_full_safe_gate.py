"""Finite coefficient gate for ALL filtered Gaussian safe phases, without purity.

Uses sealed finite Morse routines as pure functions; never invokes their mains.
The universal nilpotence/normalization theorem is a written proof, not inferred
from the finite examples. Abstract coefficient rings still need ideal oracles.
"""
import sympy as s
import verify_filtered_gaussian_critical_action as F0

t=F0.t


def filtered_normal_form(F,X,weights,B0):
    F=s.expand(s.sympify(F));X=tuple(X);weights=tuple(s.Rational(w) for w in weights);B0=s.Matrix(B0)
    if len(X)!=len(weights) or len(set(X))!=len(X):raise ValueError('distinct coordinates and matching rational weights required')
    if any(w<=0 or w>=2 for w in weights):raise ValueError('strictly positive residual weights below two required')
    if B0.shape!=(len(X),len(X)) or B0!=B0.T or B0.det()==0 or any(v.is_Rational is not True for v in B0):
        raise ValueError('fixed nonsingular rational symmetric Gaussian precision required')
    F0.polynomial_guard(F,X)
    names=sorted({str(v) for v in F.free_symbols|set(X)|{t}})
    if not X:
        invariants=F0.seed_invariants(t*F,s.S.One)
        return dict(variables=[],symbol_names=names,weights=[],phase=str(F),precision=[],critical_blocks=[],
                    critical_point=[],sigma='1',critical_action=str(t*F),quadratic_precision=[],
                    map=[],inverse_map=[],Morse_steps=[],Jacobian='1',actual_generic_rank_without_zero_positions=1,**invariants)
    for i in range(len(X)):
        for j in range(len(X)):
            if B0[i,j] and weights[i]+weights[j]!=2:raise ValueError('Gaussian precision is not weight-compatible')
    for monomial,c in s.Poly(F,*X).terms():
        if c and sum(a*w for a,w in zip(monomial,weights))>2:raise ValueError('phase exceeds weighted degree two')
    v=B0.inv()*s.Matrix([s.diff(F,x) for x in X])
    critical=[s.S.Zero]*len(X);sigma=s.S.One;blocks=[]
    for weight in sorted(set(weights)):
        I=[i for i,w in enumerate(weights) if w==weight]
        group=s.Matrix([X[i] for i in I])
        matrix=s.Matrix([[s.diff(v[i],X[j]).subs(dict.fromkeys(X,0)) for j in I] for i in I])
        rhs=(s.Matrix([v[i] for i in I])-matrix*group).applyfunc(s.expand)
        if any(X[j] in f.free_symbols and weights[j]>=weight for f in rhs for j in range(len(X))):
            raise ValueError('critical recursion is not lower-weight triangular')
        block=s.eye(len(I))-t*matrix
        values=t*block.inv()*rhs.subs(dict(zip(X,critical)),simultaneous=True)
        for i,value in zip(I,values):critical[i]=s.cancel(value)
        sigma=s.factor(sigma*block.det())
        blocks.append(dict(weight=str(weight),indices=I,matrix=[[str(a) for a in row] for row in matrix.tolist()],
                           right_side=[str(a) for a in rhs],solution=[str(s.cancel(a)) for a in values]))
    if any(s.cancel(c-t*vi.subs(dict(zip(X,critical)),simultaneous=True))!=0 for c,vi in zip(critical,v)):
        raise ValueError('critical point does not solve original gradient equations')
    q=(s.Matrix(X).T*B0*s.Matrix(X))[0]/2
    potential=s.expand(q-t*F)
    at_point=lambda f:s.cancel(f.subs(dict(zip(X,critical)),simultaneous=True))
    action=-at_point(potential)
    B=s.hessian(potential,X).subs(dict(zip(X,critical)),simultaneous=True).applyfunc(s.cancel)
    if s.cancel(B.det()/B0.det()-sigma)!=0:raise ValueError('moving Hessian determinant differs from block determinant')
    inverse=B.inv();target=(s.Matrix(X).T*B*s.Matrix(X))[0]/2
    current=F0.R.substitute(potential,X,[x+c for x,c in zip(X,critical)])+action
    current=F0.R.clean(current,X);psi=list(X);steps=[]
    for k in range(3,int(s.floor(2/min(weights)))+1):
        error=F0.R.degree_part(current,X,k)
        if error==0:continue
        delta=[s.cancel(a) for a in -inverse*s.Matrix([s.diff(error,x) for x in X])/k]
        for i,f in enumerate(delta):
            if any(x in f.free_symbols and weights[j]>=weights[i] for j,x in enumerate(X)):
                raise ValueError('Morse correction is not strictly lower-weight triangular')
        shift=[x+d for x,d in zip(X,delta)]
        current=F0.R.substitute(current,X,shift)
        psi=[F0.R.substitute(f,X,shift) for f in psi]
        steps.append(dict(degree=k,removed=str(error),correction=[str(a) for a in delta]))
    if F0.R.clean(current-target,X)!=0:raise ValueError('filtered finite Morse identity failed')
    theta=[F0.R.clean(f+c,X) for f,c in zip(psi,critical)];inv=list(X)
    for i in sorted(range(len(X)),key=lambda i:weights[i]):
        inv[i]=F0.R.clean(X[i]-(theta[i]-X[i]).xreplace(dict(zip(X,inv))),X)
    if any(F0.R.substitute(a,X,inv)!=b for a,b in zip(theta,X)) or any(F0.R.substitute(a,X,theta)!=b for a,b in zip(inv,X)):
        raise ValueError('wrong finite polynomial inverse')
    if s.cancel(s.Matrix(theta).jacobian(X).det()-1)!=0:raise ValueError('original volume not preserved')
    if any(s.cancel(a.subs(t,0)-b)!=0 for a,b in zip(theta,X)):raise ValueError('map not identity at time zero')
    if F0.R.clean(F0.R.substitute(potential,X,theta)-target+action,X)!=0:raise ValueError('wrong original-potential identity')
    return dict(variables=[str(a) for a in X],symbol_names=names,weights=[str(w) for w in weights],phase=str(F),
                precision=[[str(a) for a in row] for row in B0.tolist()],critical_blocks=blocks,
                critical_point=[str(a) for a in critical],sigma=str(s.factor(sigma)),critical_action=str(s.factor(action)),
                quadratic_precision=[[str(a) for a in row] for row in B.tolist()],map=[str(a) for a in theta],
                inverse_map=[str(a) for a in inv],Morse_steps=steps,Jacobian='1',actual_generic_rank_without_zero_positions=1,
                **F0.seed_invariants(action,sigma))


def coefficients(expressions):
    result=[]
    for expr in expressions:
        for _,c in s.Poly(s.expand(expr),t).terms():
            if c and c not in result:result.append(c)
    return result


def full_gate(H,Z,residual,X,weights,precision):
    """P=W.H(Z)+residual(Z,X); Gaussian Q0=Z.W+X^T precision X/2.

    Return polynomial generators I. Safe over ANY Q-algebra iff I is
    nilpotent after coefficient specialization. No purity premise is used.
    """
    H=tuple(s.sympify(h) for h in H);Z=tuple(Z);X=tuple(X)
    if len(H)!=len(Z) or len(set(Z+X))!=len(Z+X):
        raise ValueError('Distinct position/residual coordinates and matching H required')
    if t in Z+X or any(t in h.free_symbols or set(X)&h.free_symbols for h in H):
        raise ValueError('H must be time-independent and depend only on zero positions/parameters')
    for h in H:
        try:s.Poly(h,*sorted(h.free_symbols|set(Z),key=str))
        except s.PolynomialError as exc:raise ValueError('Polynomial H required') from exc
    model=filtered_normal_form(residual,X,weights,s.Matrix(precision))
    degree=[0 if not h else s.Poly(h,*Z).total_degree() for h in H]
    B=int(s.prod(d+1 for d in degree));q=[s.Integer(0)]*len(Z)
    for _ in range(B+1):
        q=[s.series(t*h.subs(dict(zip(Z,q)),simultaneous=True),t,0,B+1).removeO().expand() for h in H]
    substitution=dict(zip(Z,q))
    residuals=[(a-t*h.subs(substitution,simultaneous=True)).expand() for a,h in zip(q,H)]
    J=(s.eye(len(Z))-t*s.Matrix(H).jacobian(Z)).det() if Z else s.Integer(1)
    loc={str(a):a for a in (set().union(*(h.free_symbols for h in H))|s.sympify(residual).free_symbols|set(Z+X)|{t})}
    parse=lambda value:s.sympify(value,locals=loc)
    Jq=J.subs(substitution,simultaneous=True).expand()
    sigmaq=parse(model['sigma']).subs(substitution,simultaneous=True).expand()
    nuq=parse(model['action_numerator']).subs(substitution,simultaneous=True).expand()
    blocks=dict(fixedpoint_residual=coefficients(residuals),normal_Jacobian=coefficients([Jq-1]),
                residual_determinant=coefficients([sigmaq-1]),critical_action=coefficients([nuq]))
    generators=list(dict.fromkeys(c for values in blocks.values() for c in values))
    return dict(positions=[str(z) for z in Z],H=[str(h) for h in H],
        residual_model=model,polynomial_field_degree_bound=B,truncation=[str(a) for a in q],
        fixedpoint_residual=[str(a) for a in residuals],J_at_truncation=str(Jq),
        sigma_at_truncation=str(sigmaq),action_numerator_at_truncation=str(nuq),
        coefficient_blocks={k:[str(c) for c in values] for k,values in blocks.items()},
        safe_ideal_generators=[str(c) for c in generators],
        field_gate='all generators equal zero',all_ring_gate='generated ideal is nilpotent',
        normalization='Gamma_I(R)',purity_premise_required=False,
        safe_field_phases_are_pure=True,universal_claims_use_written_proof=True)
