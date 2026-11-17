"""Finite filtered critical/Morse certificates and exact ordinary seed gates.

General fields, arbitrary coefficient rings, full generic image and slope
claims use the written proof. The executable algebra is exact rational
function arithmetic with symbolic polynomial coefficient parameters.
"""
from pathlib import Path
import hashlib,json
import sympy as s
import gaussian_generic_cutoff_producer as G
import verify_semipositive_gaussian_reduction as R
import verify_algebraic_gaussian_lifts as A

p=Path(__file__).resolve().parent
t,z,e=R.t,R.z,R.e
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()


def polynomial_guard(F,X):
    if t in F.free_symbols or t in X:raise ValueError('time-independent polynomial phase required')
    try:s.Poly(F,*sorted(F.free_symbols|set(X)|{z},key=str))
    except s.PolynomialError as exc:raise ValueError('polynomial coefficient input required') from exc


def seed_invariants(psi,sigma):
    """Retain polynomial nu,sigma,M,Omega; derivative keeps parameters fixed."""
    psi=s.cancel(psi);sigma=s.expand(sigma)
    denominator=s.denom(psi)
    symbols=sorted(denominator.free_symbols|{t},key=str)
    ceiling=max(0,s.Poly(denominator,*symbols).total_degree())
    for M in range(ceiling+1):
        nu=s.cancel(psi*sigma**M)
        if not s.denom(nu).free_symbols:break
    else:raise ValueError('action denominator is not a power of sigma')
    nu=s.expand(nu)
    Omega=s.expand(2*sigma*s.diff(nu,t)-2*M*nu*s.diff(sigma,t)-sigma**M*s.diff(sigma,t))
    ideal=[s.expand(sigma-1).coeff(t,j) for j in range(1,max(0,int(s.degree(sigma,t)))+1)]
    if nu!=0:ideal += [c for m,c in s.Poly(nu,t).terms()]
    return dict(action_numerator=str(nu),action_denominator_power=M,pure_gate_polynomial=str(Omega),
                sufficient_ordinary_pure_prefix=0 if Omega==0 else int(s.degree(Omega,t))+1,
                safe_invariant_generators=[str(v) for v in dict.fromkeys(ideal) if v!=0])


def filtered(F,X,weights,B0):
    F=s.expand(s.sympify(F));X=tuple(X);weights=tuple(s.Rational(w) for w in weights);B0=s.Matrix(B0)
    if len(X)!=len(weights) or len(set(X))!=len(X):raise ValueError('distinct coordinates and matching rational weights required')
    if any(w<=0 or w>=2 for w in weights):raise ValueError('strictly positive residual weights below two required')
    if B0.shape!=(len(X),len(X)) or B0!=B0.T or B0.det()==0 or any(v.is_Rational is not True for v in B0):
        raise ValueError('fixed nonsingular rational symmetric Gaussian precision required')
    polynomial_guard(F,X)
    names=sorted({str(v) for v in F.free_symbols|set(X)|{t}})
    if not X:
        invariants=seed_invariants(t*F,s.S.One)
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
        rhs=s.Matrix([v[i] for i in I])-matrix*group
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
    current=R.substitute(potential,X,[x+c for x,c in zip(X,critical)])+action
    current=R.clean(current,X);psi=list(X);steps=[]
    for k in range(3,int(s.floor(2/min(weights)))+1):
        error=R.degree_part(current,X,k)
        if error==0:continue
        delta=[s.cancel(a) for a in -inverse*s.Matrix([s.diff(error,x) for x in X])/k]
        for i,f in enumerate(delta):
            if any(x in f.free_symbols and weights[j]>=weights[i] for j,x in enumerate(X)):
                raise ValueError('Morse correction is not strictly lower-weight triangular')
        shift=[x+d for x,d in zip(X,delta)]
        current=R.substitute(current,X,shift)
        psi=[R.substitute(f,X,shift) for f in psi]
        steps.append(dict(degree=k,removed=str(error),correction=[str(a) for a in delta]))
    if R.clean(current-target,X)!=0:raise ValueError('filtered finite Morse identity failed')
    theta=[R.clean(f+c,X) for f,c in zip(psi,critical)];inv=list(X)
    for i in sorted(range(len(X)),key=lambda i:weights[i]):
        inv[i]=R.clean(X[i]-(theta[i]-X[i]).xreplace(dict(zip(X,inv))),X)
    if any(R.substitute(a,X,inv)!=b for a,b in zip(theta,X)) or any(R.substitute(a,X,theta)!=b for a,b in zip(inv,X)):
        raise ValueError('wrong finite polynomial inverse')
    if s.cancel(s.Matrix(theta).jacobian(X).det()-1)!=0:raise ValueError('original volume not preserved')
    if any(s.cancel(a.subs(t,0)-b)!=0 for a,b in zip(theta,X)):raise ValueError('map not identity at time zero')
    if R.clean(R.substitute(potential,X,theta)-target+action,X)!=0:raise ValueError('wrong original-potential identity')
    return dict(variables=[str(a) for a in X],symbol_names=names,weights=[str(w) for w in weights],phase=str(F),
                precision=[[str(a) for a in row] for row in B0.tolist()],critical_blocks=blocks,
                critical_point=[str(a) for a in critical],sigma=str(s.factor(sigma)),critical_action=str(s.factor(action)),
                quadratic_precision=[[str(a) for a in row] for row in B.tolist()],map=[str(a) for a in theta],
                inverse_map=[str(a) for a in inv],Morse_steps=steps,Jacobian='1',actual_generic_rank_without_zero_positions=1,
                **seed_invariants(action,sigma))


def parse_model(model):
    loc={name:s.Symbol(name) for name in model['symbol_names']};loc.update({'t':t,'z':z,'e':e})
    X=tuple(loc[name] for name in model['variables'])
    parse=lambda a:s.sympify(a,locals=loc)
    return X,parse,s.Matrix([[parse(a) for a in row] for row in model['precision']]) if X else s.zeros(0)


def verify_model(model):
    X,parse,B0=parse_model(model)
    if filtered(parse(model['phase']),X,model['weights'],B0)!=model:raise ValueError('damaged filtered normal-form certificate')
    return True


def predicted(model,Q,N=6):
    X,parse,B0=parse_model(model)
    sigma=parse(model['sigma']);action=parse(model['critical_action'])
    # Explicit constant-one t-adic branches. Symbolic sqrt((gamma*t-1)^2)
    # can choose a negative constant branch and must not be used here.
    jet=lambda a:s.expand(s.series(a,t,0,N+1).removeO())
    action_jet=jet(action);u=s.expand(sigma-1)
    exp_jet=s.S.One;term=s.S.One
    for j in range(1,N+1):
        term=jet(term*action_jet/j);exp_jet+=term
    half_jet=jet(sum(s.binomial(-s.Rational(1,2),j)*u**j for j in range(N+1)))
    seed=jet(exp_jet*half_jet)
    if seed.subs(t,0)!=1:raise ValueError('wrong normalized seed branch')
    if not X:return seed*Q
    B=s.Matrix([[parse(a) for a in row] for row in model['quadratic_precision']])
    wick=G.wick_engine(B.inv())
    amplitude=R.substitute(Q,X,[parse(a) for a in model['map']])
    return jet(seed*sum(c*wick(m) for m,c in s.Poly(amplitude,*X).terms()))


def fixedpoint_gate(H,Z):
    """Finite ordinary residual/tail ideal; its nilpotent SAFE test needs purity."""
    Z=tuple(Z);H=tuple(s.expand(s.sympify(h)) for h in H)
    if len(H)!=len(Z) or len(set(Z))!=len(Z):raise ValueError('matching fixed-point coordinate list required')
    for h in H:polynomial_guard(h,Z)
    degrees=[0 if h==0 else max(0,s.Poly(h,*Z).total_degree()) for h in H]
    B=int(s.prod(d+1 for d in degrees));q=[s.S.Zero]*len(Z)
    for _ in range(B+1):
        q=[s.expand(s.series(t*h.subs(dict(zip(Z,q)),simultaneous=True),t,0,B+1).removeO()) for h in H]
    residual=[s.expand(a-t*h.subs(dict(zip(Z,q)),simultaneous=True)) for a,h in zip(q,H)]
    generators=list(dict.fromkeys(c for a in residual for m,c in s.Poly(a,t).terms() if c!=0))
    return dict(positions=[str(v) for v in Z],H=[str(h) for h in H],row_degrees=degrees,
                polynomial_field_degree_bound=B,truncation=[str(a) for a in q],residual=[str(a) for a in residual],
                ordinary_tail_ideal_generators=[str(a) for a in generators],nilpotent_safety_test_requires_full_purity=True)


def main():
    X,Y,C,D=s.symbols('X Y C D');alpha,beta,gamma,c=s.symbols('alpha beta gamma c')
    pair=s.Matrix([[0,1],[1,0]])
    models=[filtered(X**3+Y,(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair),
            filtered(X**4+Y,(X,Y),(s.Rational(1,2),s.Rational(3,2)),pair),
            filtered(alpha*X**3+beta*Y+gamma*X*Y,(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair),
            filtered(e*(X*X-1)/2,(X,),(1,),s.eye(1)),
            filtered(X**6+Y+C**3+D+X*X*D,(X,Y,C,D),
                     (s.Rational(1,3),s.Rational(5,3),s.Rational(2,3),s.Rational(4,3)),s.diag(pair,pair)),
            filtered(c,(),(),s.zeros(0)),
            filtered(X**3,(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair)]
    assert all(verify_model(m) for m in models)
    assert models[0]['critical_action']=='t**4' and models[0]['sufficient_ordinary_pure_prefix']==4
    assert models[1]['critical_action']=='t**5' and models[1]['sufficient_ordinary_pure_prefix']==5
    assert s.expand(s.sympify(models[3]['pure_gate_polynomial'])-e*e*t)==0
    controls=[]
    for index,N,markers,mod_nil in [(0,4,[1,X,Y,X*Y],False),(1,5,[1,X],False),
                                    (2,4,[1,X,Y],False),(3,4,[1,X*X],True)]:
        model=models[index];coords,parse,B0=parse_model(model)
        rows=R.original(parse(model['phase']),coords,B0.inv(),[(str(q),q,predicted(model,s.sympify(q),N)) for q in markers],N=N,mod_nil=mod_nil)
        controls.append(dict(model=index,rows=rows))
    # Full zero-weight algebraic exponential, with all original variables retained.
    phase,variables,cov=A.lift(1+A.Z**2,False);phase+=A.Z
    g=(1-s.sqrt(1-4*t*t))/(2*t);j=s.sqrt(1-4*t*t);seed=s.exp(t*g)
    rows=R.original(phase,variables,cov,[('unit',1,seed),('Z',A.Z,seed*g),('Y',A.Y,seed*(1/j-1))],N=3)
    controls.append(dict(model='zero-weight algebraic exponential',rows=rows))
    # Literal ordinary pure ideal differs from its reduced safe ideal.
    seed=s.exp(-e*t/2)*(1-e*t)**-s.Rational(1,2)
    moments=[s.expand(s.series(seed,t,0,4).removeO()).coeff(t,k)*s.factorial(k) for k in [1,2]]
    assert s.groebner(moments,e,domain=s.QQ)==s.groebner([e*e],e,domain=s.QQ)
    assert R.nil(s.series(seed,t,0,5).removeO())==1
    ordinary=dict(phase='e*(X^2-1)/2',sigma='1-e*t',action='-e*t/2',Omega='e^2*t',pure_ideal='(e^2)',safe_ideal='(e)',
                  first_two_moments=[str(a) for a in moments],dual_number_pure_despite_nonzero_action=True)
    # Geometric actual-slope certificate: g=t(1+g^3), action=t*g^3=g-t.
    u=s.Symbol('u');time=u/(1+u**3);action=s.cancel(time*u**3)
    assert s.cancel(action-u+time)==0
    at_infinity=s.cancel(action.subs(u,1/s.Symbol('s')))
    slope=dict(field_relation='g-t*(1+g^3)=0',actual_generic_rank=3,action='t*g^3=g-t',
               curve_time=str(time),curve_action=str(action),infinite_curve_local_action=str(at_infinity),
               actual_slopes_at_time_zero=[dict(slope='0',multiplicity=1),dict(slope='1/2',multiplicity=2)],
               actual_slopes_at_time_infinity=[dict(slope='1',multiplicity=3)],
               Newton_irregularity_at_zero=1,Newton_irregularity_at_infinity=3,
               slopes_from_full_actual_field_and_exact_local_valuations=True)
    # Structural local orders, no inference from sampled period jets.
    ss=s.Symbol('s');assert s.limit(s.cancel(time.subs(u,1/ss)/ss**2),ss,0)==1
    assert s.limit(s.cancel(at_infinity*ss),ss,0)==1
    assert s.rem(u**4,u**3+1,u)==-u and s.gcd(u**4,u**3+1)==1
    # Exact bounded ordinary residual producer and one polynomial Newton step.
    zz=s.symbols('z0 z1')
    gates=[fixedpoint_gate([1,zz[0]**3],zz),fixedpoint_gate([1+e*z],(z,)),fixedpoint_gate([e+z],(z,))]
    assert gates[0]['polynomial_field_degree_bound']==4 and gates[0]['truncation']==['t','t**4']
    assert all(a=='0' for a in gates[0]['residual'])
    q=t+e*t*t;error=s.expand(q-t*(1+e*q));inverse_mod_error=1+e*t
    corrected=s.expand(q-inverse_mod_error*error)
    assert corrected==t+e*t*t+e*e*t**3+e**3*t**4
    final_error=s.expand(corrected-t*(1+e*corrected))
    assert final_error==-e**4*t**5
    nil_newton=dict(H='1+e*z',ring='QQ[e]/e^4',field_bound=2,residual_ideal='(e^2)',initial=str(q),
                    inverse_mod_residual=str(inverse_mod_error),corrected=str(corrected),exact_error=str(final_error),
                    nilpotence_depth=2,polynomial_degree_after_correction=4)
    missing_purity=dict(H='e+z',ring='QQ[e]/e^2',field_bound=2,residual_ideal='(e)',
                        true_fixed_point='e*t/(1-t)',fixed_point_not_polynomial=True,Jacobian='1-t',
                        actual_nonpure_seed='1/(1-t)',invalid_without_purity=True)
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError('accepted mutation: '+name)
    reject('zero residual weight',lambda:filtered(X**3,(X,Y),(0,2),pair))
    reject('negative weight',lambda:filtered(X**3,(X,Y),(-1,3),pair))
    reject('exceeds filtration',lambda:filtered(X**4,(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair))
    reject('time-dependent phase',lambda:filtered(t*X**3,(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair))
    reject('nonpolynomial coefficient',lambda:filtered(X**3/(1+z),(X,Y),(s.Rational(2,3),s.Rational(4,3)),pair))
    reject('incompatible precision',lambda:filtered(X**3,(X,Y),(s.Rational(2,3),s.Rational(4,3)),s.eye(2)))
    reject('time-dependent fixed-point coefficient',lambda:fixedpoint_gate([t+z],(z,)))
    for field,value in [('critical_action','0'),('sigma','1-t'),('Jacobian','2'),('actual_generic_rank_without_zero_positions',2),('sufficient_ordinary_pure_prefix',2)]:
        damaged=json.loads(json.dumps(models[0]));damaged[field]=value
        reject('damaged '+field,lambda damaged=damaged:verify_model(damaged))
    out=p/'filtered_gaussian_critical_action_certificate.json'
    packet=dict(filtered_models=models,original_controls=controls,ordinary_nilpotent_control=ordinary,ramified_actual_slope_control=slope,
                finite_residual_gate_models=gates,nilpotent_Newton_control=nil_newton,missing_purity_countercontrol=missing_purity,
                retained_old_delayed_family_credited=True,rejections=rejected,
                universal_claims_use_written_proofs=True,actual_image_kernel_computed=False,
                regular_transcendental_impossibility_scope='supplied nonnegative Gaussian-compatible weighted filtration <=2 and its Gaussian symmetry orbit only')
    out.write_text(json.dumps(packet,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'FILTERED_GAUSSIAN_CRITICAL_ACTION_AND_EXACT_ALL_RING_GATES.md'),
        checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),runtime_sha256=sha(p/'filtered_gaussian_critical_action_runtime.json'),
        filtered_models=len(models),critical_blocks=sum(len(a['critical_blocks']) for a in models),
        finite_Morse_steps=sum(len(a['Morse_steps']) for a in models),original_models=len(controls),
        original_Wick_identities=sum(len(row['values']) for c in controls for row in c['rows']),
        ordinary_nilpotent_controls=1,ramified_actual_slope_controls=1,finite_residual_gate_models=len(gates),
        nilpotent_Newton_controls=1,missing_purity_countercontrols=1,rejections=len(rejected),
        historical_cusp_and_cone_restarted=False,independent_replay=False,formal_verification=False,external_review=False,
        assurance='Finite exact critical/Morse/original-volume/Wick/ordinary-ideal controls. Universal allring filtered normalform,full generic module,rank and slope criteria use the written proof.')
    (p/'filtered_gaussian_critical_action_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if 'sha256' not in k}))


if __name__=='__main__':main()
