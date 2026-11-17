"""Finite weighted polynomial Morse producer and exact Gaussian controls.

Universal arbitrary-ring/image/rank/safety claims use the canonical proof.
This producer works in rational function coefficients, returning maps,
inverse maps, determinant and original-potential identities for review.
"""
from pathlib import Path
import hashlib,json,math
import sympy as s
import gaussian_generic_cutoff_producer as G

p=Path(__file__).resolve().parent
t,z,e=s.symbols('t z e')
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()


def clean(f,X):
    return s.Add(*(s.cancel(c)*s.prod(v**a for v,a in zip(X,m)) for m,c in s.Poly(s.expand(f),*X).terms()))


def substitute(f,X,values):
    return clean(f.xreplace(dict(zip(X,values))),X)


def degree_part(f,X,k):
    return s.Add(*(c*s.prod(v**a for v,a in zip(X,m)) for m,c in s.Poly(f,*X).terms() if sum(m)==k))


def morse(F,X,weights,B0):
    weights=tuple(s.Rational(w) for w in weights)
    if len(X)!=len(weights) or len(set(X))!=len(X) or not X:
        raise ValueError('a nonempty distinct residual coordinate list is required')
    if any(w<=0 or w>=2 for w in weights):raise ValueError('residual weights must lie strictly between zero and two')
    if B0.shape!=(len(X),len(X)) or B0!=B0.T or not B0.det():raise ValueError('nonsingular symmetric precision required')
    if any(v.is_Rational is not True for v in B0):raise ValueError('fixed rational Gaussian precision required')
    for i in range(len(X)):
        for j in range(len(X)):
            if B0[i,j] and weights[i]+weights[j]!=2:raise ValueError('precision is not weight-compatible')
    F=s.expand(F)
    if t in F.free_symbols or t in X:raise ValueError('phase coefficients must be independent of the time variable')
    try:
        s.Poly(F,*sorted(F.free_symbols|set(X),key=str))
    except s.PolynomialError as exc:
        raise ValueError('polynomial phase coefficients required') from exc
    for m,c in s.Poly(F,*X).terms():
        if c and (sum(a*w for a,w in zip(m,weights))!=2 or sum(m)<2):
            raise ValueError('phase residual is not homogeneous of weight two')
    q=(s.Matrix(X).T*B0*s.Matrix(X))[0]/2
    S=s.expand(q-t*F)
    B=s.hessian(S,X).subs(dict.fromkeys(X,0));inverse=B.inv()
    target=(s.Matrix(X).T*B*s.Matrix(X))[0]/2
    sigma=s.factor(B.det()/B0.det())
    current=S;psi=list(X);steps=[]
    for k in range(3,int(s.floor(2/min(weights)))+1):
        error=degree_part(current,X,k)
        if error==0:continue
        delta=[s.cancel(v) for v in -inverse*s.Matrix([s.diff(error,v) for v in X])/k]
        shift=[v+d for v,d in zip(X,delta)]
        for i,value in enumerate(delta):
            if any(variable in value.free_symbols and weights[j]>=weights[i] for j,variable in enumerate(X)):
                raise ValueError('correction is not strictly weight-triangular')
        current=substitute(current,X,shift)
        psi=[substitute(v,X,shift) for v in psi]
        steps.append(dict(degree=k,removed=str(error),correction=[str(v) for v in delta]))
    if clean(current-target,X)!=0:raise ValueError('finite Morse reduction failed')
    inv=list(X)
    for i in sorted(range(len(X)),key=lambda i:weights[i]):
        nonlinear=psi[i]-X[i]
        inv[i]=clean(X[i]-nonlinear.xreplace(dict(zip(X,inv))),X)
    if any(substitute(value,X,inv)!=variable for value,variable in zip(psi,X)):
        raise ValueError('wrong polynomial inverse')
    if any(substitute(value,X,psi)!=variable for value,variable in zip(inv,X)):
        raise ValueError('wrong reverse polynomial inverse')
    if any(s.cancel(value.subs(t,0)-variable)!=0 for value,variable in zip(psi,X)):
        raise ValueError('coordinate change is not identity at time zero')
    if s.cancel(s.det(s.Matrix(psi).jacobian(X))-1)!=0:raise ValueError('wrong original volume')
    return dict(variables=[str(v) for v in X],weights=[str(w) for w in weights],residual_phase=str(F),
                precision=[[str(v) for v in row] for row in B0.tolist()],sigma=str(sigma),
                quadratic_precision=[[str(v) for v in row] for row in B.tolist()],
                map=[str(v) for v in psi],inverse_map=[str(v) for v in inv],steps=steps,
                determinant='1',coefficient_denominators='only powers of sigma and rational units')


def verify_model(model):
    X=s.symbols(' '.join(model['variables']),seq=True);loc={str(v):v for v in X};loc.update({'t':t,'z':z,'e':e})
    F=s.sympify(model['residual_phase'],locals=loc)
    B0=s.Matrix([[s.sympify(v,locals=loc) for v in row] for row in model['precision']])
    expected=morse(F,X,model['weights'],B0)
    if model!=expected:raise ValueError('damaged weighted Morse certificate')
    return True


def nil(f):
    return s.expand(s.series(f,e,0,2).removeO())


def original(F,X,cov,markers,N=4,mod_nil=False):
    wick=G.wick_engine(cov);phase={a:c for a,c in s.Poly(F,*X).terms()}
    power={(0,)*len(X):s.S.One};rows=[]
    for m in range(N+1):
        row={'power':m,'values':{}}
        for label,marker,predicted in markers:
            mp=s.Poly(marker,*X).terms()
            value=s.expand(sum(c*d*wick(tuple(i+j for i,j in zip(a,b))) for a,c in power.items() for b,d in mp))
            expected=s.factorial(m)*s.expand(s.series(predicted,t,0,N+1).removeO()).coeff(t,m)
            diff=s.expand(value-expected)
            if mod_nil:diff=nil(diff);value=nil(value)
            if diff!=0:raise ValueError(f'wrong original Wick identity {label},m={m}: {diff}')
            row['values'][label]=str(value)
        rows.append(row)
        nxt={}
        if m<N:
            for a,c in power.items():
                for b,d in phase.items():
                    key=tuple(i+j for i,j in zip(a,b));nxt[key]=nxt.get(key,0)+c*d
            power={a:(nil(c) if mod_nil else s.expand(c)) for a,c in nxt.items() if c!=0}
    return rows


def cutoff(r,D,L):
    if any(type(v)is not int for v in [r,D,L]) or r<0 or D<0 or L<1:raise ValueError('nonnegative size/degree and positive nilpotence exponent required')
    if r==0:return dict(slope=0,intercept=0)
    return dict(slope=int(s.ceiling(1+s.Rational((max(2,D)+1)*(r*L-1),2))),intercept=r*(L-1))


def main():
    a,b,c,d=s.symbols('a b c d');pair=s.Matrix([[0,1],[1,0]])
    models=[morse(a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),pair),
            morse(a**4,(a,b),(s.Rational(1,2),s.Rational(3,2)),pair),
            morse(a*b+2*c*d+a*a*d+c**3+a**6,(a,b,c,d),
                  (s.Rational(1,3),s.Rational(5,3),s.Rational(2,3),s.Rational(4,3)),s.diag(pair,pair)),
            morse(-(1+z)*a*b+z*a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),pair),
            morse(e*a*b+a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),pair)]
    assert all(verify_model(m) for m in models)
    controls=[]
    controls.append(dict(name='positive pure cubic',rows=original(a**3,(a,b),pair,
        [('unit',1,s.S.One),('b^3',b**3,6*t),('b^6',b**6,360*t*t)])))
    controls.append(dict(name='positive nonpure nilpotent',ring='QQ[e]/e^2',rows=original(e*a*b+a**3,(a,b),pair,
        [('unit',1,1+e*t),('a*b',a*b,1+2*e*t),('b^3',b**3,6*t+24*e*t*t)],mod_nil=True)))
    W=s.Symbol('W');phase=W*(1+z)-(1+z)*a*b+z*a**3;cov=s.diag(pair,pair)
    controls.append(dict(name='zero-weight unsafe pure',rows=original(phase,(z,W,a,b),cov,
        [('unit',1,s.S.One),('z',z,t/(1-t)),('a*b',a*b,1-t),('b^3',b**3,6*t*t*(1-t)**2)])))
    phase=W*z*z+z*a**3
    controls.append(dict(name='zero-weight safe pure',rows=original(phase,(z,W,a,b),cov,
        [('unit',1,s.S.One),('z',z,s.S.Zero),('W',W,2*t),('a*b',a*b,s.S.One)])))
    # Exact determinant purity and fixed-position gate, not a moment prefix.
    g=t/(1-t);J=1-t;sigma=(1+t*(1+g))**2
    assert s.cancel(J*J*sigma-1)==0
    assert s.cancel(1/(J*(1+t*(1+g)))-1)==0
    # Nonpure algebraic function-field degree must not replace differential rank.
    nonpure=dict(phase='T^2/2',fixed_point_field='K(t)',seed='(1-t)^(-1/2)',
                 actual_generic_rank=1,period_function_field_degree=2,
                 exact_equation='2*(1-t)*D_t(f)-f=0')
    f=(1-t)**-s.Rational(1,2)
    assert s.simplify(2*(1-t)*s.diff(f,t)-f)==0
    # Literal ordinary pure ideals: reciprocal square is triangular over QQ.
    u,v=s.symbols('u v');sig=1+u*t+v*t*t;f=s.series(sig**-s.Rational(1,2),t,0,5).removeO()
    first=[s.expand(f).coeff(t,j)*s.factorial(j) for j in [1,2]]
    assert s.groebner(first,u,v,domain=s.QQ)==s.groebner([u,v],u,v,domain=s.QQ)
    pure_ideal=dict(sigma=str(sig),first_moments=[str(a) for a in first],ordinary_ideal='(u,v)')
    # Nilpotent polynomial units and finite-binomial cutoff bounds.
    assert nil((1+e*t)*(1-e*t))==1
    bounds=[dict(r=r,D=D,L=L,**cutoff(r,D,L)) for r,D,L in [(0,0,1),(2,3,1),(4,6,1),(2,3,2),(5,9,4)]]
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError('accepted mutation: '+name)
    reject('zero residual weight',lambda:morse(a**3,(a,b),(0,2),pair))
    reject('negative residual weight',lambda:morse(a**3,(a,b),(-1,3),pair))
    reject('wrong weighted degree',lambda:morse(a**4,(a,b),(s.Rational(2,3),s.Rational(4,3)),pair))
    reject('incompatible Gaussian precision',lambda:morse(a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),s.eye(2)))
    reject('singular precision',lambda:morse(a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),s.zeros(2)))
    reject('time-dependent phase input',lambda:morse(t*a**3,(a,b),(s.Rational(2,3),s.Rational(4,3)),pair))
    reject('rational nonpolynomial phase input',lambda:morse(a**3/(1+z),(a,b),(s.Rational(2,3),s.Rational(4,3)),pair))
    for field,value in [('determinant','-1'),('map',['a','b']),('inverse_map',['a','b'])]:
        damaged=json.loads(json.dumps(models[0]));damaged[field]=value
        reject('damaged '+field,lambda damaged=damaged:verify_model(damaged))
    reject('nonpositive nilpotence depth',lambda:cutoff(2,3,0))
    reject('boolean dimension',lambda:cutoff(True,3,1))
    out=p/'semipositive_gaussian_reduction_certificate.json'
    packet=dict(Morse_models=models,original_controls=controls,nonpure_rank_boundary=nonpure,ordinary_pure_ideal=pure_ideal,
                cutoff_models=bounds,zero_weight_gate=dict(g=str(g),J=str(J),sigma=str(s.cancel(sigma)),pure=True,safe=False,bad_marker='z'),
                rejections=rejected,universal_claims_from_written_proofs=True,actual_image_kernel_computed=False,
                global_lattice_equality_under_Morse_reduction_claimed=False)
    out.write_text(json.dumps(packet,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'SEMIPOSITIVE_WEIGHTED_GAUSSIAN_REDUCTION_AND_SAFE_NORMAL_FORMS.md'),
        checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),runtime_sha256=sha(p/'semipositive_gaussian_reduction_runtime.json'),
        finite_Morse_models=len(models),nonlinear_Morse_steps=sum(len(m['steps']) for m in models),
        original_models=len(controls),original_Wick_identities=sum(len(row['values']) for control in controls for row in control['rows']),
        ordinary_pure_ideal_controls=1,nonpure_rank_controls=1,cutoff_models=len(bounds),rejections=len(rejected),
        historical_cusp_and_cone_restarted=False,independent_replay=False,formal_verification=False,external_review=False,
        assurance='Finite original identities and producer certificates;universal weighted triangular finiteness,allring safety/normalization and generic image/rank use the written proof.')
    (p/'semipositive_gaussian_reduction_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if 'sha256' not in k}))


if __name__=='__main__':main()
