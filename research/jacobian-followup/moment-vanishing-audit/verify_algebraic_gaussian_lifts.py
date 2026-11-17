"""Exact finite transcriptions for the all-ring algebraic Gaussian lift.

The canonical note proves universal full-image, normal-form, realization,
rank and monodromy assertions. These controls do not extrapolate them.
"""
from pathlib import Path
import hashlib, json, math
import sympy as s
import gaussian_generic_cutoff_producer as G

p=Path(__file__).resolve().parent
t,U,A,Z,W,Y,V,T,e,x=s.symbols('t U A Z W Y V T e x')
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()


def norm(f,nil=False):
    f=s.expand(f)
    return s.expand(f.coeff(e,0)+e*f.coeff(e,1)) if nil else f


def jet(f,N,nil=False):
    return norm(s.series(f,t,0,N+1).removeO(),nil)


def lift(H,with_time=True):
    H=s.expand(H)
    if not with_time and U in H.free_symbols:
        raise ValueError('time-dependent H needs the U,A pair')
    hz=s.diff(H,Z);C=(1+Y)*(4+6*Y+4*Y**2+Y**3)
    phase=s.expand((A if with_time else 0)+W*H+V*hz*(1+Y)-hz*C*T**2/2)
    variables=(U,A,Z,W,Y,V,T) if with_time else (Z,W,Y,V,T)
    cov=s.zeros(len(variables))
    for i in range(0,len(variables)-1,2):cov[i,i+1]=cov[i+1,i]=1
    cov[-1,-1]=1
    return phase,variables,cov


def fixed(H,N,nil=False):
    h=H.subs(U,t);g=s.S.Zero
    for _ in range(N+1):g=jet(t*h.subs(Z,g),N,nil)
    D=jet(1-t*s.diff(h,Z).subs(Z,g),N,nil)
    inv=jet(1/D,N,nil)
    return g,D,inv


def sparse_product(P,Q,nil=False):
    out={}
    for a,c in P.items():
        for b,d in Q.items():
            key=tuple(i+j for i,j in zip(a,b));out[key]=out.get(key,0)+c*d
    return {a:norm(c,nil) for a,c in out.items() if norm(c,nil)!=0}


def original_controls(H,with_time,N=3,nil=False):
    phase,variables,cov=lift(H,with_time);wick=G.wick_engine(cov)
    g,D,inv=fixed(H,N,nil);y=jet(inv-1,N,nil)
    h=H.subs(U,t);hz=s.diff(h,Z).subs(Z,g)
    hzz=s.diff(h,Z,2).subs(Z,g)
    # All momentum directions are checked, as well as the two linear markers.
    marked=[(1,s.S.One),(Z,g),(Y,y),(T**2,jet(D**4,N,nil)),
            (Z**2,jet(g**2,N,nil)),(Y**2,jet(y**2,N,nil)),
            (Z*W,inv),(Y*W,jet(t*hzz*inv**3,N,nil)),
            (V,jet(-2+s.Rational(1,2)*D**4*(4+6*y+4*y**2+y**3),N,nil))]
    if with_time:
        hu=s.diff(H,U).subs({U:t,Z:g})
        marked.extend([(U,t),(Z*A,jet(t*hu*inv,N,nil))])
    marker_maps=[(marker,{a:c for a,c in s.Poly(marker,*variables).terms()},predicted)
                 for marker,predicted in marked]
    phase_map={a:c for a,c in s.Poly(phase,*variables).terms()}
    power={(0,)*len(variables):s.S.One};rows=[]
    for m in range(N+1):
        row={'power':m,'moments':{}}
        for marker,mp,predicted in marker_maps:
            actual=norm(sum(c*d*wick(tuple(i+j for i,j in zip(a,b)))
                            for a,c in power.items() for b,d in mp.items()),nil)
            expected=norm(s.factorial(m)*s.expand(predicted).coeff(t,m),nil)
            if norm(actual-expected,nil)!=0:
                raise ValueError(f'original Wick/source identity mismatch: {H},{marker},{m},{actual},{expected}')
            row['moments'][str(marker)]=str(actual)
        rows.append(row)
        if m<N:power=sparse_product(power,phase_map,nil)
    return dict(H=str(H),dimension=len(variables),coefficient_ring='QQ[e]/e^2' if nil else 'QQ',
                phase=str(phase),fixed_point_jet=str(g),D_jet=str(D),inverse_D_jet=str(inv),
                original_Wick=rows,identities=len(rows)*len(marked))


def encode(F,a,k,v):
    """Finite specified-branch normalization; verifies every polynomial division."""
    if type(k)is not int or type(v)is not int or v<0 or k<=v:
        raise ValueError('need integers k>v>=0')
    F,a=s.expand(F),s.expand(a)
    if s.degree(a,t)>=k:raise ValueError('a must be the truncation below k')
    derivative=s.Poly(s.diff(F,x).subs(x,a),t)
    if derivative.is_zero or min(m[0] for m,c in derivative.terms())!=v:
        raise ValueError('wrong derivative valuation')
    u,w,z=s.symbols('u w z')
    raw=s.expand(F.subs(x,a+t**k*u))
    if any(m[0]<k+v for m,c in s.Poly(raw,t).terms()):
        raise ValueError('invalid branch truncation')
    normalized=s.cancel(raw/t**(k+v))
    s.Poly(normalized,t,u)
    lam=s.diff(normalized,u).subs(t,0)
    if not lam or u in lam.free_symbols:raise ValueError('not a simple normalized branch')
    u0=s.solve(normalized.subs(t,0),u)
    if len(u0)!=1:raise ValueError('branch constant is not unique')
    u0=u0[0]
    R=s.expand(w-normalized.subs(u,u0+w)/lam)
    c=R.subs(w,0).coeff(t,1)
    H=s.cancel((R.subs(w,t*(c+z))/t-c)/t)
    s.Poly(H,t,z)
    marker=s.expand(a+u0*t**k+c*t**(k+1)+t**(k+1)*z)
    relation=s.expand(z-t*H)
    if s.rem(F.subs(x,marker),relation,z)!=0:
        raise ValueError('failed exact algebraic recovery')
    return dict(F=str(F),truncation=str(a),k=k,derivative_valuation=v,
                normalized_polynomial=str(normalized),lambda_value=str(lam),u0=str(u0),c=str(c),
                H=str(s.expand(H)),original_recovery_marker=str(marker),
                full_generic_field='K(t,g)=K(t,f); written all-amplitude proof',
                dimension=7,actual_generic_rank=int(s.degree(F,x)))


def verify_encoding(model):
    symbols={'t':t,'x':x}
    if encode(s.sympify(model['F'],locals=symbols),s.sympify(model['truncation'],locals=symbols),
              model['k'],model['derivative_valuation'])!=model:
        raise ValueError('damaged branch encoding')
    return True


def rational_model(d):
    if type(d)is not int or d<2:raise ValueError('integer d>=2 required')
    z=s.Symbol('z');h=1+z**d;q=1-(d-1)*z**d
    assert s.expand(h-z*s.diff(h,z)-q)==0
    assert s.expand((d-1)*h+q)==d
    F=z-t*h
    return dict(d=d,H=str(h),q=str(q),implicit_polynomial=str(F),dimension=5,
                actual_rank=d,phase_degree=d+5,geometric_monodromy='S_'+str(d),
                global_image='K[z,1/h,1/q], t=z/h, D_t=(h^2/q)D_z',
                rational_core='K[t]' if d>=3 else 'K[t,1/(1-4*t^2)]',
                geometric_boundary_points={'origin':1,'simple_critical_values':d},
                geometric_length=d+3,characteristic_cycle={'zero':d,'origin':d-1,'each_critical':2},
                critical_value_polynomial=str(d**d*t**d-(d-1)**(d-1)),
                irreducibility_rank_monodromy_from_written_proof=True)


def verify_model(model):
    if rational_model(model['d'])!=model:raise ValueError('damaged all-rank algebraic image')
    return True


def main():
    a=s.Symbol('a');D=1-t*a;y=t*a/D
    C=lambda y:4+6*y+4*y**2+y**3
    precision=s.factor(1+t*a*(1+y)*C(y))
    assert s.cancel(precision-D**-4)==0
    assert s.cancel((1+y)*D-1)==0
    # Generic lower-triangular normal Jacobian, including cross derivatives.
    b,c,d=s.symbols('b c d')
    matrix=s.Matrix([[1,0,0],[-t*b,1-t*a,0],[-t*c*(1+Y),-t*d*(1+Y),1-t*a]])
    assert s.expand(matrix.det()-D**2)==0
    controls=[original_controls(1+Z**2,False),original_controls(1+U+Z**2,True),
              original_controls(Z**2,False),original_controls(1+e*Z,False,nil=True)]
    safe=[]
    H=U+(Z-U**2)**2*(Z+U)
    g,D0,iD=fixed(H,6)
    assert g==t**2 and D0==iD==1
    safe.append(dict(H=str(H),p=str(g),B='0',classification='safe field square remainder'))
    g,D0,iD=fixed(1+e*Z,5,True)
    assert g==t+e*t**2 and iD==1+e*t
    safe.append(dict(H='1+e*Z',ring='QQ[e]/e^2',p=str(g),B='e',inverse_D=str(iD),classification='safe nonreduced ordinary correction'))
    encodings=[encode(x**2-1-t,s.S.One,1,0),
               encode(x**2-t**2*(1+t),t,2,1),
               encode(x**2-1-t**3,s.S.One,1,0),
               encode(x-t*(1+x**5),s.S.Zero,1,0),
               encode((1-t)*x-1,s.S.One,1,0)]
    assert all(verify_encoding(m) for m in encodings)
    models=[rational_model(d) for d in [2,3,5]]
    assert all(verify_model(m) for m in models)
    coefficients=[]
    for d in [2,3,5]:
        for j in range(4):
            n=d*j+1;value=s.binomial(n,j)/n
            assert value!=0
            coefficients.append(dict(d=d,power=n,coefficient=str(value)))
    # Genuine original-marker trace failure, no guessed annihilator.
    z=s.Symbol('z');F=z-t*(1+z**3);q=1-2*z**3
    inv=s.invert(q,F,z)
    trace=s.factor(sum(s.expand(s.rem(inv**2*z**j,F,z)).coeff(z,j) for j in range(3)))
    assert trace==(9*t**3-4)/(27*t**3-4)
    assert s.cancel(q*(1+s.Rational(2,3)*((1+z**3)/q-1))-1)==0
    trace_control=dict(original_marker='(1+2*Y/3)^2',actual_period='1/(1-2*g^3)^2',
                       trace=str(trace),normalized_trace=str(trace/3),
                       fails_global_lattice='finite poles whereas full rational core is K[t]')
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError('accepted mutation: '+name)
    reject('omit time pair for time-dependent H',lambda:lift(1+U+Z,False))
    reject('wrong branch valuation',lambda:encode(x*x-t*t*(1+t),t,2,0))
    reject('insufficient k',lambda:encode(x*x-t*t*(1+t),t,1,1))
    reject('wrong branch truncation',lambda:encode(x*x-1-t,s.S.Zero,1,0))
    reject('degree one substituted in all-rank tower',lambda:rational_model(1))
    reject('boolean degree',lambda:rational_model(True))
    damaged=dict(encodings[1]);damaged['original_recovery_marker']='0'
    reject('drop actual branch recovery',lambda:verify_encoding(damaged))
    for field,value in [('actual_rank',1),('rational_core','K(t)'),('geometric_monodromy','trivial'),('geometric_length',1)]:
        damaged=dict(models[1]);damaged[field]=value
        reject('damaged '+field,lambda damaged=damaged:verify_model(damaged))
    out=p/'algebraic_gaussian_lift_certificate.json'
    packet=dict(compensator_precision=str(precision),normal_jacobian=str(matrix.det()),
                original_controls=controls,safe_normal_forms=safe,branch_encodings=encodings,
                rational_cover_models=models,all_rank_marker_coefficients=coefficients,
                genuine_Gaussian_trace_failure=trace_control,rejections=rejected,
                universal_claims_from_written_proof_not_controls=True,actual_image_kernel_computed=False)
    out.write_text(json.dumps(packet,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md'),
                 checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),
                 runtime_sha256=sha(p/'algebraic_gaussian_lift_runtime.json'),
                 original_models=len(controls),original_Wick_identities=sum(c['identities'] for c in controls),
                 safe_normal_forms=len(safe),branch_encodings=len(encodings),rational_cover_models=len(models),
                 all_rank_marker_coefficients=len(coefficients),genuine_trace_controls=1,rejections=len(rejected),
                 historical_cusp_and_cone_restarted=False,independent_replay=False,formal_verification=False,external_review=False,
                 assurance='Finite transcriptions; universal determinant cancellation,all-amplitude image,all-ring normalform,algebraic realization and all-d monodromy/lattice proofs are written internally.')
    (p/'algebraic_gaussian_lift_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if 'sha256' not in k}))


if __name__=='__main__':main()
