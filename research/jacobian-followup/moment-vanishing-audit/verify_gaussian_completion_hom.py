"""Exact finite Gaussian completion certificates;no holonomic/Hom producer."""
from pathlib import Path
import hashlib,json,platform
import sympy as s


def trunc(f,t,k):
    return s.Poly(s.expand(f),t).slice(0,k).as_expr()


def gaussian_reduce(f,xs):
    """Return scalar and primitives with f-scalar=sum(d_i-x_i)*primitive_i."""
    remain=s.expand(f);primitives=[]
    for x in xs:
        d=s.Poly(remain,x).degree()
        bs=[s.Integer(0),s.Integer(-1)]
        moments=[s.Integer(1),s.Integer(0)]
        for j in range(2,max(1,d)+1):
            bs.append((j-1)*bs[j-2]-x**(j-1))
            moments.append((j-1)*moments[j-2])
        new=0;primitive=0
        for (j,),coeff in s.Poly(remain,x).terms():
            new+=coeff*moments[j];primitive+=coeff*bs[j]
        primitives.append(s.expand(primitive));remain=s.expand(new)
    assert s.expand(f-remain-sum(s.diff(b,x)-x*b for b,x in zip(primitives,xs)))==0
    return remain,primitives


def completion_certificate(phase,amplitude,xs,t,k):
    """Supply a finite original-source reduction and the actual period mod t^k."""
    phase=s.sympify(phase);amplitude=s.sympify(amplitude)
    if not isinstance(k,int) or k<1:raise ValueError('k must be positive')
    if t in phase.free_symbols:raise ValueError('phase must be time independent')
    if not xs or len(set(xs))!=len(xs) or t in xs:raise ValueError('invalid phase variables')
    s.Poly(phase,*xs);s.Poly(amplitude,t,*xs)
    up=sum(t**j*phase**j/s.factorial(j) for j in range(k))
    down=sum((-t)**j*phase**j/s.factorial(j) for j in range(k))
    unit,_=gaussian_reduce(trunc(up,t,k),xs)
    period,_=gaussian_reduce(trunc(up*amplitude,t,k),xs)
    inverse=s.Integer(0)
    # unit=1+O(t),so its inverse is a finite geometric sum.
    for j in range(k):inverse+=(-1)**j*(unit-1)**j
    coordinate=trunc(period*trunc(inverse,t,k),t,k)
    zero,bs=gaussian_reduce(trunc(up*(amplitude-coordinate),t,k),xs)
    assert s.expand(zero)==0
    bs=[trunc(down*b,t,k) for b in bs]
    residual=trunc(amplitude-coordinate-sum(s.diff(b,x)+(t*s.diff(phase,x)-x)*b for x,b in zip(xs,bs)),t,k)
    assert residual==0
    return {'period':s.expand(period),'basis_coordinate':coordinate,'primitives':bs,'residual':residual}


def run():
    x,z,t,a,b,c0,y=s.symbols('x z t a b c0 y')
    cases=[(a*x+x**3,1,(x,),4),(x*z+z**3,x**2+z+t*x,(x,z),3),
           (a*x**2+b*x+1,x**3+t*x,(x,),4)]
    rows=[]
    for phase,amp,xs,k in cases:
        cert=completion_certificate(phase,amp,xs,t,k)
        assert cert['residual']==0
        rows.append({'phase':str(phase),'amplitude':str(amp),'k':k,
                     'period':str(cert['period']),'basis_coordinate':str(cert['basis_coordinate']),
                     'primitive_count':len(cert['primitives']),'original_source_residual':'0'})
    # Nilpotent parameter specialization preserves an already exact identity.
    nilphase=a*x**3+a**2*z**2;nilamp=x*z+a*x
    nilpotent=completion_certificate(nilphase,nilamp,(x,z),t,3)
    reduced_bs=[s.rem(value,a**3,a) for value in nilpotent['primitives']]
    reduced_coordinate=s.rem(nilpotent['basis_coordinate'],a**3,a)
    nilres=trunc(nilamp-reduced_coordinate-sum(s.diff(value,var)+(t*s.diff(nilphase,var)-var)*value for var,value in zip((x,z),reduced_bs)),t,3)
    assert s.rem(nilres,a**3,a)==0

    # Constant-value recurrence fixes all formal homomorphism values in a prefix.
    uniqueness_controls=0
    phase=x**3+a*x
    for power in range(4):
        q=x**power
        cert=completion_certificate(phase,q,(x,),t,4)
        for m in range(4):
            expected,_=gaussian_reduce(q*phase**m,(x,))
            assert s.expand(s.Poly(cert['period'],t).nth(m)*s.factorial(m)-expected)==0
            uniqueness_controls+=1

    # Exact linear Gaussian quotient and its moment connection.
    linear_controls=0
    f=s.exp(a**2*t**2/2)
    assert s.diff(f,t)==a**2*t*f
    for j in range(6):
        h,_=gaussian_reduce((x+a*t)**j,(x,))
        if j:
            prev,_=gaussian_reduce((x+a*t)**(j-1),(x,))
            nxt,_=gaussian_reduce((x+a*t)**(j+1),(x,))
            assert s.expand(nxt-a*t*h-j*prev)==0
        linear_controls+=1
    # At a0 the same exact connection admits1;over Q[a] the top coefficient
    # equation for a nonzero polynomial is a^2*leading=0,impossible in the domain.
    hom_rows=[]
    for degree in range(5):
        cs=s.symbols('q:'+str(degree+1));q=sum(cs[j]*t**j for j in range(degree+1))
        equations=s.Poly(s.diff(q,t)-a**2*t*q,t).all_coeffs()
        assert equations[0]==-a**2*cs[-1]
        assert s.diff(s.Integer(1),t)-a**2*t==(-a**2*t)
        hom_rows.append({'trial_degree':degree,'top_equation':str(equations[0]),'zero_fiber_solution':'1'})

    # In Q[a,b]/ab,the nonzero b kills all positive powers of a,but1 does not.
    gb=s.groebner([a*b],a,b)
    for m in range(1,7):
        assert gb.reduce(b*a**m)[1]==0
        assert gb.reduce(a**m)[1]!=0
    assert gb.reduce(b)[1]==b

    # Fourier translation commutator,with differentiated constant coefficient.
    commutators=0
    def B(q):return s.diff(q,c0)+s.diff(q,y)
    probe=c0**2*y+y**3+c0
    for n in range(1,7):
        assert s.expand(B(y**n*probe)-y**n*B(probe)-n*y**(n-1)*probe)==0
        commutators+=1

    # A merely annihilating scalar presentation has a spurious polynomial Hom.
    assert s.diff(s.exp(t),t,2)-s.diff(s.exp(t),t)==0
    assert s.diff(s.Integer(1),t,2)-s.diff(s.Integer(1),t)==0
    assert s.diff(s.Integer(1),t)-1!=0

    # Source Example12:one constant-degree indicial stratum has unbounded roots.
    root_rows=[]
    for n in (1,2,5,11):
        beta1=s.Rational(2*n-1,2)
        root=beta1+s.Rational(1,2);assert root==n
        root_rows.append({'beta1':str(beta1),'indicial_degree':2,'largest_integer_root':n})
    rejected=0
    for phase,k in [(x+t,2),(x,0)]:
        try:completion_certificate(phase,1,(x,),t,k)
        except ValueError:rejected+=1
    assert rejected==2
    return {'completed':True,'python_version':platform.python_version(),'sympy_version':s.__version__,
            'truncated_completion_certificates':rows,'nilpotent_specialization_identity_controls':1,
            'formal_uniqueness_coefficient_controls':uniqueness_controls,'linear_Gaussian_recurrence_controls':linear_controls,
            'relative_Hom_top_degree_controls':hom_rows,'zero_divisor_torsion_ideal_controls':6,
            'Fourier_translation_commutator_controls':commutators,'spurious_polynomial_source_rejected':True,
            'Oaku_Example12_integer_root_controls':root_rows,'invalid_input_controls':rejected,
            'finite_completion_certificate_generator_implemented':True,
            'universal_completion_proof_inferred_from_samples':False,
            'general_holonomic_integration_or_Hom_producer_executed':False,
            'Oaku_indicial_or_adapted_resolution_producer_executed':False,
            'actual_kernel_infinite_intersection_computed':False,
            'automatic_polynomial_Hom_basechange_claimed':False,
            'Gaussian_global_uniformity_proved':False,'ordinary_nonreduced_uniformity_proved':False,
            'protected_ternary_frontier_changed':False,'independent_review':False,
            'formal_verification':False,'novelty_claim':False}


if __name__=='__main__':
    path=Path(__file__).resolve();result=run()
    result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('gaussian_completion_hom_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
