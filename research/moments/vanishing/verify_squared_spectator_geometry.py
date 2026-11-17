"""Exact transcriptions for squared-factor invariance and arbitrary-rank images.

Constructors accept every positive ell and any distinct nonzero rational
pole parameters. Small controls check formulas; the canonical note proves
all-ring, all-index, irreducibility and unbounded-rank statements.
"""
from pathlib import Path
import hashlib,json,math
import sympy as sp
import gaussian_generic_cutoff_producer as G
import gaussian_ordinary_cutoff_producer as O

p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()


def ore_product(P,Q):
    """Laurent t / Euler normal order: theta*t^b=t^b*(theta+b)."""
    out={}
    for (a,i),c in P.items():
        for (b,j),d in Q.items():
            for k in range(i+1):
                key=(a+b,k+j);out[key]=out.get(key,sp.S.Zero)+c*d*math.comb(i,k)*sp.Integer(b)**(i-k)
    return {e:sp.expand(c) for e,c in out.items() if c}


def initial(P,ell):
    if not P:raise ValueError('zero weighted operator')
    w=min(ell*a-i for (a,i),c in P.items() if c)
    return w,{e:c for e,c in P.items() if ell*e[0]-e[1]==w and c}


def coefficients(ell,m):
    k=ell*m;return sp.Integer(math.factorial(2*k)//(2**k*math.factorial(k)))


def constructor(ell,poles=(1,)):
    if type(ell) is not int or ell<1:raise ValueError('positive integer ell required')
    poles=tuple(sp.Rational(v) for v in poles)
    if any(not v for v in poles) or len(set(poles))!=len(poles):raise ValueError('distinct nonzero poles required')
    z,t=sp.symbols('z t');b=sp.Poly(sp.prod(2*ell*z+2*j+1 for j in range(ell)),z)
    lead=b.LC();bs={int(e[0]):c for e,c in b.terms()}
    matrix=sp.zeros(1+ell*len(poles));basis=['1'];operators=[]
    for h,lam in enumerate(poles):
        start=1+h*ell;basis += ['theta^'+str(j)+' h_'+str(lam) for j in range(ell)]
        A={(0,0):sp.S.One};A.update({(1,j):-lam*c for j,c in bs.items()})
        w,ini=initial(A,ell)
        if w!=0 or ini!={(0,0):sp.S.One,(1,ell):-lam*(2*ell)**ell}:raise ValueError('wrong prime weighted symbol')
        for j in range(ell-1):matrix[start+j+1,start+j]=1
        col=start+ell-1;matrix[0,col]=-1/(lam*t*lead)
        for j in range(ell):matrix[start+j,col]=-bs.get(j,0)/lead
        matrix[start,col]+=1/(lam*t*lead)
        operators.append(dict(lambda_value=str(lam),inhomogeneous_A=[[list(e),str(c)] for e,c in sorted(A.items())],
            homogeneous_theta_A=[[list(e),str(c)] for e,c in sorted(ore_product({(0,1):sp.S.One},A).items())],
            weighted_initial=[[list(e),str(c)] for e,c in sorted(ini.items())]))
    return dict(ell=ell,poles=[str(v) for v in poles],actual_rank=matrix.rows,quotient_rank=ell*len(poles),
        positive_slope=str(sp.Rational(1,ell)) if poles else None,positive_slope_multiplicity=ell*len(poles),Newton_irregularity=len(poles),
        theta_basis=basis,theta_connection_columns=[[str(v) for v in row] for row in matrix.tolist()],
        b_terms=[[j,str(c)] for j,c in sorted(bs.items())],operators=operators,
        global_image='sum D_t*h_lambda; empty pole set gives K[t]',
        exact_rank_uses_written_prime_symbol_and_nonisomorphism_proof=True,actual_image_kernel_computed=False)


def verify_constructor(model):
    ell=model['ell'];poles=tuple(sp.Rational(v) for v in model['poles'])
    if constructor(ell,poles)!=model:raise ValueError('damaged full-image model')
    t=sp.Symbol('t');K=7
    values=[sp.S.One]
    for lam in poles:
        for j in range(ell):values.append(sum(coefficients(ell,m)*lam**m*sp.Integer(m)**j*t**m for m in range(K+2)))
    matrix=sp.Matrix([[sp.sympify(c) for c in row] for row in model['theta_connection_columns']])
    for col,f in enumerate(values):
        rhs=sum(matrix[row,col]*values[row] for row in range(len(values)))
        remainder=sp.expand(t*sp.diff(f,t)-rhs)
        if any(remainder.coeff(t,m)!=0 for m in range(-1,K+1)):raise ValueError('wrong original connection action')
    T,xi=sp.symbols('T xi')
    for op,lam in zip(model['operators'],poles):
        F=1-lam*(2*ell)**ell*T*xi**ell
        if sp.cancel(F.subs(T,1/(lam*(2*ell)**ell*xi**ell)))!=0:raise ValueError('wrong prime-symbol quotient')
        if sp.rem(xi*(lam*(2*ell)**ell*T*xi**(ell-1))-1,F,T)!=0:raise ValueError('wrong inverse xi in prime quotient')
        for m in range(7):
            a=coefficients(ell,m);next_a=coefficients(ell,m+1)
            if next_a!=sp.prod(2*ell*m+2*j+1 for j in range(ell))*a:raise ValueError('wrong full-index recurrence transcription')
    for i,lam in enumerate(poles):
        for mu in poles[i+1:]:
            Fl=1-lam*(2*ell)**ell*T*xi**ell;Fm=1-mu*(2*ell)**ell*T*xi**ell
            if sp.expand(mu/(mu-lam)*Fl-lam/(mu-lam)*Fm)!=1:raise ValueError('wrong nonisomorphism Bezout certificate')
    return True


def pole_marker_operator(ell,pole_order,half_marker):
    if type(pole_order) is not int or pole_order<1 or type(half_marker) is not int or half_marker<0:raise ValueError('invalid pole/marker index')
    z=sp.Symbol('z')
    return sp.Poly(sp.Rational(1,math.factorial(pole_order-1))*sp.prod(z+j for j in range(1,pole_order))*
                   sp.prod(2*ell*z+2*j+1 for j in range(half_marker)),z)


def phase_data(ell,deformation=False):
    constructor(ell)
    # Z,W,S,Y, followed by parameter e when requested.
    terms=[([0,1,0,2*ell],0,sp.S.One),([1,1,0,2*ell],1,sp.S.One),
           ([0,0,2,2*ell],1,-sp.S.One),([1,0,2,2*ell],2,-sp.Rational(3,2)),
           ([2,0,2,2*ell],3,-sp.Rational(1,2))]
    return dict(dimension=4,parameter_count=1,
        phase_terms=[[e+[power if deformation else 0],str(c)] for e,power,c in terms],
        precision_matrix=[[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])


def main():
    models=[]
    for ell,poles in [(1,(1,)),(2,(1,)),(5,(1,)),(2,(1,2))]:
        model=constructor(ell,poles);assert verify_constructor(model);models.append(model)
    # Commutators increase weight; arbitrary finite Laurent examples test order.
    P={(0,2):sp.S.One,(1,0):sp.Integer(2)};Q={(-1,1):sp.Integer(3),(0,0):sp.S.One}
    multiplicative=[]
    for ell in [1,2,5]:
        wP,inP=initial(P,ell);wQ,inQ=initial(Q,ell);wPQ,inPQ=initial(ore_product(P,Q),ell)
        commutative={}
        for (a,i),c in inP.items():
            for (b,j),d in inQ.items():commutative[(a+b,i+j)]=commutative.get((a+b,i+j),0)+c*d
        assert wPQ==wP+wQ and inPQ==commutative
        multiplicative.append(dict(ell=ell,weight=wPQ,initial_terms=[[list(e),str(c)] for e,c in inPQ.items()]))
    kernel=[];z=sp.Symbol('z')
    for ell,pole,k in [(1,1,0),(2,3,2),(5,2,1)]:
        op=pole_marker_operator(ell,pole,k)
        for m in range(6):
            left=sp.binomial(m+pole-1,pole-1)*sp.factorial2(2*ell*m+2*k-1)
            assert left==op.eval(m)*coefficients(ell,m)
        kernel.append(dict(ell=ell,pole_order=pole,half_marker=k,theta_polynomial=str(op.as_expr())))
    actual=[]
    for ell in [1,2,5]:
        data=phase_data(ell,True);e=sp.Symbol('c0')
        for m in range(4):
            pure=O.original_moment(data,(0,0,0,0),m)
            marked=O.original_moment(data,(1,0,0,0),m)
            assert pure==(1 if m==0 else 0)
            assert marked==(0 if m==0 else sp.factorial(m)*coefficients(ell,m)*e**(m-1))
            actual.append(dict(ell=ell,power=m,pure=str(pure),Z_marker=str(marked)))
    # Unit-valued rational moments preserve ordinary ideals, not just radicals.
    a,b,X,Y=sp.symbols('a b X Y');base=a+b*X;profile=[]
    wick=G.wick_engine(sp.eye(1))
    def expectation(expr,var):return sp.expand(sum(c*wick((e[0],)) for e,c in sp.Poly(sp.expand(expr),var).terms()))
    for square in [Y**2,(1+Y)**2]:
        for d,m in [(0,1),(1,2),(2,3)]:
            old=[expectation(X**i*base**m,X) for i in range(d+1)]
            new=[expectation(X**i*base**m,X)*expectation(Y**j*square**m,Y) for i in range(d+1) for j in range(d+1-i)]
            assert sp.groebner(old,a,b,domain=sp.QQ)==sp.groebner(new,a,b,domain=sp.QQ)
            assert expectation(square**m,Y)>0
            constraint=[a*a*b*b]
            assert sp.groebner(old+constraint,a,b,domain=sp.QQ)==sp.groebner(new+constraint,a,b,domain=sp.QQ)
            profile.append(dict(square=str(square),degree=d,power=m,beta_zero_unit=str(expectation(square**m,Y))))
    # All-N,E formula is proved by the retained charge identity, checked here
    # on its FIVE defining monomials rather than an expanding moment prefix.
    data=phase_data(2,True)
    assert all(es[-1]==es[0]-es[1]+1 for es,c in data['phase_terms'])
    Z,W,S,V,e=sp.symbols('Z W S Y e');phase=O.polynomial(data['phase_terms'],(Z,W,S,V,e))
    charged=sp.expand(Z*sp.diff(phase,Z)-W*sp.diff(phase,W));assert charged!=0
    radial_polynomials=[(Z*W)**2,(Z*W)+S*S+V*V]
    assert all(sp.expand(Z*sp.diff(f,Z)-W*sp.diff(f,W))==0 for f in radial_polynomials)
    c,delta,offset=sp.symbols('c delta offset')
    assert sp.expand((c-1)*delta+offset-(c+offset-1))==sp.expand((c-1)*(delta-1))
    sharp=[]
    for L,C,E in [(1,1,0),(3,1,2),(3,2,0),(3,sp.Rational(5,2),sp.Rational(1,2)),(7,4,3),(7,sp.Rational(1,2),8)]:
        C,E=sp.Rational(C),sp.Rational(E)
        if C<1:
            m=int(sp.floor(E/(1-C)))+1;marker=m;exponent=0;valid=False
        else:
            m=int(sp.floor(C+E))+1;marker=1;exponent=m-1;valid=exponent>=L
        assert m>C*marker+E
        assert valid==(C>=1 and C+E>=L)
        sharp.append(dict(nilpotence_depth=L,slope=str(C),intercept=str(E),valid=valid,
            boundary_marker_Z_degree=marker,boundary_power=m,coefficient_e_exponent=exponent,
            coefficient_nonzero_mod_eL=exponent<L))
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError('accepted invalid geometry: '+name)
    for value in [0,-1,True,sp.Rational(1,2)]:reject('invalid ell '+str(value),lambda:constructor(value))
    for poles in [(0,),(1,1)]:reject('invalid pole packet '+str(poles),lambda:constructor(1,poles))
    damaged=json.loads(json.dumps(models[1]));damaged['actual_rank']=1
    reject('source rank one substituted for actual rank',lambda:verify_constructor(damaged))
    damaged=json.loads(json.dumps(models[-1]));damaged['theta_connection_columns'][0][-1]='0'
    reject('split unit extension substituted',lambda:verify_constructor(damaged))
    damaged=json.loads(json.dumps(models[0]));damaged['positive_slope']='0'
    reject('regularity substituted for irregular quotient',lambda:verify_constructor(damaged))
    out=p/'squared_spectator_geometry_certificate.json'
    certificate=dict(models=models,weighted_products=multiplicative,pole_marker_operators=kernel,
        original_Wick_controls=actual,ordinary_profile_controls=profile,sharp_affine_controls=sharp,
        nonradial_charge_derivative=str(charged),rejections=rejected,
        all_index_claims_from_written_proof_not_finite_extrapolation=True,
        arbitrary_nilpotent_depth_and_arbitrary_rank_from_formulas=True,
        global_quotient_semisimplicity_claimed=False,actual_image_kernel_computed=False)
    out.write_text(json.dumps(certificate,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'SQUARED_SPECTATOR_INVARIANCE_AND_UNBOUNDED_GAUSSIAN_IMAGES.md'),
        checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),
        runtime_sha256=sha(p/'squared_spectator_geometry_runtime.json'),
        imports_sha256={name:sha(p/name) for name in ['gaussian_generic_cutoff_producer.py','gaussian_ordinary_cutoff_producer.py',
            'FINITE_MULTIPLIER_MIXED_DECISION.md','RANK_ONE_IMAGES_AND_NILPOTENT_CUTOFFS.md']},
        model_controls=len(models),weighted_product_controls=len(multiplicative),pole_marker_controls=len(kernel),
        original_Wick_controls=len(actual),ordinary_profile_controls=len(profile),sharp_affine_controls=len(sharp),rejections=len(rejected),
        historical_cusp_and_cone_restarted=False,independent_replay=False,formal_verification=False,external_review=False,
        assurance='Finite exact formula transcriptions; general rank irreducibility, all-ring unit transfer and all-real affine geometry use the complete written proofs. No universal assurance promotion.')
    (p/'squared_spectator_geometry_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k.endswith('controls') or k in ['completed','rejections']}))


if __name__=='__main__':main()
