"""Exact original Gaussian time-change controls; universal claims use the proof."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,sys
import sympy as s
p=Path(__file__).resolve().parent
sys.path.insert(0,str(p))
import gaussian_generic_cutoff_producer as G
t,z,u,A,x,Z,W,eta=s.symbols('t z u A x Z W eta')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
pair=s.Matrix([[0,1],[1,0]])

def expectation(poly,X,cov):
    poly=s.expand(poly)
    if not X:return poly
    wick=G.wick_engine(cov)
    return s.expand(sum(c*wick(m) for m,c in s.Poly(poly,*X).terms()))

def clock_data(B,X,cov,a):
    X=tuple(X);a=s.sympify(a);B=s.sympify(B)
    if not X or len(set(X))!=len(X) or any(not isinstance(v,s.Symbol) for v in X):raise ValueError('distinct original Gaussian variables required')
    if any(v in X for v in [u,A,t,z]):raise ValueError('reserved clock variable')
    if cov.shape!=(len(X),len(X)) or cov!=cov.T or cov.det()==0:raise ValueError('nonsingular symmetric covariance required')
    if any(not c.is_Rational for c in cov):raise ValueError('fixed rational covariance required')
    if s.denom(s.cancel(B))!=1 or any(v in B.free_symbols for v in [u,A,t,z]):raise ValueError('original polynomial phase required')
    s.Poly(B,*X)
    if s.denom(s.cancel(a))!=1:raise ValueError('polynomial clock required')
    s.Poly(a,t)
    if a==0 or a.subs(t,0)!=0:raise ValueError('nonconstant zero-origin clock required')
    if any(v in a.free_symbols for v in X+(u,A,z)):raise ValueError('clock independent of original variables required')
    return A+s.cancel(a.subs(t,u)/u)*B,(u,A)+X,s.diag(pair,cov)

def direct_moment(B,X,cov,a,alpha,beta,Q,m):
    phase,Y,V=clock_data(B,X,cov,a)
    return expectation(u**alpha*A**beta*Q*phase**m,Y,V)

def source_moment(B,X,cov,a,alpha,beta,Q,m):
    C=s.cancel(a.subs(t,z)/z)
    f=sum(expectation(Q*B**n,X,cov)*(t*C)**n/s.factorial(n) for n in range(m+1))
    return s.expand(s.diff(z**alpha*f,z,beta).subs(z,t)).coeff(t,m)*s.factorial(m)

def monomial_moment(B,X,cov,r,alpha,beta,Q,m):
    if type(r)!=int or r<1:raise ValueError('positive integer clock degree required')
    j=m-alpha+beta
    if j<0 or j%r:return s.S.Zero
    n=j//r;k=alpha+(r-1)*n
    if k<beta:return s.S.Zero
    return s.factorial(m)/s.factorial(n)*s.ff(k,beta)*expectation(Q*B**n,X,cov)

def verify_model(model):
    loc={str(v):v for v in [t,u,A,x,Z,W,eta]}
    X=tuple(loc[name] for name in model['variables'])
    B=s.sympify(model['phase'],locals=loc);a=s.sympify(model['clock'],locals=loc)
    cov=s.Matrix([[s.sympify(v) for v in row] for row in model['covariance']])
    for row in model['rows']:
        Q=s.sympify(row['marker'],locals=loc)
        for m,value in enumerate(row['values']):
            actual=direct_moment(B,X,cov,a,row['alpha'],row['beta'],Q,m)
            predicted=source_moment(B,X,cov,a,row['alpha'],row['beta'],Q,m)
            if s.expand(actual-predicted)!=0 or str(actual)!=value:raise ValueError('damaged original/source moment identity')
            if model['monomial_degree'] is not None:
                mono=monomial_moment(B,X,cov,model['monomial_degree'],row['alpha'],row['beta'],Q,m)
                if s.expand(actual-mono)!=0:raise ValueError('damaged individual monomial-clock moment')
    return True

def predicted_ideal(r,d,C,E):
    C,E=Fraction(C),Fraction(E)
    if min(C,E)<0:raise ValueError('nonnegative cutoff coordinates required')
    if C<1:return 0
    delta=r-C*d
    return None if delta<=0 else E//delta+1

def enumerated_ideal(r,d,C,E,N=8):
    C,E=Fraction(C),Fraction(E)
    cap=max(4,int(E/(1-C))+2) if C<1 else 4
    exponents=[]
    for n in range(N+1):
        q=d*n
        for alpha in range(cap+1):
            for beta in range(4):
                m=alpha+r*n-beta
                if alpha+(r-1)*n>=beta and m>=1 and m>C*(alpha+beta+q)+E:
                    exponents.append(n)
    return min(exponents) if exponents else None

def original_height(d,L,D):return min(D//d,L-1) if L else D//d
def direct_height(r,d,L,D):
    vals=[]
    for n in range(D//d+1):
        if L and n>=L:continue
        q=d*n
        for alpha in range(D-q+1):
            for beta in range(D-q-alpha+1):
                if alpha+(r-1)*n>=beta:vals.append(alpha+r*n-beta)
    return max(vals)

def main():
    models=[]
    inputs=[(Z**3,(Z,W),pair,t**2,2,[(0,0,1),(1,0,W**3),(0,1,W**3),(0,2,W**6)]),
            (x*x,(x,),s.eye(1),t**3,3,[(0,0,1),(0,1,1),(1,2,x*x)]),
            (Z*W+Z*Z,(Z,W),pair,t+t*t,None,[(0,0,1),(1,1,Z),(0,2,W)]),
            (eta*x,(x,),s.eye(1),t**2,2,[(0,0,1),(0,1,x),(2,2,x)]),
            (s.Integer(2),(x,),s.eye(1),t,1,[(0,0,1),(2,1,1)]),
            (x**3,(x,),s.eye(1),t*t+t**3,None,[(0,0,1),(0,1,x)])]
    for B,X,cov,a,r,markers in inputs:
        rows=[]
        for alpha,beta,Q in markers:
            values=[]
            for m in range(7):
                actual=direct_moment(B,X,cov,a,alpha,beta,s.sympify(Q),m)
                assert s.expand(actual-source_moment(B,X,cov,a,alpha,beta,s.sympify(Q),m))==0
                if r:assert s.expand(actual-monomial_moment(B,X,cov,r,alpha,beta,s.sympify(Q),m))==0
                values.append(str(actual))
            rows.append(dict(alpha=alpha,beta=beta,marker=str(Q),values=values))
        model=dict(phase=str(B),variables=[str(v) for v in X],covariance=[[str(v) for v in row] for row in cov.tolist()],clock=str(a),monomial_degree=r,rows=rows)
        assert verify_model(model);models.append(model)
    ideals=[]
    for r,d in [(1,2),(2,3),(5,2)]:
        for C,E in [(Fraction(1,2),0),(Fraction(1,2),2),(1,0),(1,3),(Fraction(3,2),1),(3,0)]:
            expected=predicted_ideal(r,d,C,E);got=enumerated_ideal(r,d,C,E)
            assert got==expected
            ideals.append(dict(r=r,d=d,C=str(C),E=str(E),ordinary_eta_ideal_exponent=expected,enumeration=got))
    heights=[]
    for r,d,L in [(1,2,0),(2,3,0),(5,2,0),(5,2,4),(2,3,2)]:
        for D in range(13):
            rhs=D+max(r*original_height(d,L,q)-q for q in range(D+1))
            assert direct_height(r,d,L,D)==rhs
            heights.append(dict(r=r,d=d,nilpotence_order=L or None,marker_degree=D,exact_height=rhs))
    polynomial_ideals=[]
    for a in [t+t**5,2*t*t-t**3]:
        degree=s.degree(a,t);d=2
        for C,E in [(1,0),(1,3),(Fraction(3,2),1),(3,0)]:
            exponents=[]
            # Exact source polynomials for eta*Z^2, marker W^(2n).
            for n in range(9):
                for alpha in range(3):
                    for beta in range(3):
                        poly=s.expand(s.diff(z**alpha*(t*s.cancel(a.subs(t,z)/z))**n,z,beta).subs(z,t))
                        for (m,),coefficient in s.Poly(poly,t).terms():
                            if coefficient and m>=1 and m>Fraction(C)*(2*n+alpha+beta)+E:exponents.append(n)
            actual=min(exponents) if exponents else None
            expected=predicted_ideal(int(degree),d,C,E);assert actual==expected
            polynomial_ideals.append(dict(clock=str(a),degree=int(degree),base_power=d,C=str(C),E=str(E),ordinary_eta_ideal_exponent=expected))
    # One series polynomializes; its marked derivative retains a finite pole.
    a=t*t+2*t;f=1+t;derivative=1/(2*(1+t))
    assert s.expand(f*f-1-a)==0 and s.cancel(2*f*derivative-1)==0
    assert s.denom(derivative)==2*t+2
    scalar_boundary=dict(series='sqrt(1+s)',clock=str(a),selected_composition=str(f),derivative_composition=str(derivative),
                         single_series_polynomiality_insufficient=True,abstract_scalar_control_not_claimed_new_Gaussian_example=True)
    e=s.Symbol('e')
    nil=lambda f:s.expand(f).subs(e**2,0)
    for Q in [1,x,x*x]:
        values=[direct_moment(x,(x,),s.eye(1),e*t,0,0,s.sympify(Q),m) for m in range(5)]
        assert all(s.Poly(v,e).coeff_monomial(e**k)==0 for v in values[2:] for k in [0,1])
    nonunit_boundary=dict(ring='QQ[e]/e^2',base='X',clock='e*t',lift='A+e*X',
        base_pure=False,base_safe=False,lift_pure=True,lift_safe=True,leading_unit_hypothesis_essential=True)
    # Monomial cover has distinct Laurent residue classes and does not kill this quadratic field.
    residue=[]
    for r in [1,2,3,5]:
        support={j:{j+r*k for k in range(-3,4)} for j in range(r)}
        assert all(not support[i]&support[j] for i in range(r) for j in range(i))
        assert s.gcd(1+t**r,s.diff(1+t**r,t))==1
        residue.append(dict(r=r,Laurent_residue_classes_disjoint=True,one_plus_t_to_r_squarefree=True))
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError(name)
    reject('clock constant term',lambda:clock_data(x,(x,),s.eye(1),1+t))
    reject('zero clock',lambda:clock_data(x,(x,),s.eye(1),0))
    reject('nonpolynomial clock',lambda:clock_data(x,(x,),s.eye(1),t/(1-t)))
    reject('original-variable clock',lambda:clock_data(x,(x,),s.eye(1),t*x))
    reject('singular covariance',lambda:clock_data(x,(x,),s.zeros(1),t))
    reject('nonpolynomial phase',lambda:clock_data(1/x,(x,),s.eye(1),t))
    reject('nonpositive monomial degree',lambda:monomial_moment(x,(x,),s.eye(1),0,0,0,1,1))
    damaged=json.loads(json.dumps(models[0]));damaged['rows'][0]['values'][0]='2'
    reject('damaged original moment',lambda:verify_model(damaged))
    packet=dict(original_models=models,ordinary_ideal_controls=ideals,polynomial_clock_ordinary_ideal_controls=polynomial_ideals,exact_height_controls=heights,
                single_series_boundary=scalar_boundary,nonunit_clock_boundary=nonunit_boundary,
                monomial_generic_faithfulness_controls=residue,rejections=rejected,
                universal_claims_use_written_proofs=True,actual_image_kernel_computed=False,all_ring_nilpotents_retained=True,
                arbitrary_polynomial_clock_actual_rank_preservation_claimed=False)
    out=p/'gaussian_time_change_certificate.json';out.write_text(json.dumps(packet,indent=2)+'\n')
    runtime=dict(Python=sys.version,Sympy=s.__version__,sealed_import_sha256={str(G.__file__):sha(Path(G.__file__))})
    run=p/'gaussian_time_change_runtime.json';run.write_text(json.dumps(runtime,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'POLYNOMIAL_GAUSSIAN_TIME_CHANGES_AND_EXACT_CUTOFF_REGIONS.md'),
        checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),runtime_sha256=sha(run),original_models=len(models),
        original_Wick_identities=sum(len(row['values']) for model in models for row in model['rows']),
        ordinary_ideal_controls=len(ideals),polynomial_clock_ordinary_ideal_controls=len(polynomial_ideals),exact_height_controls=len(heights),single_series_boundary_controls=1,nonunit_clock_boundary_controls=1,
        monomial_generic_faithfulness_controls=len(residue),rejections=len(rejected),historical_cusp_and_cone_restarted=False,
        independent_replay=False,formal_verification=False,external_review=False,
        assurance='Exact finite Wick/source/individual-moment and ordinary ideal/height controls;universal all-phase/all-ring transport,descent and generic faithfulness use the written proofs.')
    (p/'gaussian_time_change_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if 'sha256' not in k}))

if __name__=='__main__':main()
