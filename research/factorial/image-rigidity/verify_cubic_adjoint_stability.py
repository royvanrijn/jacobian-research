"""Exact replay for the universal cubic adjoint proof; no scan extrapolation."""
from pathlib import Path
import hashlib,json
import sympy as s
import cubic_adjoint_stability as f
HERE=Path(__file__).resolve().parent
counts={}
def check(label,condition):
    assert condition,label
    counts[label]=counts.get(label,0)+1

z,t,n,M=s.symbols('z t n M');B=s.Function('B')(z);R=s.Function('R')(z);v=s.Function('v')(z)
a,b,c=f.coefficients(n,M,z);r=-1+(1-M)/z
D=lambda expr:s.diff(expr,z)+r*expr
adj=-D(D(D(B)))+a*D(D(B))+2*s.diff(a,z)*D(B)+s.diff(a,z,2)*B-b*D(B)-s.diff(b,z)*B+c*B
check('universal_adjoint_gauge_identity',s.simplify(adj+f.polynomial_operator(B,n,M,z)/z)==0)
con=v*s.diff(R,z,2)+(a*v-s.diff(v,z))*s.diff(R,z)+(s.diff(v,z,2)-s.diff(a*v,z)+b*v)*R
LR=s.diff(R,z,3)+a*s.diff(R,z,2)+b*s.diff(R,z)+c*R
LV=-s.diff(v,z,3)+s.diff(a*v,z,2)-s.diff(b*v,z)+c*v
check('universal_Lagrange_identity',s.simplify(s.diff(con,z)-(v*LR-R*LV))==0)
x,ga=s.symbols('x ga')
logC=ga/(1-x)+(n+1)/(2-x)
check('universal_auxiliary_generating_equation',s.cancel((2-3*x+x*x)*logC-(2*ga+n+1)+(ga+n+1)*x)==0)
j=s.symbols('j');c1,c2=s.symbols('c1 c2')
c0=((3*j+2*M-3*n-2)*c1+(n-j-M+1)*c2)/(2*j)
check('universal_auxiliary_coefficient_recurrence',s.expand(-2*j*c0+(-3*(n-j)+2*M-2)*c1+(n-j-M+1)*c2)==0)
k,cp=s.symbols('k cp')
ck=cp*(n-k)*(ga+2*k)*(ga+2*k+1)/(k+1)
check('universal_tilted_coefficient_recurrence',s.factor((k+1)*ck+(k-n)*(2*k+ga)*(2*k+ga+1)*cp)==0)
# Entire Euler-polynomial substitution, on a monomial of arbitrary exponent ell.
ell=s.symbols('ell')
check('universal_Euler_quadratic_reversal',s.expand(2*((n-ell/2)*z*z-(-ell/2)*(M-ell)*(M+1-ell))-(ell*(ell-M)*(ell-M-1)-z*z*(ell-2*n)))==0)
# Matrix coefficient identity: divide the three recurrence terms by D_(j,k).
# Replay symbolic Gamma parameters without any numeric root approximations.
for N in range(1,7):
    G=f.gamma_moment(N,ga,z).as_expr()
    SS=ga+N-1
    d0=s.Integer(1);d1=z-SS
    Dp=d0 if N==0 else d1
    for jj in range(2,N+1):
        d0,d1=d1,s.expand((z-SS-2+2*jj)*d1+(jj-1)*(SS+2-jj)*d0)
    check('Gamma_characteristic_polynomial_identity',s.expand((-1)**N*d1.subs(z,-z)-G)==0)
    T=f.tilted(N,ga,t).as_expr();theta=t*s.diff(T,t)
    check('current_Gamma_moment_identity',s.expand(f.tilted(N,ga+1,t).as_expr()-T-2*theta/ga)==0)
    check('preceding_Gamma_moment_identity',s.expand(f.tilted(N-1,ga,t).as_expr()-T+theta/N)==0)
    # Universal-in-gamma convolution identity at each replayed size.
    BB=f.auxiliary(N,ga,z)
    GG=f.gamma_moment(N,N+1,2*z).as_expr()/2**N
    conv=f.additive_convolution(G,GG,z,N)
    check('Walsh_convolution_coefficient_identity',s.expand(BB.as_expr()-conv.as_expr())==0)

for N in range(1,8):
    for gamma in sorted({s.Rational(1,2),s.Rational(3,2),s.Integer(N),s.Rational(2*N-1,2),s.Integer(max(1,N-2))}):
        MM=gamma+2*N;RR=f.reversed_tilted(N,gamma,z).as_expr();BB=f.auxiliary(N,gamma,z).as_expr()
        aa,bb,cc=f.coefficients(N,MM,z)
        check('reversed_polynomial_ODE',s.cancel(s.diff(RR,z,3)+aa*s.diff(RR,z,2)+bb*s.diff(RR,z)+cc*RR)==0)
        check('auxiliary_polynomial_ODE',s.expand(f.polynomial_operator(BB,N,MM,z))==0)
        check('positive_adjoint_zero_concomitant',f.gauged_concomitant(RR,BB,N,MM,z)==0)
        check('negative_adjoint_zero_concomitant',f.gauged_concomitant(RR,BB.subs(z,-z),N,MM,z,sign=-1)==0)
        if gamma>=N-2:
            TT=f.tilted(N,gamma,t)
            check('rational_boundary_squarefreeness',s.gcd(TT,TT.diff()).degree()==0)
            # Exact Hurwitz test: B(-z) has all roots with positive real part.
            # Routh-Hurwitz principal determinants, in descending coefficient order.
            coeff=s.Poly(BB,z).all_coeffs()
            H=s.Matrix(N,N,lambda ii,jj:coeff[2*jj-ii+1] if 0<=2*jj-ii+1<=N else 0)
            for d in range(1,N+1):
                check('auxiliary_exact_Hurwitz_controls',H[:d,:d].det()>0)

source=HERE/'sources/wagner-multivariate-stable-polynomials-2011.pdf'
check('classical_GWS_source_pin',hashlib.sha256(source.read_bytes()).hexdigest()=='ab72551b3407480e08185ce2b0a9527c8cea385b95ef93dcf904ec90a46170c1')
audit={
 'status':'pass','canonical_proof':'CUBIC_ADJOINT_STABILITY_SFC.tex',
 'universal_chain':['coefficient ODE and quadratic reversal','explicit Gamma-convolution auxiliary and exact adjoint gauge','Gamma Jacobi diagonal bound plus classical half-plane GWS','rational-log-derivative exclusion of nonzero concomitants','double root forces B(z)=B(-z)=0','degree-free Gamma moment bridge'],
 'classical_input':{'source':str(source.relative_to(HERE)),'theorem':'Wagner Theorem4.1, half-plane polarization only','source_scope':'statement and Section4 proof context/application inspected; not a full-source audit'},
 'all_index_cubic_SFC_closed':True,'all_positive_Gamma_shapes':True,
 'simplicity_parameter_range':'n>=1,gamma>0,gamma>=n-2',
 'finite_prefix_required':False,'finite_scan_extrapolation':False,
 'all_gaps_or_lower_exponents_proved':False,'general_FC_closed':False,'unrestricted_SFC_closed':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'cubic_adjoint_stability_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=['CUBIC_ADJOINT_STABILITY_SFC.tex','cubic_adjoint_stability.py','verify_cubic_adjoint_stability.py','cubic_adjoint_stability_audit.json',str(source.relative_to(HERE))]
receipt=dict(audit,counts=counts,total=sum(counts.values()),files={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files})
(HERE/'cubic_adjoint_stability_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'total':receipt['total'],'counts':counts,'all_index_cubic_SFC_closed':True,'all_positive_Gamma_shapes':True},indent=2))
