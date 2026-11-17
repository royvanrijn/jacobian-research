"""Check the algebra supporting a universal shared-parabola separation.

Finite replays verify bookkeeping. The universal all-a, all-n and
all-positive-shape assertions are proved in the companion TeX note.
"""
from pathlib import Path
import hashlib,json
import sympy as s
import pade_parabola_sfc as p
HERE=Path(__file__).resolve().parent
z=p.Z
counts={}
def check(category,condition):
    assert condition,category
    counts[category]=counts.get(category,0)+1

# Exact universal geometric identities. Strict inequalities follow from
# the stated domains, not a numerical sample of eigenvalues.
N,S,tau,alpha,u,w,x,y=s.symbols('N S tau alpha u w x y',real=True)
check('universal_parabola_gap',s.expand(N*(S-2*tau+N)-tau*(S-tau)-(N-tau)*S-(N-tau)**2)==0)
check('universal_invariant_completion',s.expand(u**2+w**2+u*x+w*y-((u-alpha)*(u+alpha+x)+(w+y/2)**2+alpha**2+alpha*x-y**2/4))==0)
M=s.symbols('M',real=True)
check('parabola_x_boundary',s.expand((M+N)-2*(N-1)+N-(M+2))==0)
check('Cauchy_upper_sum_boundary',s.expand(S+1-N-(M+1)).subs(S,M+N)==0)

# Characteristic coefficients checked at a symbolic recurrence index.
j,k=s.symbols('j k',integer=True,positive=True)
def coefficient(j,k):
    return (-1)**(j-k)*s.factorial(j)/s.factorial(k)*s.factorial(S-k)/(s.factorial(j-k)*s.factorial(S-j))
check('matrix_coefficient_recurrence',s.combsimp(coefficient(j,k)-coefficient(j-1,k-1)+(S+2-2*j)*coefficient(j-1,k)-(j-1)*(S+2-j)*coefficient(j-2,k))==0)
check('matrix_constant_recurrence',s.combsimp(coefficient(j,0)+(S+2-2*j)*coefficient(j-1,0)-(j-1)*(S+2-j)*coefficient(j-2,0))==0)
check('matrix_top_coefficient',s.simplify(coefficient(j,j)-1)==0)

formal,t=s.symbols('formal t')
G=s.exp(z*formal)*(1+formal)**(-t-1)
check('forward_parameter_recurrence',s.simplify((t+1)*G/(1+formal)-(t+1+z)*G-formal*s.diff(G,formal)+z*(1+formal)*G)==0)
check('degree_lowering_recurrence',s.simplify((1+formal)*s.diff(G,formal)-(z*(1+formal)-t-1)*G)==0)
check('fixed_shape_backward_recurrence',s.simplify(s.diff(G,formal)-(z-t-1)*G+formal*s.diff(G,formal)-z*formal*G)==0)

l=s.symbols('l',real=True)
last_b=(M+N)*(M+N-1)/N
last_c=(M+N-1)*(M+N-2)/(M-1)
check('last_parabola_parameter',s.cancel(last_b*(1-(M+N-2)/last_c)-(M+N))==0)
interior_b=l+N+j;interior_c=interior_b*(interior_b-1)/(l+j)
check('interior_parabola_parameter',s.cancel(interior_b*(1-(interior_b-1)/interior_c)-N)==0)
check('boundary_elimination',s.cancel(M+N-1+z-M*z/(M+N)-(M+N-1+N*z/(M+N)))==0)
# The single-row case has M=l+2>1, so this is positive.
check('single_row_first_margin',s.expand(((M+N)*(M+N-1)-N*N)-(M*(M-1)+N*(2*M-1)))==0)

# Representative exact rational parameters include the boundary M close
# to -1 for the inside lemma, and noninteger Gamma shapes for transfer.
for MM in [s.Rational(-9,10),s.Rational(-1,2),0,s.Rational(7,3),5]:
    for NN in range(1,8):
        SS=MM+NN
        T=s.diag(*[SS+2-2*j for j in range(1,NN+1)])
        for j in range(1,NN):T[j-1,j]=T[j,j-1]=s.I*s.sqrt(j*(SS+1-j))
        check('matrix_characteristic_replays',s.expand(T.charpoly(z).as_expr()-p.R(MM,NN).as_expr())==0)
        check('matrix_Hermitian_replays',T+T.conjugate().T==2*s.diag(*[SS+2-2*j for j in range(1,NN+1)]))

for NN in range(2,8):
    for ll in map(s.sympify,[s.Rational(-1,2),0,1,3,s.Rational(7,3)]):
        for aa in [1,2,4,7]:
            MM=ll+aa
            check('common_root_diagnostics',s.gcd(p.R(MM,NN),p.R(ll,NN-1)).degree()==0)
            if aa>=2:
                A,B=p.pencil(ll,aa,NN);rows=p.normalized_continuants(ll,aa,NN)
                constant=s.prod(ll+NN+j for j in range(1,aa))
                check('independent_pencil_continuants',s.expand((A+z*B).det(method='domain-ge')-constant*rows[-1].as_expr())==0)

for beta in map(s.Rational,[s.Rational(1,7),1,s.Rational(7,2)]):
    for aa in [0,1,2,8,100]:
        for NN in [1,2,5]:
            cert=p.adjacent_monomial_certificate(aa,NN,beta)
            check('structural_certificate_domains',cert['accepted'] and not cert['unrestricted_SFC_closed'])
            Bmoment=sum(s.binomial(NN,k)*s.rf(beta+aa*NN,k)*(-1)**k*z**(NN-k) for k in range(NN+1))
            check('Gamma_moment_reciprocal_identity',s.expand(p.R(beta+aa*NN-1,NN).as_expr()-Bmoment)==0)
            if NN>=2:
                check('Gamma_adjacent_exact_replays',s.gcd(p.R(beta+aa*NN-1,NN),p.R(beta+aa*(NN-1)-1,NN-1)).degree()==0)

# Exact next-boundary identities, with t still symbolic. These do not
# upgrade the proof scope to nonadjacent supports.
for gap in [1,2,3]:
    for NN in range(1,6):
        V=p.gap_R(t,gap,NN).as_expr();W=p.gap_R(t,gap,NN-1).as_expr()
        check('gap_Pascal_identity',s.expand(V-z*W+s.rf(t+1,gap)*W.subs(t,t+gap))==0)
        check('gap_parameter_raising',s.expand((t+1)*V.subs(t,t+1)-(t+1+gap*NN)*V+gap*z*NN*W)==0)
        check('gap_sparse_parameter_recurrence',s.expand(z*(W-W.subs(t,t+1))-s.rf(t+2,gap-1)*((t+gap*NN+1)*W.subs(t,t+gap)-(t+gap+1)*W.subs(t,t+gap+1)))==0)
        if gap==1:check('gap_one_reconciliation',s.expand(V-p.R(t,NN).as_expr())==0)

proof_audit={
    'proof':'UNIFORM_PADE_PARABOLA_SFC.tex',
    'universal_scope':'ALL integer a>=0, terminal n>=1, Gamma beta>0: two consecutive moments of nonzero X^a(A+B*X) never both vanish',
    'generalized_coprimality':'R_N(l+a,z) and R_(N-1)(l,z) coprime for all N>=2, real l>-1 and integer a>=1',
    'inside_parabola':'For every real M>-1, every zero of R_N(M,z) has y^2<4N(x+N)',
    'outside_parabola':'The necessary degree-(a-1) continuant is nonzero in the same parabola, by a proved invariant continued-fraction inequality',
    'nonadjacent_gap_identities_only':True,'finite_scan_extrapolation':False,'vertical_bound_required':False,
    'independence':'No contour, MR density, finite-prefix squarefreeness, parity, irreducibility or prime-gap input',
    'source_comparison':'Elementary parabola principle compared with Saff-Varga1976 Theorem1.1; fully reproved here. Only that displayed theorem inspected, no whole-paper or priority audit.',
    'arbitrary_monomial_gaps_closed':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'pade_parabola_sfc_audit.json').write_text(json.dumps(proof_audit,indent=2)+'\n')
files=['UNIFORM_PADE_PARABOLA_SFC.tex','pade_parabola_sfc.py','verify_pade_parabola_sfc.py','pade_parabola_sfc_audit.json','sources/saff-varga-1976-sharpness.pdf']
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'all_a_all_n_adjacent_monomial_SFC':True,'all_positive_Gamma_shapes':True,
    'generalized_Pade_coprimality':True,'nonadjacent_gap_identities_only':True,'finite_scan_extrapolation':False,'vertical_bound_required':False,
    'arbitrary_monomial_gaps_closed':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'pade_parabola_sfc_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'all_a':True,'all_n':True,'all_positive_Gamma_shapes':True,'unrestricted_SFC_closed':False}))
