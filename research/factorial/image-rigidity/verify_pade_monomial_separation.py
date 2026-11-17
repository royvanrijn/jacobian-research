"""Verify all-n algebraic strip and shifted positivity certificates.

The universal proof is PADE_ADJACENT_MONOMIALS.tex. The finite polynomial
replays below check bookkeeping; they do not extrapolate an index scan.
"""
from pathlib import Path
from math import factorial
import hashlib,json
import sympy as s
import pade_monomial_separation as p
HERE=Path(__file__).resolve().parent
z,n,u=p.Z,p.N,p.U
counts={}
def check(category,condition):
    assert condition,category
    counts[category]=counts.get(category,0)+1

# The generic matrix recurrence is checked in symbolic integer parameters.
S,j,k=s.symbols('S j k',integer=True,positive=True)
def A(j,k):return (-1)**(j-k)*s.factorial(j)/s.factorial(k)*s.factorial(S-k)/(s.factorial(j-k)*s.factorial(S-j))
check('symbolic_matrix_coefficient_recurrence',s.combsimp(A(j,k)-A(j-1,k-1)+(S+2-2*j)*A(j-1,k)-(j-1)*(S+2-j)*A(j-2,k))==0)
check('symbolic_matrix_constant_recurrence',s.combsimp(A(j,0)+(S+2-2*j)*A(j-1,0)-(j-1)*(S+2-j)*A(j-2,0))==0)
check('symbolic_matrix_subleading_recurrence',s.expand(-j*(S-j+1)+(j-1)*(S-j+2)+(S+2-2*j))==0)
a,nn=s.symbols('a nn',integer=True,positive=True)
check('symbolic_strip_difference',s.expand((a-1)*nn+2-(a+1)*(nn-1)-(a+3-2*nn))==0)

for M in range(0,7):
    for N in range(0,6):
        P,Q=p.pade_pair(M,N)
        for degree in range(M+N+1):
            coefficient=sum(Q.nth(t)/s.factorial(degree-t) for t in range(min(N,degree)+1))
            check('Pade_coefficient_replays',coefficient==P.nth(degree))
        if N:
            d0=s.Poly(1,z);d1=s.Poly(z-(M+N),z)
            for row in range(2,N+1):
                d0,d1=d1,s.Poly((z-M-N-2+2*row)*d1.as_expr()+(row-1)*(M+N+2-row)*d0.as_expr(),z)
            expected=(-1)**N*s.Rational(factorial(M),factorial(M+N))*d1.as_expr()
            check('matrix_polynomial_replays',s.expand(Q.as_expr()-expected)==0)

closed={
 1:s.Integer(2),
 2:z+9*n-3,
 3:z**2+2*(4*n-1)*z+8*(2*n-1)*(4*n-1),
 4:z**3+(11*n-3)*z**2+3*(5*n-2)*(5*n-1)*z+5*(5*n-3)*(5*n-2)*(5*n-1),
 5:z**4+2*(7*n-2)*z**3+6*(3*n-1)*(7*n-2)*z**2+24*(2*n-1)*(3*n-1)*(6*n-1)*z+72*(2*n-1)*(3*n-2)*(3*n-1)*(6*n-1),
 6:z**5+(17*n-5)*z**4+10*(19*n**2-12*n+2)*z**3+2*(7*n-3)*(7*n-2)*(17*n-5)*z**2+5*(7*n-4)*(7*n-3)*(7*n-2)*(7*n-1)*z+7*(7*n-5)*(7*n-4)*(7*n-3)*(7*n-2)*(7*n-1)}
for aa,H in closed.items():
    check('symbolic_auxiliary_templates',s.expand(p.auxiliary_symbolic(aa).as_expr()-H)==0)
    for index in range(2,7):
        P,Q=p.pade_pair(aa*index,index)
        P0,Q0=p.pade_pair(aa*(index-1),index-1)
        ell=(aa+1)*(index-1)+1
        C=s.Rational(factorial(index),factorial((aa+1)*index))*s.Rational((-1)**(index-1)*factorial(aa*(index-1)),factorial((aa+1)*(index-1)))
        check('cross_determinant_replays',s.expand((P*Q0-P0*Q).as_expr()-C*z**ell*H.subs(n,index))==0)
        B=sum(s.binomial(index,t)*s.Rational(factorial(aa*index+t),factorial(aa*index))*(-1)**t*z**(index-t) for t in range(index+1))
        factor=s.Rational((-1)**index*factorial(aa*index),factorial((aa+1)*index))
        check('moment_Pade_reciprocal_replays',s.expand(Q.as_expr()-factor*B)==0)

certificates=[]
for aa in range(2,7):
    cert=p.shifted_all_index_certificate(aa)
    check('all_n_shift_certificates',cert['status']=='all_indices_certified' and cert['all_n_ge2_certified'])
    H=p.auxiliary_symbolic(aa);shifted=s.Poly(s.expand(H.as_expr().subs(z,z+(aa-1)*n+2)),z)
    params=[s.sympify(t,locals={'n':n}) for t in cert['b_parameters']]
    check('symbolic_Schwarz_reconstruction',s.simplify(p.reconstruct(params).as_expr()-shifted.as_expr())==0)
    for b,ratio in zip(params,cert['positive_ratio_certificates']):
        num=[s.Rational(t) for t in ratio['numerator']];den=[s.Rational(t) for t in ratio['denominator']]
        check('positive_polynomial_ratios',all(t>=0 for t in num+den) and num[0]>0 and den[0]>0)
        exact=(sum(c*u**i for i,c in enumerate(num))/sum(c*u**i for i,c in enumerate(den))).subs(u,n-2)
        check('positive_ratio_identity',s.cancel(exact-b)==0)
    certificates.append(cert)

# Independent determinant and Hermitian-part checks of the dissipative
# matrix construction, in abstract positive parameters, degree1..5.
for degree in range(1,6):
    b0=s.Symbol('b0',positive=True);r=s.symbols('r1:'+str(degree),positive=True)
    T=s.zeros(degree);T[0,0]=-b0
    for i in range(degree-1):T[i,i+1]=-r[i];T[i+1,i]=r[i]
    check('dissipative_Hermitian_identity',T+T.T==s.diag(-2*b0,*([0]*(degree-1))))
    chars=s.Poly(T.charpoly(z).as_expr(),z)
    check('dissipative_characteristic_identity',s.expand(chars.as_expr()-p.reconstruct([b0]+[t**2 for t in r]).as_expr())==0)

check('failed_real_rooted_shortcut',s.discriminant(closed[3].subs(n,2),z)<0)
q=s.Poly(closed[5].subs(n,2),z).all_coeffs();hurwitz3=q[1]*q[2]*q[3]-q[3]**2-q[1]**2*q[4]
check('failed_unshifted_Hurwitz_shortcut',hurwitz3<0)
# The generic local test remains explicitly sufficient; zero or negative
# indexed answers do not claim mathematical failure.
for aa,index,expected in [(1,2,True),(3,3,True),(6,4,True),(7,6,False),(8,6,False)]:
    row=p.indexed_strip_certificate(aa,index)
    check('indexed_strip_boundaries',row['accepted']==expected and not row['counterexample_claimed'])

for aa in range(1,7):
    for index in range(2,5):
        check('literal_integer_auxiliary_identity',s.expand(p.auxiliary_at_integer(aa,index).as_expr()-closed[aa].subs(n,index))==0)
finite_diagnostics=[]
for aa,index in [(10,20),(20,50),(40,100)]:
    row=p.finite_shift_certificate(aa,index)
    check('finite_shift_diagnostics',row['accepted'] and not row['all_indices_claimed'] and not row['all_a_claimed'])
    decoded=[s.Rational(int(v['numerator'],16),int(v['denominator'],16)) for v in row['b_parameters_hex']]
    check('finite_witness_encoding',all(v>0 for v in decoded) and s.expand(p.reconstruct(decoded).as_expr()-p.auxiliary_at_integer(aa,index).as_expr().subs(z,z+(aa-1)*index+2))==0)
    finite_diagnostics.append({'a':aa,'terminal_index':index,'accepted':row['accepted'],
        'outside_registered_a_range':aa>6,'outside_indexed_strip_family':aa<2*index-3,
        'b_signs':['positive']*len(row['b_parameters_hex']),
        'unbounded_extrapolation_claimed':False})
(HERE/'pade_separation_diagnostics.json').write_text(json.dumps({'status':'three_finite_terminal_certificates',
    'controls':finite_diagnostics,'all_a_or_all_n_claimed':False,
    'replay':'finite_shift_certificate(a,n) in pade_monomial_separation.py'},indent=2)+'\n')

t,formal=s.symbols('t formal')
G=s.exp(z*formal)*(1+formal)**(-t-1)
check('generating_recurrence_forward',s.simplify((t+1)*G/(1+formal)-(t+1+z)*G-formal*s.diff(G,formal)+z*(1+formal)*G)==0)
check('generating_recurrence_lowering',s.simplify((1+formal)*s.diff(G,formal)-(z*(1+formal)-t-1)*G)==0)
for aa in range(2,7):
    A,B=p.pencil_symbolic(aa)
    check('pencil_determinant',s.factor(B.det())==s.Rational(1,aa+1))
    # A tridiagonal continuant supplies a separate reconstruction.
    mat=A+z*B;f0=s.Integer(1);f1=mat[0,0]
    for j in range(1,aa-1):f0,f1=f1,s.expand(mat[j,j]*f1-mat[j-1,j]*mat[j,j-1]*f0)
    check('pencil_auxiliary_identity',s.expand((aa+1)*f1-closed[aa])==0)

payload={'status':'all_n_certificates_verified','parameter':'u=n-2>=0',
         'a_values':[2,3,4,5,6],'certificates':certificates,
         'universal_evidence':'PADE_ADJACENT_MONOMIALS.tex: strip, cross determinant and dissipative matrices',
         'all_a_claimed':False,'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'pade_shift_positive_certificates.json').write_text(json.dumps(payload,indent=2)+'\n')
audit={'proved_all_starting_indices':'ordinary X^a(A+B*X), a=1,2,3,4,5,6; a0 previously covered',
       'indexed_unbounded_family':'all n>=2,a>=1 with a>=2*n-3',
       'denominator_strip':'M-N+2<=Re(z)<=M+N; strict for N>=2',
       'explicit_pencil_target':'Every eigenvalue of -B^-1*A has real part<(a-1)*n+2; necessary pencil proved, uniform spectral inequality open',
       'generic_criterion':'positive Schwarz reconstruction of H_(a,n)(z+(a-1)*n+2) excludes every adjacent common root',
       'shift_certificates_cover_all_n_ge2':True,'finite_index_extrapolation_used':False,
       'failed_shortcuts':{'H3_n2_discriminant':str(s.discriminant(closed[3].subs(n,2),z)),
                           'H5_n2_Hurwitz3':str(hurwitz3)},
       'provenance':'Elementary Padé coefficient, finite-difference and matrix derivations. Rocks2015 author draft scope compared; no full final-publication or priority audit.',
       'all_a_proved':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
       'other_Gamma_shapes_claimed':False,'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'pade_monomial_separation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=('PADE_ADJACENT_MONOMIALS.tex','pade_monomial_separation.py','verify_pade_monomial_separation.py',
       'pade_shift_positive_certificates.json','pade_monomial_separation_audit.json','pade_separation_diagnostics.json',
       'sources/rocks-2015-author-preprint.pdf','sources/rocks-2015-author-preprint.txt')
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'all_n_a1_through6_proved':True,'unbounded_indexed_strip_family_proved':True,
         'positive_symbolic_certificates_all_n':True,'all_a_proved':False,
         'unrestricted_SFC_closed':False,'general_FC_closed':False,
         'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'pade_monomial_separation_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'all_n_a_values':[1,2,3,4,5,6],
                 'indexed_family':'a>=2n-3','all_a_proved':False,'unrestricted_SFC_closed':False}))
