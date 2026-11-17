"""Exact replay of the uniform polar theorem, backward certificates, and limits."""
from pathlib import Path
from math import factorial,comb
import hashlib,json,copy
import sympy as s
import gap_three_polar_continuant as f
import gap_three_wronskian as wr
import higher_gap_adjoint_rank as old
from prime_power_windows import newton_residuals
HERE=Path(__file__).resolve().parent
counts={}
def check(label,condition):
    assert condition,label
    counts[label]=counts.get(label,0)+1
z,y,x,u,g=s.symbols('z y x u g',real=True)
for n in range(1,7):
    B=old.auxiliary(n,g,3,z);P=f.parameter_product(n,g,x)
    # P(k) has backwards sample values. Its finite-difference transform is exactly B.
    diffs=[sum((-1)**k*s.binomial(j,k)*P.eval(k) for k in range(j+1)) for j in range(n+1)]
    check('uniform_parameter_finite_difference_identity',s.expand(B.as_expr()-sum(diffs[j]*z**j/s.factorial(j) for j in range(n+1)))==0)
    J,K=f.companions(B,n,g,z,y)
    pol=(1+y)**n
    for j in range(n):pol=s.expand((g+2+3*j)*pol-y*s.diff(pol,y))
    check('polar_product_equals_binomial_companion',s.expand(pol-J.as_expr())==0)
    check('Möbius_companion_identity',s.cancel(J.as_expr()-(1+y)**n*K.as_expr().subs(z,-y/(1+y)))==0)
    check('nondegenerate_Möbius_endpoints',s.expand(J.eval(-1)-s.factorial(n))==0 and s.expand(K.eval(-1)-P.eval(n))==0)
    for beta in [s.Rational(1,2),s.Integer(1),s.Integer(7)]:
        gamma=n+beta-1
        JJ=s.Poly(J.as_expr().subs(g,gamma),y);KK=s.Poly(K.as_expr().subs(g,gamma),z)
        check('finite_companion_roots_strictly_below_minus1',int(JJ.count_roots(-s.oo,-1))==n and int(KK.count_roots(-s.oo,-1))==n)
        check('finite_companion_simplicity',s.gcd(JJ,JJ.diff()).degree()==0 and s.gcd(KK,KK.diff()).degree()==0)

# Exact inverse-companion counterprofile, outside the arithmetic-progression Gamma family.
n=8;qs=[s.Rational(90+j,108) for j in range(n)]
check('counter_companion_distinct_roots_below_minus1',all(0<q<1 for q in qs) and len(set(qs))==n)
Kstar=s.Poly(s.prod(1+q*z for q in qs),z)
Bstar=s.Poly(sum(Kstar.nth(j)/s.ff(n,j)*z**j for j in range(n+1)),z).monic()
W=old.gap_three_W(Bstar.as_expr(),z);parts=wr.cyclic_sections(W,z,u)
realcounts=[int(p.count_roots(-s.oo,s.oo)) for p in parts]
check('inverse_companion_does_not_make_any_section_real_rooted',[p.degree() for p in parts]==[5,5,4] and realcounts==[3,3,0])
check('generic_counter_sections_are_squarefree',all(s.gcd(p,p.diff()).degree()==0 for p in parts))

symbolic_profiles=[]
for n in [2,3]:
    parts=wr.cyclic_sections(wr.wronskian(n,g,z),z,u)
    cert=f.backward_certificate(parts,2*n,u)
    check('symbolic_certificate_retains_nonzero_conditions',cert['status']=='formal_nonzero_conditions')
    for j,pivot in cert['weights'].items():
        num,den=s.fraction(s.factor(pivot))
        check('finite_index_all_positive_shape_pivot_certificate',all(c>0 for c in s.Poly(s.expand(num),g).all_coeffs()) and all(c>0 for c in s.Poly(s.expand(den),g).all_coeffs()))
    symbolic_profiles.append({'n':n,'scope':'ALL gamma>0 at this FIXED index only','weights':{str(j):str(w) for j,w in sorted(cert['weights'].items())}})

numeric_profiles=[]
for n in [8,24]:
    parts=wr.cyclic_sections(wr.centered(n,1,z),z,u)
    cert=f.backward_certificate(parts,2*n,u)
    check('ordinary_backward_certificate_replayed',f.replay_certificate(cert,2*n,u))
    check('ordinary_backward_certificate_positive_weights',all(w>0 for w in cert['weights'].values()))
    numeric_profiles.append({'n':n,'gamma':n,'weights':{str(j):str(w) for j,w in sorted(cert['weights'].items())},'q':{str(j):[str(c) for c in p.all_coeffs()] for j,p in sorted(cert['q'].items())}})
    tampered=copy.deepcopy(cert);tampered['weights'][2*n]+=1
    check('tampered_backward_certificate_rejected',not f.replay_certificate(tampered,2*n,u))

zero=[s.Poly(1+u+u*u,u),s.Poly(1+u,u),s.Poly(1+2*u,u)]
result=f.backward_certificate(zero,6,u)
check('zero_pivot_is_inconclusive_not_failure',result['status']=='inconclusive_zero_pivot' and result['at']==6 and s.gcd(s.gcd(zero[0],zero[1]),zero[2]).degree()==0)

large_controls=[]
for n in [24,48]:
    pp=wr.cyclic_sections(wr.centered(n,1,z),z,u)
    negative=[int(p.count_roots(-s.oo,0)) for p in pp]
    check('larger_index_diagnostic_only',negative==[p.degree() for p in pp])
    large_controls.append({'n':n,'gamma':n,'degrees':[p.degree() for p in pp],'negative_root_counts':negative,'uniform_claim':False})

native=[]
for n,expected in [(22,[([19,27],[22,33],[1,0,0,1],1)]),
                    (63,[([54,80],[60,89],[1,0,0,2],2),([60,89],[63,95],[1,0,0,1],1)])]:
    coeff=[comb(n,k)*factorial(n+3*k-1)//factorial(n-1) for k in range(n+1)]
    bad=[row for row in newton_residuals(coeff,3) if not row['residual_squarefree']]
    observed=[(row['left'],row['right'],row['residual_coefficients'],row['e']) for row in bad]
    check('native_prime_residual_boundary',observed==expected)
    native.append({'n':n,'bad_edges':bad,'SFC_counterexample':False})
check('native_prime_and_parity_still_permit_quadratic',native[1]['bad_edges'][0]['e']==2 and native[1]['bad_edges'][0]['residual_coefficients']==[1,0,0,2])

source='sources/iskander-jain-2020-delta-laguerre-polya.pdf'
check('finite_difference_source_pin',hashlib.sha256((HERE/source).read_bytes()).hexdigest()=='4d9a13762ee08160e0441c44bdbaec6eacbc6f8f714af60eb812017f4643d709')
audit={'status':'pass','canonical_proof':'GAP_THREE_POLAR_AND_CONTINUANT.tex',
 'uniform_real_rooted_polar_companion_proved':True,'companion_roots_simple_below_minus1_proved':True,
 'generic_inverse_transfer_to_real_rooted_individual_section_refuted':True,
 'nonzero_pivot_backward_coprimality_criterion_proved':True,'positive_weights_required_for_coprimality':False,
 'zero_pivot_fail_closed':True,'native_prime_squarefree_residual_shortcut_refuted':True,
 'uniform_pivot_nonvanishing_proved':False,'uniform_section_existence_proved':False,
 'all_index_gap_three_SFC_closed':False,'general_FC_closed':False,'unrestricted_SFC_closed':False,
 'finite_extrapolation':False,'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False,
 'source_scope':'Iskander-Jain arXiv v1 Theorems1.1/1.2 and finite/infinite zero-free hypotheses inspected. The needed polar lemma is reproved independently. No full-paper/dependency or priority audit.',
 'boundary':'The arithmetic-progression Gamma structure still needs a valid transfer or uniform nonzero pivots. Generic companion counterexample is not a Gamma moment pair or SFC counterexample.'}
certificates={'symbolic_fixed_index_profiles':symbolic_profiles,'numeric_fixed_index_profiles':numeric_profiles,'uniform_index_claim':False}
diagnostics={'generic_inverse_companion':{'n':8,'Kstar_factors':[str(q) for q in qs],'section_degrees':[5,5,4],'section_real_root_counts':realcounts,'is_actual_Gamma_family':False,'SFC_counterexample':False},'larger_index_controls':large_controls,'native_prime_boundaries':native,'source_imported_for_proof':False}
for name,data in [('gap_three_polar_continuant_audit.json',audit),('gap_three_polar_continuant_certificates.json',certificates),('gap_three_polar_continuant_diagnostics.json',diagnostics)]:
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
files=['GAP_THREE_POLAR_AND_CONTINUANT.tex','gap_three_polar_continuant.py','verify_gap_three_polar_continuant.py','gap_three_polar_continuant_audit.json','gap_three_polar_continuant_certificates.json','gap_three_polar_continuant_diagnostics.json','gap_three_wronskian.py','higher_gap_adjoint_rank.py','prime_power_windows.py',source]
receipt=dict(audit,total=sum(counts.values()),counts=counts,files={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files})
(HERE/'gap_three_polar_continuant_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'total':receipt['total'],'counts':counts,'all_index_gap_three_SFC_closed':False},indent=2))
