"""Exact replay of the Wronskian recurrence, real-section criterion inputs and fifth order."""
from pathlib import Path
import hashlib,json
import sympy as s
import gap_three_wronskian as f
import higher_gap_adjoint_rank as old
HERE=Path(__file__).resolve().parent
counts={}
def check(key,condition):
    assert condition,key
    counts[key]=counts.get(key,0)+1
z,u,g,beta=s.symbols('z u g beta',real=True)
# Universal arbitrary-function transformation of a pair Wronskian.
a=s.symbols('a');v=s.Function('v')(z);w=s.Function('w')(z)
Wr=lambda v,w:v*s.diff(w,z)-s.diff(v,z)*w
lhs=Wr(a*v-z*s.diff(v,z),a*w-z*s.diff(w,z))
rhs=a*(a-1)*Wr(v,w)-a*z*s.diff(Wr(v,w),z)+z*z*Wr(s.diff(v,z),s.diff(w,z))
check('universal_pair_Wronskian_transformation',s.expand(lhs-rhs)==0)
for n in range(0,7):
    B=old.auxiliary(n,g,3,z).as_expr();Bnext=old.auxiliary(n+1,g,3,z).as_expr()
    check('all_parameter_shape_shift_replay',s.expand(old.auxiliary(n,g-1,3,z).as_expr())==s.expand(B-s.diff(B,z)))
    U=f.wronskian(n,g,z)
    V=f.wronskian(n,g-1,z)
    check('symbolic_gamma_fixed_shape_recurrence',s.expand(f.next_fixed_shape(U.as_expr(),V.as_expr(),n,g,z).as_expr()-f.wronskian(n+1,g,z).as_expr())==0)
    check('exact_constant_coefficient_product',s.factor(U.nth(0)-s.prod((g+3*j+1)*(g+3*j+2) for j in range(n)))==0)
    if n:
        check('top_two_coefficient_identity',s.expand(U.nth(2*n-1)-n*(g+n+1))==0)
    if n<=5:
        direct=old.gap_three_W(B,z)
        check('rational_recurrence_matches_radical_definition',s.expand(U.as_expr()-direct.as_expr())==0)
        alpha=(1+s.sqrt(-3))/2;alpha_bar=s.conjugate(alpha)
        cd=s.expand(Bnext.subs(z,alpha*z)*B.subs(z,alpha_bar*z)-B.subs(z,alpha*z)*Bnext.subs(z,alpha_bar*z))
        check('consecutive_auxiliary_divided_difference',s.simplify(s.expand(cd-(alpha-alpha_bar)*z*U.as_expr()))==0)
    coeff=f.centered_coefficients(n,beta)
    for c in coeff:
        check('finite_bivariate_coefficient_positivity_controls',all(v>=0 for v in s.Poly(c,beta).all_coeffs()) and c.subs(beta,0)>0)

# The fifth-order polynomial and its strictly positive discriminant are symbolic in gamma.
A,B,C=f.fifth_section_coefficients(g);p=f.cyclic_sections(f.wronskian(5,g,z),z,u)[2]
check('fifth_section_symbolic_identity',s.expand(p.as_expr()-(A*u*u+B*u+C))==0)
disc=s.expand(B*B-4*A*C)
check('fifth_discriminant_exact_coefficients',tuple(s.Poly(disc,g).all_coeffs())==f.FIFTH_DISCRIMINANT_DESCENDING)
check('fifth_discriminant_all_positive_coefficients',all(c>0 for c in s.Poly(disc,g).all_coeffs()))
for expr in [A,B,C]:check('fifth_section_all_positive_coefficients',all(c>0 for c in s.Poly(s.expand(expr),g).all_coeffs()))
# Algebraic shape outside the allowed range has a genuine repeated root and common NONREAL section root.
gamma=(-15+s.sqrt(-15))/6
W=f.wronskian(2,gamma,z);pp=f.cyclic_sections(W,z,u)
GG=s.gcd(s.gcd(s.Poly(pp[0].as_expr(),u,extension=True),s.Poly(pp[1].as_expr(),u,extension=True)),s.Poly(pp[2].as_expr(),u,extension=True))
check('nonpositive_complex_shape_control_has_common_root',GG.degree()==1)
check('excluded_shape_control_is_not_real',s.simplify(s.im(-GG.nth(0)/GG.nth(1)))!=0)

profiles=[]
for n,gamma in [(6,s.Rational(501,100)),(6,60),(10,s.Rational(901,100)),(10,100),(12,s.Rational(1101,100)),(12,120)]:
    W=f.wronskian(n,gamma,z);pp=f.cyclic_sections(W,z,u)
    intervals=[p.intervals(eps=s.Rational(1,10**10)) for p in pp]
    neg=[int(p.count_roots(-s.oo,0)) for p in pp]
    for p,r in zip(pp,neg):check('parameter_edge_finite_real_rooted_sections',r==p.degree() and s.gcd(p,p.diff()).degree()==0)
    pair_results=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        tagged=[(iv,k) for k in [i,j] for iv,m in intervals[k]];tagged.sort(key=lambda v:sum(v[0]))
        alt=all(tagged[k][0][1]<tagged[k+1][0][0] and tagged[k][1]!=tagged[k+1][1] for k in range(len(tagged)-1))
        check('parameter_edge_finite_strict_interlacing',alt)
        pair_results.append(bool(alt))
    profiles.append({'n':n,'gamma':str(gamma),'degrees':[p.degree() for p in pp],'negative_root_counts':neg,'pairwise_strict_interlacing':pair_results})

U=z**4+11*z**3+111*z**2+1100*z+10000;V=(z*z+z+1)**2
for val in [U,V]:
    cc=s.Poly(val,z).all_coeffs()
    check('failed_induction_inputs_monic_positive_strict_logconcave',cc[0]==1 and all(c>0 for c in cc) and all(cc[k]**2>cc[k-1]*cc[k+1] for k in range(1,len(cc)-1)))
    parts=f.cyclic_sections(s.Poly(val,z),z,u)
    for p in parts:check('failed_induction_input_sections_real_negative',p.count_roots(-s.oo,0)==p.degree())
    root0=-parts[0].nth(0)/parts[0].nth(1);root1=-parts[1].nth(0)/parts[1].nth(1)
    check('failed_induction_input_sections_strictly_ordered',root1<root0<0)
output=f.next_fixed_shape(U,V,2,2,z);bad=f.cyclic_sections(output,z,u)[0]
check('failed_induction_output_section',bad.as_expr()==u*u+1772*u+900000)
check('failed_induction_output_nonreal_roots',s.discriminant(bad.as_expr(),u)==-460016)

audit={'status':'pass','canonical_proof':'GAP_THREE_WRONSKIAN_RECURRENCE.tex',
 'uniform_Wronskian_shape_shift_and_recurrence_proved':True,'bivariate_coefficient_positivity_proved':True,
 'no_common_real_section_root_proved':True,'one_real_rooted_section_suffices_proved':True,
 'fifth_order_all_positive_Gamma_shapes_proved':True,
 'uniform_section_real_rootedness_proved':False,'all_index_gap_three_SFC_closed':False,
 'separate_input_interlacing_induction_refuted':True,
 'general_FC_closed':False,'unrestricted_SFC_closed':False,'finite_extrapolation':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False,
 'boundary':'Recurrences and weaker conditional criterion all parameters; actual uniform simplicity only n<=5 at gap3, reusing the all-gap n<=4 theorem. Uniform higher-index real-section proof remains open.'}
diagnostics={'status':'exact','parameter_edge_profiles':profiles,'profiles_are_finite_evidence_only':True,
 'failed_separate_input_induction':{'U':str(U),'V':str(V),'a':10,'output_section':str(bad.as_expr()),'discriminant':-460016,'inputs_are_actual_common_Gamma_family':False,'SFC_counterexample':False},
 'fifth_discriminant_descending_coefficients':list(f.FIFTH_DISCRIMINANT_DESCENDING)}
for name,data in [('gap_three_wronskian_audit.json',audit),('gap_three_wronskian_diagnostics.json',diagnostics)]:
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
files=['GAP_THREE_WRONSKIAN_RECURRENCE.tex','gap_three_wronskian.py','verify_gap_three_wronskian.py','gap_three_wronskian_audit.json','gap_three_wronskian_diagnostics.json','higher_gap_adjoint_rank.py','HIGHER_GAP_ADJOINT_RANK.tex']
receipt=dict(audit,total=sum(counts.values()),counts=counts,files={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files})
(HERE/'gap_three_wronskian_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'total':receipt['total'],'counts':counts,'all_index_gap_three_SFC_closed':False},indent=2))
