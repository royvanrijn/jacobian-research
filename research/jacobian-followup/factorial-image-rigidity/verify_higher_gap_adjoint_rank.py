"""Exact replay of general adjoints, both rank outcomes, and failed Hurwitz shortcut."""
from pathlib import Path
import hashlib,json
import sympy as s
import higher_gap_adjoint_rank as f
import cubic_adjoint_stability as cubic
HERE=Path(__file__).resolve().parent
counts={}
def check(key,condition):
    assert condition,key
    counts[key]=counts.get(key,0)+1

z,x,t,n,M,g,d,q=s.symbols('z x t n M g d q')
u=1-x;A=u*(q-1);Aprime=1-(d+1)*q
b=(M-1)*(q-1)+d*n
logC=g/u+(n+1)*d*q/(u*(q-1))+(n+1)/x
check('universal_residue_generating_identity',s.cancel(((n+1)*A/x-Aprime+b-A*logC).subs(M,g+d*n))==0)
# Full Euler coefficient recurrence and derivative-degree transformation.
ell=s.symbols('ell')
check('universal_Euler_degree_reversal',s.expand(d*((n-ell/d)*z**d+ell*q/d)-(ell*q-z**d*(ell-d*n)))==0)
# Verify the adjoint identity symbolically, with arbitrary n,M and arbitrary B.
B=s.Function('B')(z)
for gap in range(2,5):
    v=s.Function('v')(z)
    L=f.scalar_operator(v,n,M,gap,z)
    coeff=[s.expand(L).coeff(v if j==0 else s.diff(v,z,j)) for j in range(gap+2)]
    check('monic_scalar_equation',s.simplify(coeff[-1]-1)==0)
    check('next_Euler_coefficient',s.simplify(coeff[-2]-gap*(1-M)/z)==0)
    rr=-1+(1-M)/z
    gaugeD=lambda val:s.diff(val,z)+rr*val
    adj=0
    for j,coef in enumerate(coeff):
        val=coef*B
        for _ in range(j):val=gaugeD(val)
        adj+=(-1)**j*val
    check('arbitrary_parameter_adjoint_gauge',s.simplify(adj-f.polynomial_operator(B,n,M,gap,z)/z)==0)

for gap in range(2,6):
    for N in range(1,4):
        BB=f.auxiliary(N,g,gap,z).as_expr()
        C=s.series((1-x)**(-g)*((1-(1-x)**gap)/(gap*x))**(-N-1),x,0,N+1).removeO().expand()
        expected=sum(s.ff(N,j)*C.coeff(x,j)*z**(N-j) for j in range(N+1))
        check('generating_series_and_recurrence',s.expand(BB-expected)==0)
        check('symbolic_gamma_auxiliary_ODE',s.expand(f.polynomial_operator(BB,N,g+gap*N,gap,z))==0)
        TT=f.tilted(N,g,gap,t).as_expr();theta=t*s.diff(TT,t)
        check('all_positive_shape_current_bridge',s.simplify(f.tilted(N,g+1,gap,t).as_expr()-TT-gap*theta/g)==0)
        check('all_positive_shape_previous_bridge',s.simplify(f.tilted(N-1,g,gap,t).as_expr()-TT+theta/N)==0)
        if gap==2:check('retained_cubic_auxiliary_recovered',s.expand(BB-cubic.auxiliary(N,g,z).as_expr())==0)
        for gamma in [s.Rational(1,2),s.Integer(N)]:
            R=f.reversed_tilted(N,gamma,gap,z).as_expr()
            check('ordinary_point_reversed_ODE',f.scalar_operator(R,N,gamma+gap*N,gap,z)==0)
        if gap<=4:
            gamma=s.Integer(N);BB=f.auxiliary(N,gamma,gap,z).as_expr();R=f.reversed_tilted(N,gamma,gap,z).as_expr()
            roots=f.roots_of_unity(gap)
            norm=s.prod(w**N for w in roots)*s.prod(roots[j]-roots[i] for i in range(gap) for j in range(i+1,gap))
            det=f.adjoint_matrix(BB,gap,z,full=True).det(method='domain-ge')
            check('full_adjoint_Wronskian_reconstruction',s.simplify(s.expand(det-norm*R))==0)

profiles=[]
for N,gamma,gap in [(2,s.Rational(1,2),3),(3,3,3),(4,4,3),(2,2,4),
                     (2,s.Rational(-3,2),2),(2,(-15+s.sqrt(-15))/6,3)]:
    BB=f.auxiliary(N,gamma,gap,z).as_expr();R=f.reversed_tilted(N,gamma,gap,z)
    G=f.exact_rank_gcd(BB,gap,z)
    # Zero itself is excluded from the rank theorem; remove it before comparing.
    while G.nth(0)==0:G=G.exquo(s.Poly(z,z,domain=G.domain))
    RR=s.Poly(R.as_expr(),z,extension=True);bad=s.gcd(RR,RR.diff())
    check('rank_gcd_detects_both_squarefree_and_multiple',G.degree()==bad.degree())
    if gamma in [s.Rational(-3,2),(-15+s.sqrt(-15))/6]:check('genuine_degenerate_outside_positive_shapes',G.degree()==gap)
    else:check('positive_shape_rank_controls',G.degree()==0)
    profiles.append({'n':N,'gamma':str(gamma),'gap':gap,'rank_minor_gcd_degree':G.degree(),'repeated_root_gcd_degree':bad.degree(),'inside_positive_Gamma_scope':s.sympify(gamma).is_positive is True})

BB=f.auxiliary(4,4,3,z).as_expr();W=f.gap_three_W(BB,z)
check('failed_Hurwitz_example_B',BB==z**4+36*z**3+520*z**2+3560*z+9720)
expected=[1,36,740,10560,113800,933920,5768000,24494400,59875200]
check('failed_Hurwitz_example_W',W.all_coeffs()==expected)
cc=W.all_coeffs();H=s.Matrix(8,8,lambda i,j:cc[2*j-i+1] if 0<=2*j-i+1<=8 else 0)
hurwitz=[H[:j,:j].det(method='domain-ge') for j in range(1,9)]
check('exact_negative_seventh_Hurwitz_determinant',hurwitz[6]==-169684039754781133897728000000)
u=s.symbols('u');parts=f.cyclic_parts(W,3,z,u);G=s.gcd(s.gcd(parts[0],parts[1]),parts[2])
check('failed_sufficient_target_not_SFC_failure',G.degree()==0)
# Replay one nontrivial complex parameter where the cubic orbit really is repeated.
gamma=(-15+s.sqrt(-15))/6;BB=f.auxiliary(2,gamma,3,z).as_expr();Wbad=f.gap_three_W(BB,z)
partsbad=f.cyclic_parts(Wbad,3,z,u)
partsgcd=s.gcd(s.gcd(s.Poly(partsbad[0].as_expr(),u,extension=True),s.Poly(partsbad[1].as_expr(),u,extension=True)),s.Poly(partsbad[2].as_expr(),u,extension=True))
check('cyclic_W_rank_detects_actual_multiple_root',partsgcd.degree()==1)

section_controls=[]
for N in [4,6,8]:
    WW=f.gap_three_W(f.auxiliary(N,N,3,z).as_expr(),z)
    pp=f.cyclic_parts(WW,3,z,u)
    intervals=[p.intervals(eps=s.Rational(1,10**8)) for p in pp]
    root_counts=[int(p.count_roots(-s.oo,0)) for p in pp]
    for p,root_count in zip(pp,root_counts):
        check('finite_section_simple_negative_roots',root_count==p.degree() and s.gcd(p,p.diff()).degree()==0)
    pair_controls=[]
    for i in range(3):
        for j in range(i+1,3):
            tagged=[(interval,idx) for idx in [i,j] for interval,mult in intervals[idx]]
            tagged.sort(key=lambda row:sum(row[0]))
            strict=all(tagged[k][0][1]<tagged[k+1][0][0] and tagged[k][1]!=tagged[k+1][1] for k in range(len(tagged)-1))
            check('finite_section_strict_interlacing',strict)
            pair_controls.append({'sections':[i,j],'strict_interlacing':bool(strict)})
    section_controls.append({'n':N,'gamma':N,'gap':3,'degrees':[p.degree() for p in pp],
       'negative_real_root_counts':root_counts,'pairs':pair_controls,
       'isolating_intervals':[[[str(a),str(b)] for (a,b),mult in row] for row in intervals]})

aa,bb,uu=s.symbols('aa bb uu')
check('universal_positive_third_central_moment_bound',s.expand(aa*bb*(2*bb-aa+uu)-3*aa**2*bb+2*aa**3-2*aa*(bb-aa)**2-aa*bb*uu)==0)
ss,tt=s.symbols('ss tt')
mu1=ss/2;mu2=(ss**2+2*tt)/6;mu3=ss*tt/2
check('universal_quartic_square_forces_zero_third_central_moment',s.expand(mu3-3*mu1*mu2+2*mu1**3)==0)
for gap in range(1,7):
    rr=s.Poly(s.rf(x,gap),x)
    check('rising_factorial_increasing_convex_controls',all(c>=0 for c in rr.diff().all_coeffs()) and any(c>0 for c in rr.diff().all_coeffs()) and all(c>=0 for c in rr.diff().diff().all_coeffs()))

audit={'status':'pass','canonical_proof':'HIGHER_GAP_ADJOINT_RANK.tex',
 'uniform_higher_gap_adjoint_proved':True,'all_gap_initial_four_orders_all_positive_shapes_proved':True,'all_positive_shape_lower_one_bridge_proved':True,
 'exact_repeated_root_rank_equivalence_proved':True,'full_Wronskian_reconstruction_proved':True,
 'gap_three_closed_sector_sufficiency_proved':True,'gap_three_sector_zero_freeness_proved':False,
 'Hurwitz_shortcut_refuted':True,'higher_gap_all_index_SFC_closed':False,
 'general_FC_closed':False,'unrestricted_SFC_closed':False,'finite_extrapolation':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False,
 'boundary':'All d>=2 adjoints and rank equivalence. Uniform rank/sector exclusion at d>=3 remains open. Complex/negative parameter failures are outside positive Gamma scope.'}
diagnostics={'status':'exact','profiles':profiles,'failed_Hurwitz_shortcut':{'n':4,'gamma':4,'gap':3,'B':str(f.auxiliary(4,4,3,z).as_expr()),'W_descending_coefficients':expected,'Hurwitz_leading_determinants':[str(v) for v in hurwitz],'cyclic_parts':[str(p.as_expr()) for p in parts],'cyclic_gcd_degree':G.degree(),'SFC_counterexample':False},'finite_section_controls':section_controls,'uniform_section_interlacing_proved':False,'floating_root_angles_not_used_as_evidence':True}
for name,data in [('higher_gap_adjoint_rank_audit.json',audit),('higher_gap_adjoint_rank_diagnostics.json',diagnostics)]:
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
files=['HIGHER_GAP_ADJOINT_RANK.tex','higher_gap_adjoint_rank.py','verify_higher_gap_adjoint_rank.py','higher_gap_adjoint_rank_audit.json','higher_gap_adjoint_rank_diagnostics.json']
receipt=dict(audit,counts=counts,total=sum(counts.values()),files={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files})
(HERE/'higher_gap_adjoint_rank_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'total':receipt['total'],'counts':counts,'higher_gap_all_index_SFC_closed':False},indent=2))
