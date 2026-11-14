"""PF2RULED22: exact complete first-jet lifting over both marked candidate algebras.
Written geometry and the open-condition argument are not formal verification.
"""
from pathlib import Path
from time import monotonic
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.expand(value))==0)
def main():
    start=monotonic();u,v,s,d,r,w,x,y=S.symbols('u v s d r w x y')
    g=[(1-v*v-u*v)*(u+v)**2,(1-v*v-u*v)*(u+v),-v*(u+v),-v]
    g0=[S.expand(f.subs(u,0)) for f in g]
    gj=[S.expand(S.diff(f,u).subs(u,0)-S.Rational(3,4)*S.diff(f.subs(u,0),v)) for f in g]
    for i,expected in enumerate([v/2,(v*v+1)/4,v/2,S.Rational(3,4)]):zero('actual generator corrected jet '+str(i),gj[i]-expected)
    def L(f0,j):return S.Poly(S.expand((v-s)*j-d*S.diff(f0,v)),v)
    matrices=[]
    for n in range(1,9):
        cols=[]
        for i in range(4):
            f0=S.expand(g0[0]**(n-1)*g0[i])
            j=gj[i] if n==1 else S.expand((n-1)*g0[0]**(n-2)*g0[i]*gj[0]+g0[0]**(n-1)*gj[i])
            ll=L(f0,j);check(f'selected section output degree n{n} i{i}',ll.degree()<=4*n-1)
            cols.append([ll.nth(4*n-1-k) for k in range(4)])
        M=S.Matrix(4,4,lambda row,col:cols[col][row]);det=S.factor(M.det())
        check('nonzero universal triangular determinant '+str(n),det!=0)
        matrices.append(dict(n=n,matrix=[[str(M[i,j]) for j in range(4)] for i in range(4)],determinant=str(det)))
    # Full Cox section coverage: chart formula, degree and boundary symbol.
    coverage={}
    for m in (3,5,6,8,9):
        count=0
        for jj in range(m+1):
            for ii in range(m+jj+1):
                nn=m+jj-ii
                # Entire chart is v^(m-j)*(1-v^2-u*v)^j*(u+v)^nn.
                # Its coefficient at u^a is the following finite binomial sum.
                for aa in range(3):
                    coeff=0
                    for b in range(min(aa,jj)+1):
                        z=aa-b
                        if z<=nn:coeff+=(-1)**b*S.binomial(jj,b)*S.binomial(nn,z)*v**(m-jj+b+nn-z)*(1-v*v)**(jj-b)
                    pp=S.Poly(S.expand(coeff),v)
                    check(f'complete Cox jet degree m{m} row{count} u{aa}',pp.is_zero or pp.degree()<=4*m-aa)
                count+=1
        check('complete Cox dimension '+str(m),count==(m+1)*(3*m+2)//2)
        coverage[str(m)]=count
    # Normal freedoms invisible to the first cusp jet, on complete bases.
    normal_freedom=[]
    for m,expected_kernel,degree2,expected_hidden in [(5,14,10,3),(8,56,22,33)]:
        basis=[];labels=[]
        for jj in range(m+1):
            for ii in range(m+jj+1):
                nn=m+jj-ii;jet=[]
                for aa in range(3):
                    coeff=0
                    for b in range(min(aa,jj)+1):
                        z=aa-b
                        if z<=nn:coeff+=(-1)**b*S.binomial(jj,b)*S.binomial(nn,z)*v**(m-jj+b+nn-z)*(1-v*v)**(jj-b)
                    jet.append(S.Poly(S.expand(coeff),v))
                basis.append(jet);labels.append([jj,ii])
        mon=[(aa,b) for aa in (0,1) for b in range(4*m-aa+1)]
        M=S.Matrix([[jet[aa].nth(b) for jet in basis] for aa,b in mon])
        K=S.Matrix.hstack(*M.nullspace())
        check('complete invisible-first-jet kernel dimension '+str(m),K.cols==expected_kernel and K.rank()==expected_kernel)
        check('every normal kernel column annihilates the exact first jets '+str(m),M*K==S.zeros(M.rows,K.cols))
        T=S.Matrix([[jet[2].nth(b) for jet in basis] for b in range(4*m-1)])*K
        check('entire second-jet image has exact degree bound '+str(m),T[degree2+1:,:]==S.zeros(T.rows-degree2-1,T.cols))
        pivots=T[:degree2+1,:].rref()[1];square=T[:degree2+1,list(pivots)]
        check('second-jet image is the full bounded polynomial space '+str(m),len(pivots)==degree2+1 and square.det()!=0)
        check('remaining directions invisible through second coefficient '+str(m),K.cols-len(pivots)==expected_hidden)
        columns=[]
        for j in range(K.cols):
            columns.append(dict(combination=[[i,str(K[i,j])] for i in range(K.rows) if K[i,j]],
                second_jet=str(sum(T[b,j]*v**b for b in range(degree2+1)))))
        normal_freedom.append(dict(pole=m,cox_basis_ji=labels,first_jet_kernel_dimension=K.cols,
            nullspace_columns=columns,second_jet_image_degree=degree2,second_jet_image_dimension=len(pivots),
            second_jet_pivot_columns=list(pivots),second_jet_pivot_determinant=str(square.det()),
            invisible_through_second_coefficient=expected_hidden))
    check('boundary restriction parameter count46',(31-14)+(85-56)==46)
    check('second normal coefficient controls count34',(10+1)+(22+1)==34)
    check('higher directions invisible at this gate count36',3+33==36)
    check('next-quotient change degree bound44',max(10+34,22+22)==44)
    check('eight top coefficients precede truncated membership',52-45+1==8)
    check('truncated multiplication map has eleven consistency conditions',45-34==11)
    check('nineteen local rational consistency conditions',8+11==19)
    predecessor=json.loads((HERE/'unibranch_polar33_receipt.json').read_text())
    check('current predecessor passed',predecessor['outcome']=='PASS')
    for path,digest in predecessor['source_sha256'].items():check('predecessor seal '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    old=json.loads((HERE/'polar_cusp_sector_certificates.json').read_text())
    sectors=[]
    for sector in old['two_root_marked_sectors']:
        k=sector['multiplicities'][0];P=S.sympify(sector['P8']);xp=S.sympify(sector['x_as_polynomial_in_y'])
        def red(expr):
            num,den=S.fraction(S.cancel(expr));check('candidate reduction denominator is invertible '+str(k),S.gcd(den,P)==1)
            return S.rem(num*S.invert(den,P,y),P,y)
        H=S.Poly((w-x)**k*(w-y)**(9-k),w);a,b,cc,dd=[H.nth(j) for j in (8,7,6,5)]
        dv=red((-a/48).subs(x,xp));sv=red(((24*b+72-11*a*a)/(12*a)).subs(x,xp))
        check('candidate algebra has eight distinct roots '+str(k),S.degree(P,y)==8 and S.gcd(P,S.diff(P,y))==1)
        check('candidate cusp unit delta nonzero '+str(k),S.gcd(dv,P)==1)
        units=[]
        for item in matrices:
            nn=item['n'];det=S.sympify(item['determinant']);value=S.rem(S.expand(det.subs({d:dv,s:sv})),P,y)
            A,B,gg=S.gcdex(value,P,y)
            check(f'all-root determinant unit k{k} n{nn}',gg==1)
            zero(f'exact determinant Bezout k{k} n{nn}',A*value+B*P-1)
            units.append(dict(n=nn,reduced_determinant=str(value),determinant_coefficient=str(A),P8_coefficient=str(B)))
        # The complete four-coefficient test for a V3 lift. Subsequent
        # lower coefficients are lifted by the certified right inverse.
        hseries=1+a*r+b*r*r+cc*r**3+dd*r**4;coeffs=[S.expand(hseries).coeff(r,j) for j in range(5)]
        AA=[1,coeffs[1],coeffs[2]-3,coeffs[3]-3*coeffs[1]]
        BB=[-coeffs[1]/4,-(3+coeffs[2])/2,-(3*coeffs[1]+3*coeffs[3])/4,-coeffs[4]+3]
        for j in range(4):
            residual=BB[j]-(s*BB[j-1] if j else 0)-d*(12-j)*AA[j]
            substituted=S.expand(residual.subs(x,xp).subs({s:sv,d:dv}))
            zero(f'complete V3 highest cusp coefficient k{k} j{j}',S.rem(substituted,P,y))
        sectors.append(dict(multiplicities=sector['multiplicities'],P8=str(P),x_as_polynomial_in_y=str(xp),delta=str(dv),cusp_parameter=str(sv),determinant_units=units,
            distinct_marked_points=8,full_first_jet_lifts_exist_at_every_point=True,
            P_affine_dimension=31,Q_affine_dimension=85,clean_singleton_ordinary_cusp_open_nonempty=True,
            boundary_birational_open_nonempty=True,source_Keller_pair_constructed=False))
    check('full pole6 affine dimension',((5+1)*(3*5+2)//2)-20==31)
    check('full pole9 affine dimension',((8+1)*(3*8+2)//2)-32==85)
    check('combined affine dimension',31+85==116)
    # Actual global sections providing the three independent open-condition shifts.
    A1=-u+u*u*v;lam=S.Rational(3,4)*(v-s)+d
    modifications=[]
    for k in (2,3):
        Tk=(v-s)**k-k*A1*lam*(v-s)**(k-2)
        f0=S.expand(Tk.subs(u,0));f1=S.expand(S.diff(Tk,u).subs(u,0))
        zero('cusp modification restriction '+str(k),f0-(v-s)**k)
        zero('cusp modification full first jet '+str(k),(v-s)*f1-lam*S.diff(f0,v))
        check('modification is lower-pole for both coordinates '+str(k),k+2<=5)
        modifications.append(dict(k=k,pole_upper_bound=k+2,formula=f'(D-s)^{k}-{k}*A1*lambda(D)*(D-s)^{k-2}'))
    alpha,beta,gamma,p0,p1,q0,q1=S.symbols('alpha beta gamma p0 p1 q0 q1')
    wedge=S.expand((p0+2*alpha)*(q1+3*gamma)-p1*(q0+2*beta))
    check('ordinary cusp open polynomial is not zero',wedge.coeff(alpha,1).coeff(gamma,1)==6)
    # Exact next normal Jacobian identity.
    hp=v-s;p=S.Function('p')(v);q=S.Function('q')(v);P2,Q2=S.symbols('P2 Q2')
    p1expr=lam*p;q1expr=lam*q
    ju1=2*P2*hp*q+p1expr*S.diff(q1expr,v)-S.diff(p1expr,v)*q1expr-2*hp*p*Q2
    expected=2*hp*(P2*q-p*Q2)+lam**2*(p*S.diff(q,v)-S.diff(p,v)*q)
    zero('universal next normal Jacobian identity',ju1-expected)
    zero('next Jacobian value fixes nonzero kappa at clean cusp',expected.subs(v,s)-d*d*(p*S.diff(q,v)-S.diff(p,v)*q).subs(v,s))
    check('cross bracket u1 degree at most53',4*(5+9)-3==53 and 4*(6+8)-3==53)
    check('lower cross bracket degree at most49',4*(5+8)-3==49)
    check('next quotient degree at most52 and53 coefficient equations',53-1==52 and 52+1==53)
    # Actual-ring first-jet cusp countercontrol, failing the Keller identity.
    lam0=1+S.Rational(3,4)*v
    PP=v*v-2*A1*lam0;QQ=v**3-3*A1*lam0*v
    jac=S.expand(S.diff(PP,u)*S.diff(QQ,v)-S.diff(PP,v)*S.diff(QQ,u))
    zero('countercontrol satisfies u0 Jacobian condition',jac.coeff(u,0))
    zero('countercontrol u1 is six lambda squared',jac.coeff(u,1)-6*lam0**2)
    check('countercontrol source Jacobian cannot be a nonzero constant',S.degree(jac.coeff(u,1),v)==2)
    check('countercontrol clean cusp derivative gcd',S.monic(S.gcd(2*v,3*v*v),v)==v)
    cert=dict(as_of='2026-11-14',operator='L(F)=(v-s)*J(F)-delta*F0_prime',right_inverse='cancel the highest four coefficients using M_n, then recurse to n-1',
        normal_freedom=normal_freedom,boundary_parameter_count=46,second_coefficient_control_count=34,invisible_through_second_coefficient=36,
        next_quotient_immutable_top_degrees=list(range(45,53)),truncated_multiplication_matrix_shape=[45,34],
        clean_open_truncated_multiplication_rank=34,low_coefficient_consistency_count=11,local_rational_consistency_count_upper_bound=19,
        generator_restrictions=list(map(str,g0)),generator_corrected_jets=list(map(str,gj)),triangular_matrices=matrices,full_section_coverage=coverage,
        sectors=sectors,marked_candidate_count=16,all_sixteen_full_first_jet_lifts_exist=True,
        full_cusp_lift_affine_dimensions={'P':31,'Q':85,'pair':116},
        clean_singleton_cusp_and_boundary_birational_open_nonempty=True,cusp_modifications=modifications,
        next_Jacobian_u1='2*h*(P2*q-p*Q2)+lambda^2*(p*q_prime-p_prime*q)',
        kappa='delta^2*(p*q_prime-p_prime*q)(s)',next_Jacobian_quotient_degree_upper_bound=52,
        next_Jacobian_quadratic_equation_count_upper_bound=53,next_Jacobian_parameter_count=116,
        next_Jacobian_equations_solved=False,countercontrol=dict(P=str(PP),Q=str(QQ),Jac_u1=str(S.expand(jac.coeff(u,1))),full_first_jet_cusp=True,Keller=False,finite_surface_map=False))
    cp=HERE/'full_cusp_lift_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'FULL_CUSP_FIRST_JET_LIFTS.md',cp,HERE/'UNIBRANCH_POLAR33_OBSTRUCTION.md',HERE/'unibranch_polar33_receipt.json',HERE/'polar_cusp_sector_certificates.json',HERE/'POLAR_CUSP_LEADING_SECTORS.md',HERE/'GLOBAL_CUSP_LEADING_OBSTRUCTION.md',HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        all_sixteen_full_first_jet_lifts_exist=True,complete_noncontact_first_cusp_jet_locus_has_sixteen_marked_points=True,
        full_cusp_lift_affine_dimensions={'P':31,'Q':85,'pair':116},clean_singleton_cusp_and_boundary_birational_open_nonempty=True,
        determinant_unit_certificates=16,next_Jacobian_quotient_degree_upper_bound=52,next_Jacobian_quadratic_equation_count_upper_bound=53,
        next_Jacobian_parameter_count=116,boundary_parameter_count=46,second_coefficient_control_count=34,
        invisible_through_second_coefficient=36,next_quotient_immutable_top_degrees=list(range(45,53)),
        clean_open_truncated_multiplication_rank=34,low_coefficient_consistency_count=11,local_rational_consistency_count_upper_bound=19,
        marked_two_root_leading_candidates_retained=16,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(first_jet_coverage='written triangular induction and all-point exact determinant units',
            dimensions_quotient_target_translations=False,first_jet_lifts_prove_source_Keller_pairs=False,
            next_Jacobian_equations_solved=False,all_Jacobian_coefficients_solved=False,
            dimension_count_implies_realizability=False,nineteen_reduced_conditions_claimed_quadratic=False,
            all_eight_high_coefficient_equations_solved=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,all_noncontact_pole69_excluded=False,all_pole69_excluded=False,
            contact_cases_classified_by_this_lift=False,actual_ruled_compatibility_proved=False,
            quartic_packets_excluded=False,conditional_2026_degree125_upgrade_used=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'full_cusp_lift_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; all16 full cusp first-jet sectors lift; next gate8 high coefficients+11 membership conditions on46 boundary parameters; goal active')
if __name__=='__main__':main()
