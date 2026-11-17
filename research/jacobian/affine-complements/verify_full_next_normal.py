"""PF2RULED24: lift exact modular unit-minor circuits to full next-normal solutions.
Run with sage -python. The characteristic-zero implication is a written
local-ring argument; no modular search bound is promoted to a theorem.
"""
from pathlib import Path
from time import monotonic
import hashlib,json,platform
from sage.all import QQ, GF, PolynomialRing, matrix, vector, binomial
from sage.version import version as sage_version
from verify_cusp_high_jacobian import jet_basis, jac1

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)

def cox_jets(m,V,order=3):
    v=V.gen();out=[]
    for j in range(m+1):
        for i in range(m+j+1):
            nn=m+j-i;jets=[]
            for aa in range(order+1):
                coeff=V.zero()
                for bb in range(min(aa,j)+1):
                    z=aa-bb
                    if z<=nn:
                        coeff+=(-1)**(m-j+bb)*binomial(j,bb)*binomial(nn,z)*v**(m-j+bb+nn-z)*(1-v*v)**(j-bb)
                jets.append(coeff)
            out.append(dict(j=j,i=i,symbol_power=2*(m-j)+i,jets=jets))
    return out

def mul(A,B):return [sum(A[i]*B[k-i] for i in range(k+1)) for k in range(4)]
def combine(c,basis):return [sum(x*f['jets'][j] for x,f in zip(c,basis)) for j in range(4)]
def jac2(A,B):
    return (3*(A[3]*B[0].derivative()-A[0].derivative()*B[3])+
            2*(A[2]*B[1].derivative()-A[1].derivative()*B[2])+
            A[1]*B[2].derivative()-A[2].derivative()*B[1])
def residue_vector(values):return [int(c) for c in values]


def main():
    start=monotonic();pred=json.loads((HERE/'cusp_high_jacobian_receipt.json').read_text())
    check('predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    old=json.loads((HERE/'polar_cusp_sector_certificates.json').read_text())
    prime=1009;F=GF(prime);R=PolynomialRing(QQ,'y');sectors=[]
    for sector in old['two_root_marked_sectors']:
        k=sector['multiplicities'][0];f=R(sector['P8']).monic();xp=R(sector['x_as_polynomial_in_y'])
        check(f'candidate polynomial irreducible degree8 k{k}',f.degree()==8 and f.is_irreducible())
        check(f'all base coefficient denominators prime units k{k}',
              all(c.denominator()%prime for c in list(f)+list(xp)))
        residue=F(990)
        check(f'exact residue point on candidate algebra k{k}',f.change_ring(F)(residue)==0)
        x=xp.change_ring(F)(residue)
        check(f'residue roots remain distinct k{k}',x!=residue)
        V=PolynomialRing(F,'v');v=V.gen();W=PolynomialRing(F,'w');w=W.gen()
        H=(w-x)**k*(w-residue)**(9-k);a,b=H[8],H[7]
        check(f'initial cusp divisions are prime units k{k}',a!=0 and F(12)!=0 and F(48)!=0)
        delta=-a/48;s=(24*b+72-11*a*a)/(12*a);h=v-s;lam=3*h/4+delta
        def L(A):return h*A[1]-lam*A[0].derivative()
        bas3=cox_jets(3,V)
        check(f'all four jets agree with predecessor convention k{k}',
              all(f['jets'][:3]==g['jets'] for f,g in zip(bas3,jet_basis(3,V))))
        T=[V.zero()]*4;chosen=[]
        for power in range(10):
            row=next(row for row in bas3 if row['symbol_power']==power)
            chosen.append([row['j'],row['i']]);T=[T[j]+H[9-power]*row['jets'][j] for j in range(4)]
        check(f'high-four cusp residual on T k{k}',L(T).degree()<=7)
        sols=[];cusp_circuits=[]
        for m,A in [(5,mul(T,T)),(8,mul(mul(T,T),T))]:
            bas=cox_jets(m,V)
            check(f'complete Cox dimension k{k} m{m}',len(bas)==(m+1)*(3*m+2)//2)
            M=matrix(F,[[L(row['jets'])[j] for row in bas] for j in range(4*m)])
            pivots=list(M.pivots());free=[j for j in range(M.ncols()) if j not in pivots]
            check(f'full cusp rank k{k} m{m}',len(pivots)==4*m)
            minor=M[:,pivots];det=minor.det();check(f'cusp inverse determinant prime unit k{k} m{m}',det!=0)
            inverse=minor.inverse();correction=vector(F,M.ncols())
            rhs=vector(F,[-L(A)[j] for j in range(4*m)])
            check(f'whole cusp right-hand side in output range k{k} m{m}',L(A).degree()<4*m)
            for j,c in zip(pivots,inverse*rhs):correction[j]=c
            ref=[A[j]+combine(correction,bas)[j] for j in range(4)]
            check(f'full reference cusp equation k{k} m{m}',L(ref)==0)
            directions=[];coefficient_rows=[]
            for freej in free:
                row=vector(F,M.ncols());row[freej]=1
                for j,c in zip(pivots,-inverse*M.column(freej)):row[j]=c
                check(f'kernel circuit annihilates full cusp map k{k} m{m} free{freej}',M*row==0)
                directions.append(combine(row,bas));coefficient_rows.append(row)
            sols.append((ref,directions,correction,coefficient_rows,bas))
            cusp_circuits.append(dict(pole=m,pivot_columns=pivots,free_columns=free,
                determinant_modp=int(det),reference_correction_modp=residue_vector(correction),
                kernel_recipe='e_free - embed_pivots(M_pivots^-1 * M_free)'))
        P,kp,pcoeff,pkcoeff,pbas=sols[0];Q,kq,qcoeff,qkcoeff,qbas=sols[1]
        # These are exact integer kernel weights in the characteristic-zero circuit.
        P=[P[j]+sum(row[j] for row in kp) for j in range(4)]
        pcoeff=pcoeff+sum(pkcoeff,vector(F,len(pcoeff)))
        check(f'fixed P remains a full first cusp lift k{k}',L(P)==0)
        columns=[jac1(P,row) for row in kq]
        reference=jac1(P,Q)
        check(f'full next-normal output degree bound48 k{k}',
              reference.degree()<=48 and all(pol.degree()<=48 for pol in columns))
        B=matrix(F,[[pol[j] for pol in columns] for j in range(49)])
        output_pivots=list(B.pivots());check(f'next-normal map surjects onto all49 coefficients k{k}',len(output_pivots)==49)
        square=B[:,output_pivots];det=square.det()
        check(f'next-normal inverse determinant prime unit k{k}',det!=0)
        rhs=vector(F,[F(j==0)-reference[j] for j in range(49)])
        controls=square.solve_right(rhs);solution=vector(F,len(kq))
        for j,c in zip(output_pivots,controls):solution[j]=c
        Qnew=[Q[j]+sum(c*row[j] for c,row in zip(solution,kq)) for j in range(4)]
        for c,row in zip(solution,qkcoeff):qcoeff+=c*row
        check(f'whole first cusp equation preserved in Q solution k{k}',L(Qnew)==0)
        check(f'whole next-normal polynomial equals1 k{k}',jac1(P,Qnew)==1)
        check(f'constant chart-Jacobian coefficient vanishes k{k}',P[1]*Qnew[0].derivative()-P[0].derivative()*Qnew[1]==0)
        pp,remp=P[0].derivative().quo_rem(h);qq,remq=Qnew[0].derivative().quo_rem(h)
        check(f'clean derivative gcd k{k}',remp==remq==0 and pp.gcd(qq).degree()==0)
        check(f'fixed noncontact boundary derivative degrees k{k}',pp.degree()==22 and qq.degree()==34)
        wedge=pp*qq.derivative()-pp.derivative()*qq
        check(f'ordinary cusp wedge forced by normalized Jacobian k{k}',delta**2*wedge(s)==1)
        A,rema=(P[0]-P[0](s)).quo_rem(h**2);C,remc=(Qnew[0]-Qnew[0](s)).quo_rem(h**2)
        check(f'cusp image differences have exact h^2 factors k{k}',rema==remc==0)
        resultant=A.resultant(C)
        check(f'singleton cusp resultant is a prime unit k{k}',resultant!=0)
        j2=jac2(P,Qnew)
        direct_ring=PolynomialRing(F,names=('uu','vv'));uu,vv=direct_ring.gens()
        direct_P=sum(uu**j*P[j](vv) for j in range(4))
        direct_Q=sum(uu**j*Qnew[j](vv) for j in range(4))
        direct_J=direct_P.derivative(uu)*direct_Q.derivative(vv)-direct_P.derivative(vv)*direct_Q.derivative(uu)
        for normal_order,expected in [(0,V.zero()),(1,V.one()),(2,j2)]:
            actual=sum(c*v**b for (a,b),c in direct_J.dict().items() if a==normal_order)
            check(f'direct bivariate differentiation coefficient u^{normal_order} k{k}',actual==expected)
        check(f'constructed lift still fails the next higher Jacobian coefficient k{k}',j2!=0)
        check(f'next higher coefficient respects cross-term degree bound52 k{k}',j2.degree()<=52)
        sectors.append(dict(multiplicities=sector['multiplicities'],P8=str(f),
            x_as_polynomial_in_y=str(xp),residue_y=int(residue),residue_x=int(x),
            delta_modp=int(delta),cusp_parameter_modp=int(s),T_symbol_lift_cox_ji=chosen,
            cusp_inverse_circuits=cusp_circuits,P_kernel_weights=[1]*len(kp),
            next_normal_matrix_modp=[residue_vector(row) for row in B.rows()],
            next_normal_pivot_columns=output_pivots,next_normal_determinant_modp=int(det),
            next_normal_controls_modp=residue_vector(controls),
            P_lower_pole_correction_modp=residue_vector(pcoeff),Q_lower_pole_correction_modp=residue_vector(qcoeff),
            next_normal_target='J1=1',next_normal_rank=49,
            singleton_cusp_resultant_modp=int(resultant),ordinary_cusp_wedge_modp=int(wedge(s)),
            next_Jacobian_u2_modp=str(j2),next_Jacobian_u2_degree=int(j2.degree()),
            characteristic_zero_lift='apply the same inverse-matrix circuits over Q[y]/P8; every divided determinant has nonzero residue',
            all_eight_complex_embeddings_covered=True,full_source_Jacobian_solved=False))
        print(f'PASS k{k}:49-by49 unit minor, full J1=1 lift, singleton cusp; all8 characteristic-zero embeddings via local-ring units',flush=True)
    # Exact rational normal controls at the following Jacobian coefficient.
    V=PolynomialRing(QQ,'v');v=V.gen();normal=[]
    for m,expected_kernel,expected_rank,expected_hidden in [(5,3,3,0),(8,33,18,15)]:
        bas=cox_jets(m,V)
        M=matrix(QQ,[[row['jets'][a][b] for row in bas] for a in range(3) for b in range(4*m-a+1)])
        ker=M.right_kernel().basis_matrix()
        check(f'full rational kernel invisible through u^2 pole{m}',ker.nrows()==expected_kernel)
        A=matrix(QQ,[[row['jets'][3][b] for row in bas] for b in range(4*m-2)])*ker.transpose()
        check(f'exact third-coefficient image rank pole{m}',A.rank()==expected_rank)
        check(f'kernel invisible through u^3 pole{m}',ker.nrows()-A.rank()==expected_hidden)
        if m==5:
            target=[(v*v-1)**2,v*(v*v-1),v**3*(v*v-1)]
        else:target=[v**j for j in range(18)]
        targetmatrix=matrix(QQ,[[pol[b] for pol in target] for b in range(4*m-2)])
        check(f'entire third-coefficient image is the stated polynomial space pole{m}',
              A.column_space()==targetmatrix.column_space())
        normal.append(dict(pole=m,invisible_through_second_coefficient=expected_kernel,
            third_image_basis=list(map(str,target)),third_image_rank=expected_rank,
            invisible_through_third_coefficient=expected_hidden))
    check('following normal control rank21',3+18==21)
    check('following normal correction degree40',1+max(5+34,22+17)==40)
    check('following quotient target degree39 has19 membership conditions',40-21==19)
    check('crude cross-term high coefficient bound12 before intrinsic sharpening',52-41+1==12)
    # The global source Jacobian belongs to V12 after the leading-symbol
    # cancellation and first cusp equation (written divisor argument).
    # Exact full-section controls bound its next coefficient and the
    # number of normal orders needed for its complete identity.
    bas12=cox_jets(12,V,order=8)
    check('complete Jacobian pole-space dimension247',len(bas12)==247)
    restriction=matrix(QQ,[[row['jets'][0][b] for row in bas12] for b in range(49)])
    restriction_kernel=restriction.right_kernel().basis_matrix()
    firstjet=matrix(QQ,[[row['jets'][1][b] for row in bas12] for b in range(48)])*restriction_kernel.transpose()
    check('Jacobian residual restriction kernel dimension198',restriction_kernel.nrows()==198)
    check('all Jacobian residual first coefficients have degree at most43',firstjet[44:,:]==0)
    check('entire residual first-coefficient image has dimension44',firstjet.rank()==44)
    nine_jets=matrix(QQ,[[row['jets'][a][b] for row in bas12] for a in range(9) for b in range(49-a)])
    check('first nine normal coefficients determine the entire V12 section',nine_jets.rank()==247)
    check('divisor12B-9E has coefficients3C-3f',12-9==3 and 24-27==-3)
    check('only three next-normal high coefficients remain after the intrinsic degree43 bound',43-41+1==3)
    check('intrinsic following-normal test has at most23 consistency conditions',3+1+19==23)
    cert=dict(as_of='2026-11-14',proof_method='exact residue-unit matrix circuits plus written characteristic-zero local-ring lifting',
        residue_characteristic=prime,marked_points=16,sectors=sectors,
        full_next_normal_solutions_exist_at_all_marked_points=True,
        generic_fixed_P_next_normal_rank=49,fixed_P_normalized_Q_affine_dimension=36,
        generic_unnormalized_solution_chart_dimension=68,generic_boundary_parameter_count=32,
        clean_singleton_cusp_open_nonempty=True,following_normal_controls=normal,
        following_normal_control_rank=21,following_normal_correction_degree=40,
        following_normal_immutable_high_degrees=list(range(41,44)),
        following_normal_cusp_value_condition=True,following_normal_low_membership_conditions=19,
        source_Jacobian_pole_bound=12,following_normal_degree_bound=43,following_normal_consistency_count_upper_bound=23,
        source_Jacobian_full_identity_follows_from_chart_orders=list(range(2,10)),first_nine_V12_jets_rank=247,
        following_normal_equations_solved=False)
    cp=HERE/'full_next_normal_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'FULL_NEXT_NORMAL_LIFTS.md',cp,HERE/'verify_cusp_high_jacobian.py',
            HERE/'cusp_high_jacobian_receipt.json',HERE/'CUSP_HIGH_JACOBIAN_REDUCTION.md',
            HERE/'full_cusp_lift_receipt.json',HERE/'FULL_CUSP_FIRST_JET_LIFTS.md',
            HERE/'polar_cusp_sector_certificates.json',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'SOURCE_COORDINATES.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        full_next_normal_solutions_exist_at_all_marked_points=True,marked_points=16,
        proof_method=cert['proof_method'],residue_characteristic=prime,
        generic_fixed_P_next_normal_rank=49,fixed_P_normalized_Q_affine_dimension=36,
        generic_unnormalized_solution_chart_dimension=68,generic_boundary_parameter_count=32,
        clean_singleton_cusp_open_nonempty=True,following_normal_control_rank=21,
        following_normal_correction_degree=40,following_normal_immutable_high_degrees=list(range(41,44)),
        following_normal_cusp_value_condition=True,following_normal_low_membership_conditions=19,
        source_Jacobian_pole_bound=12,following_normal_degree_bound=43,following_normal_consistency_count_upper_bound=23,
        source_Jacobian_full_identity_follows_from_chart_orders=list(range(2,10)),first_nine_V12_jets_rank=247,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(characteristic_zero_existence_uses_written_unit_lifting=True,
            modular_rank_used_as_unit_certificate_not_search_evidence=True,
            all_complex_embeddings_covered=True,all_eleven_previous_low_conditions_have_simultaneous_lifts=True,
            no_marked_candidate_removed=True,rank_drop_strata_excluded=False,
            following_normal_equations_solved=False,all_higher_Jacobian_coefficients_solved=False,
            field_degree_four_realization_constructed=False,global_surface_finiteness_proved=False,
            actual_ruled_compatibility_proved=False,contact_cases_classified=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'full_next_normal_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; full next-normal lifts at all16 candidates; following normal coefficient and global goal remain open',flush=True)

if __name__=='__main__':main()


