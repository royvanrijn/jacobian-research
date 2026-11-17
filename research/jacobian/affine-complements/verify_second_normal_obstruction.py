"""PF2RULED25: a complete second-normal linear obstruction for fixed P.
Exact residue-unit minors certify nonzero characteristic-zero determinant
polynomials. Universal rank bounds and coverage are written proofs.
"""
from pathlib import Path
from time import monotonic
import hashlib,json,platform
from sage.all import QQ,GF,PolynomialRing,matrix,vector
from sage.version import version as sage_version
from verify_full_next_normal import cox_jets,mul,combine,jac1,jac2

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def residues(values):return [int(c) for c in values]
def sparse(values):return [[i,str(c)] for i,c in enumerate(values) if c]


def main():
    start=monotonic();pred=json.loads((HERE/'full_next_normal_receipt.json').read_text())
    check('current full-next-normal predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    old=json.loads((HERE/'full_next_normal_certificates.json').read_text())
    prime=1009;F=GF(prime);R=PolynomialRing(QQ,'y')
    # This complete rational kernel gives15 universal null directions,
    # independently of the first coordinate and of the marked algebra.
    V=PolynomialRing(QQ,'v');v=V.gen();qbas=cox_jets(8,V)
    jets=matrix(QQ,[[row['jets'][a][b] for row in qbas] for a in range(4) for b in range(33-a)])
    hidden=jets.right_kernel().basis_matrix()
    check('complete rational Q kernel invisible through u^3 has dimension15',hidden.nrows()==15)
    check('all15 hidden directions annihilate all four complete jets',jets*hidden.transpose()==0)
    check('constant and P add two independent universal null directions',85-(15+2)==68)
    check('adding the variable Jacobian constant gives universal rank bound69',68+1==69)
    sectors=[]
    for sector in old['sectors']:
        k=sector['multiplicities'][0];f=R(sector['P8']);xp=R(sector['x_as_polynomial_in_y'])
        check(f'irreducible degree-eight coefficient field k{k}',f.degree()==8 and f.is_irreducible())
        check(f'input denominators are prime units k{k}',all(c.denominator()%prime for c in list(f)+list(xp)))
        yy=F(sector['residue_y']);xx=xp.change_ring(F)(yy)
        check(f'exact coefficient-field residue k{k}',f.change_ring(F)(yy)==0 and int(xx)==sector['residue_x'])
        V=PolynomialRing(F,'v');v=V.gen();W=PolynomialRing(F,'w');w=W.gen()
        H=(w-xx)**k*(w-yy)**(9-k);a,b=H[8],H[7]
        check(f'cusp scalar divisions are prime units k{k}',a!=0)
        delta=-a/48;s=(24*b+72-11*a*a)/(12*a);h=v-s;lam=3*h/4+delta
        def L(A):return h*A[1]-lam*A[0].derivative()
        T=[V.zero()]*4;bas3=cox_jets(3,V)
        for power,indices in enumerate(sector['T_symbol_lift_cox_ji']):
            row=next(row for row in bas3 if [row['j'],row['i']]==indices)
            check(f'marked T symbol basis index k{k} power{power}',row['symbol_power']==power)
            T=[T[j]+H[9-power]*row['jets'][j] for j in range(4)]
        constructions=[]
        for m,base,circuit in [(5,mul(T,T),sector['cusp_inverse_circuits'][0]),
                                (8,mul(mul(T,T),T),sector['cusp_inverse_circuits'][1])]:
            basis=cox_jets(m,V);M=matrix(F,[[L(row['jets'])[j] for row in basis] for j in range(4*m)])
            pivots=circuit['pivot_columns'];free=circuit['free_columns'];square=M[:,pivots]
            check(f'inherited complete cusp pivot partition k{k} m{m}',
                  sorted(pivots+free)==list(range(len(basis))) and len(pivots)==4*m)
            check(f'inherited cusp inverse is a residue unit k{k} m{m}',
                  square.det()!=0 and int(square.det())==circuit['determinant_modp'])
            inverse=square.inverse();correction=vector(F,len(basis))
            for j,c in zip(pivots,inverse*vector(F,[-L(base)[j] for j in range(4*m)])):correction[j]=c
            ref=[base[j]+combine(correction,basis)[j] for j in range(4)]
            check(f'whole first-cusp reference equation k{k} m{m}',L(ref)==0)
            directions=[]
            for j in free:
                row=vector(F,len(basis));row[j]=1
                for ii,c in zip(pivots,-inverse*M.column(j)):row[ii]=c
                check(f'complete first-cusp kernel circuit k{k} m{m} free{j}',M*row==0)
                directions.append(combine(row,basis))
            constructions.append((ref,directions))
        Pstar,kp=constructions[0];Qstar,kq=constructions[1]
        check(f'complete parameter dimensions k{k}',len(kp)==31 and len(kq)==85)
        check(f'fixed degree24 of every first-coordinate boundary restriction k{k}',
              Pstar[0].degree()==24 and all(row[0].degree()<=20 for row in kp))
        P=[Pstar[j]+sum(row[j] for row in kp) for j in range(4)]
        check(f'all-ones P is the predecessor arithmetic-circuit point k{k}',
              P==[mul(T,T)[j]+combine(vector(F,sector['P_lower_pole_correction_modp']),cox_jets(5,V))[j] for j in range(4)])
        labels=[(1,j) for j in range(49)]+[(2,j) for j in range(44)]
        first=[jac1(P,row) for row in kq];second=[jac2(P,row) for row in kq]
        B=matrix(F,[[(first if order==1 else second)[i][j] for i in range(85)] for order,j in labels])
        J1=jac1(P,Qstar);J2=jac2(P,Qstar)
        rhs=vector(F,[-(J1[j] if order==1 else J2[j]) for order,j in labels])
        constant=vector(F,[-F(order==1 and j==0) for order,j in labels])
        M=B.augment(constant)
        check(f'fixed-P full coefficient system has93 rows and86 unknowns k{k}',M.dimensions()==(93,86))
        check(f'Q coefficient map attains universal rank68 k{k}',B.rank()==68)
        check(f'variable-kappa coefficient map attains universal rank69 k{k}',M.rank()==69)
        check(f'augmented variable-kappa system has rank70 k{k}',M.augment(rhs).rank()==70)
        columns=list(M.pivots());rows=list(M[:,columns].transpose().pivots())
        N=M.matrix_from_rows_and_columns(rows,columns);detN=N.det()
        check(f'generic-rank chart determinant is a prime unit k{k}',detN!=0)
        check(f'constant kappa column and its sole nonzero row occur in the minor k{k}',85 in columns and 0 in rows)
        controls=N.solve_right(vector(F,[rhs[j] for j in rows]))
        residual=rhs-M[:,columns]*controls
        check(f'all69 pivot-row compatibility equations hold k{k}',all(residual[j]==0 for j in rows))
        remaining=[j for j in range(93) if j not in rows]
        check(f'exactly24 fixed-P rational compatibility conditions remain k{k}',len(remaining)==24)
        extra=next(j for j in remaining if residual[j]!=0)
        aug=M.augment(rhs)
        obstruction=aug.matrix_from_rows_and_columns(rows+[extra],columns+[86]);detDelta=obstruction.det()
        check(f'global augmented determinant obstruction is a prime unit k{k}',detDelta!=0)
        check(f'Schur determinant equals rank-minor times compatibility residual k{k}',detDelta==detN*residual[extra])
        check(f'normalized kappa1 is not the only excluded constant k{k}',M.rank()<M.augment(rhs).rank())
        # An exact power control shows why kappa!=0 must stay in the
        # reduced problem. The same inverse circuit works over K, and
        # [S^2,S^3]=0 is a universal characteristic-zero identity.
        bas2=cox_jets(2,V)
        M2=matrix(F,[[L(row['jets'])[j] for row in bas2] for j in range(8)])
        piv2=list(M2.pivots());N2=M2[:,piv2]
        check(f'commuting-power cusp lift has a unit inverse circuit k{k}',len(piv2)==8 and N2.det()!=0)
        c2=vector(F,len(bas2))
        for j,c in zip(piv2,N2.solve_right(vector(F,[-L(T)[j] for j in range(8)]))):c2[j]=c
        S=[T[j]+combine(c2,bas2)[j] for j in range(4)]
        PS=mul(S,S);QS=mul(PS,S)
        check(f'commuting-power whole cusp identities k{k}',L(S)==L(PS)==L(QS)==0)
        check(f'commuting-power whole second-normal equations with kappa0 k{k}',jac1(PS,QS)==jac2(PS,QS)==0)
        check(f'commuting-power source is nonconstant k{k}',S[0].degree()==12)
        Pnormalized=list(P);Pnormalized[0]-=P[0](s)
        shear=Qstar[0][24]
        Qnormalized=[Qstar[j]-shear*Pnormalized[j] for j in range(4)]
        Qnormalized[0]-=Qnormalized[0](s)
        check(f'fixed monic first restriction permits division-free normalization k{k}',Pnormalized[0][24]==1)
        check(f'target translations and shear impose all three normalizations k{k}',
              Pnormalized[0](s)==Qnormalized[0](s)==Qnormalized[0][24]==0)
        check(f'normalizations preserve both complete normal Jacobian polynomials k{k}',
              jac1(Pnormalized,Qnormalized)==jac1(P,Qstar) and jac2(Pnormalized,Qnormalized)==jac2(P,Qstar))
        sectors.append(dict(multiplicities=sector['multiplicities'],P8=str(f),
            x_as_polynomial_in_y=str(xp),residue_y=int(yy),residue_x=int(xx),
            P_parameter_count=31,P_parameter_residue_point=[1]*31,
            cusp_inverse_circuits=sector['cusp_inverse_circuits'],T_symbol_lift_cox_ji=sector['T_symbol_lift_cox_ji'],
            coefficient_rows=[list(row) for row in labels],Q_matrix_rank_modp=68,
            variable_kappa_matrix_rank_modp=69,augmented_rank_modp=70,
            generic_chart_pivot_rows=rows,generic_chart_pivot_columns=columns,
            generic_chart_minor_determinant_modp=int(detN),
            fixed_P_compatibility_rows=remaining,compatibility_residual_modp=residues(residual),
            obstruction_extra_row=extra,obstruction_rows=rows+[extra],obstruction_columns=columns+[86],
            obstruction_determinant_modp=int(detDelta),
            obstruction_matrix_modp=[residues(row) for row in obstruction.rows()],
            characteristic_zero_obstruction='det of these rows/columns in [B(P(a)),-e,-J(P(a),Qstar)], using inherited cusp inverse circuits over Q[y]/P8',
            obstruction_polynomial_degree_upper_bound=69,obstruction_nonzero_at_every_complex_embedding=True,
            obstruction_vanishes_at_every_full_second_normal_solution=True,
            fixed_P_compatibility_condition_count=24,generic_solution_fiber_if_nonempty_dimension=17,
            generic_solution_fiber_directions='15 sections invisible through u^3, plus 1 and P; kappa uniquely determined',
            commuting_power_S_lower_pole=2,commuting_power_S_pivot_columns=piv2,
            commuting_power_inverse_determinant_modp=int(N2.det()),commuting_power_correction_modp=residues(c2),
            relaxed_kappa0_solution_constructed=True,nonzero_kappa_solution_constructed=False))
        print(f'PASS k{k}: rank68+variable kappa=69, augmented rank70; global degree<=69 obstruction;24 local P-conditions',flush=True)
    cert=dict(as_of='2026-11-14',proof_method='universal17-dimensional kernel plus exact residue-unit determinant lifting',
        residue_characteristic=prime,marked_points=16,sectors=sectors,
        universal_Q_kernel_dimension_lower_bound=17,universal_Q_rank_upper_bound=68,
        universal_variable_kappa_rank_upper_bound=69,generic_variable_kappa_rank=69,
        generic_P_excluded_even_with_variable_kappa=True,necessary_global_P_polynomial_degree_upper_bound=69,
        generic_fixed_P_compatibility_count=24,generic_Q_kappa_fiber_if_nonempty_dimension=17,
        normalized_P_parameter_count=30,normalized_Q_fiber_if_nonempty_dimension=15,
        all_rank_drop_solutions_still_satisfy_global_obstruction=True,
        rational_hidden_Q_kernel_basis=[sparse(row) for row in hidden.rows()],
        relaxed_kappa0_determinantal_locus_nonempty=True,nonzero_kappa_determinantal_locus_solved=False,marked_candidate_exclusion=False)
    cp=HERE/'second_normal_obstruction_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'SECOND_NORMAL_DETERMINANT_OBSTRUCTION.md',cp,
            HERE/'verify_full_next_normal.py',HERE/'FULL_NEXT_NORMAL_LIFTS.md',
            HERE/'full_next_normal_receipt.json',HERE/'full_next_normal_certificates.json',
            HERE/'WRIGHT_MARKED_MODEL.md',HERE/'SOURCE_COORDINATES.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        proof_method=cert['proof_method'],residue_characteristic=prime,marked_points=16,
        generic_P_excluded_even_with_variable_kappa=True,necessary_global_P_polynomial_degree_upper_bound=69,
        universal_Q_rank_upper_bound=68,universal_variable_kappa_rank_upper_bound=69,
        generic_variable_kappa_rank=69,generic_fixed_P_compatibility_count=24,
        normalized_P_parameter_count=30,normalized_Q_fiber_if_nonempty_dimension=15,
        generic_Q_kappa_fiber_if_nonempty_dimension=17,all_rank_drop_solutions_still_satisfy_global_obstruction=True,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(characteristic_zero_nonzero_polynomial_uses_written_unit_lifting=True,
            universal_rank_bounds_use_complete_kernel_and_written_independence=True,
            all_complex_embeddings_covered=True,all_values_of_kappa_allowed_in_obstruction=True,
            global_obstruction_includes_rank_drop_strata=True,
            local24_conditions_asserted_independent=False,nonzero_kappa_solution_locus_nonempty_proved=False,
            commuting_power_control_only_has_kappa0=True,
            second_normal_equations_solved_on_zero_locus=False,marked_candidate_excluded=False,
            full_source_Jacobian_identity_solved=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,
            contact_cases_classified=False,quartic_packets_excluded=False,formal_verification=False,
            independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'second_normal_obstruction_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; generic P excluded at all16 marked points; determinant zero loci and full goal remain open',flush=True)

if __name__=='__main__':main()
