"""PF2RULED27: global polynomial parametrization and linear constant elimination.
No full degree-eight Groebner search or constant-remainder certificate search
is run. All algebraic identities are checked over the exact marked fields.
"""
from pathlib import Path
from time import monotonic
import gzip,hashlib,json,platform
from sage.all import QQ,NumberField,PolynomialRing
from sage.version import version as sage_version

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def polynomials(values):return list(map(str,values))

def main():
    start=monotonic();pred=json.loads((HERE/'square_boundary_gap_receipt.json').read_text())
    check('current boundary-gap predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    old=json.loads(gzip.decompress((HERE/'square_boundary_degree8_checkpoint.json.gz').read_bytes()))
    Y=PolynomialRing(QQ,'y');sectors=[]
    for sec in old['sectors']:
        k=sec['multiplicities'][0];p8=Y(sec['P8'])
        check(f'irreducible degree-eight marked field k{k}',p8.degree()==8 and p8.is_irreducible())
        K=NumberField(p8,'y');A=PolynomialRing(K,names=sec['variables'],order='degrevlex');a=A.gens()
        G=[A(g) for g in sec['high_groebner_basis']];inputs=[A(g) for g in sec['high_generators']]
        check(f'complete original high ideal has five coefficients and six GB generators k{k}',len(inputs)==5 and len(G)==6)
        for i,coeff in enumerate(sec['GB_in_input_coefficients']):
            check(f'input-to-GB exact source identity k{k} {i}',sum(A(c)*f for c,f in zip(coeff,inputs))==G[i])
        for i,coeff in enumerate(sec['input_in_GB_coefficients']):
            check(f'GB-to-input exact source identity k{k} {i}',sum(A(c)*f for c,f in zip(coeff,G))==inputs[i])
        g1,g2=G[1],G[2]
        m11,m12=g1.derivative(a[1]),g1.derivative(a[2]);m21,m22=g2.derivative(a[1]),g2.derivative(a[2])
        determinant=m11*m22-m12*m21
        check(f'global high-ideal elimination determinant is a nonzero FIELD scalar k{k}',determinant.degree()==0 and determinant!=0)
        r1=g1.subs({a[1]:0,a[2]:0});r2=g2.subs({a[1]:0,a[2]:0})
        sol={a[1]:(-r1*m22+r2*m12)/determinant,a[2]:(-r2*m11+r1*m21)/determinant}
        sol[a[0]]=-G[3].subs({a[0]:0}).subs(sol)
        sol[a[4]]=-G[4].subs({a[4]:0}).subs(sol)
        sol[a[5]]=-G[5].subs({a[5]:0}).subs(sol)
        free=[a[3],a[6],a[7],a[8],a[9]]
        check(f'five eliminated variables are polynomials only in five free variables k{k}',
              all(set(A(f).variables())<=set(free) for f in sol.values()))
        check(f'global polynomial graph has degree at most3 k{k}',max(A(f).degree() for f in sol.values())<=3)
        vectors={a[1]:[A.zero() for _ in G],a[2]:[A.zero() for _ in G]}
        vectors[a[1]][1]=m22/determinant;vectors[a[1]][2]=-m12/determinant
        vectors[a[2]][1]=-m21/determinant;vectors[a[2]][2]=m11/determinant
        vectors[a[0]]=[-G[3].derivative(a[1])*vectors[a[1]][j]-G[3].derivative(a[2])*vectors[a[2]][j]+A(j==3) for j in range(6)]
        vectors[a[4]]=[A(j==4) for j in range(6)];vectors[a[5]]=[A(j==5) for j in range(6)]
        for variable in sol:
            check(f'graph generator explicit membership in high ideal k{k} {variable}',
                  sum(c*g for c,g in zip(vectors[variable],G))==variable-sol[variable])
        check(f'all high equations vanish identically on the global polynomial graph k{k}',all(g.subs(sol)==0 for g in G))
        # First parameterize by q4..q7 and the old constant coordinate t0.
        B0=PolynomialRing(K,names=['q4','q5','q6','q7','t0'],order='degrevlex');q4,q5,q6,q7,t0=B0.gens()
        free_map=A.hom([0,0,0,t0,0,0,q4,q5,q6,q7],B0)
        images=[free_map(A(sol.get(v,v))) for v in a];phi=A.hom(images,B0)
        V0=PolynomialRing(B0,'v');v0=V0.gen();VA=PolynomialRing(A,'v')
        def transfer(text):return V0([phi(c) for c in VA(text)])
        root0=transfer(sec['reduced_root']);M0=transfer(sec['reduced_critical_remainder']);L0=transfer(sec['reduced_Rprime_over_h'])
        Wring=PolynomialRing(K,'w');w=Wring.gen();y=K.gen();x=K(Y(sec['x_as_polynomial_in_y']));H=(w-x)**k*(w-y)**(9-k)
        aa,bb=H[8],H[7];s=(24*bb+72-11*aa*aa)/(12*aa)
        mu=root0[0].derivative(t0)
        check(f'old free coordinate shifts S0 by a fixed nonzero FIELD scalar k{k}',mu.degree()==0 and mu!=0 and
              all(root0[j].derivative(t0)==0 for j in range(1,root0.degree()+1)))
        check(f'critical remainder and L are independent of root translation k{k}',
              all(c.derivative(t0)==0 for c in list(M0)+list(L0)))
        base_value=root0(s).subs({t0:0})
        B=PolynomialRing(K,names=['q4','q5','q6','q7','z'],order='degrevlex');q=B.gens();z=q[4]
        inject=B0.hom([*q[:4],0],B);beta=inject(base_value)
        normalize=B0.hom([*q[:4],(z-beta)/K(mu)],B)
        V=PolynomialRing(B,'v');v=V.gen();h=v-B(s)
        root=V([normalize(c) for c in root0]);M=V([normalize(c) for c in M0]);L=V([normalize(c) for c in L0])
        W=root-z;Ms=M(s)
        check(f'global constant normalization is exactly S0(s)=z k{k}',root(s)==z)
        check(f'W has no z and vanishes doubly at the cusp k{k}',W(s)==W.derivative()(s)==0 and all(c.derivative(z)==0 for c in W))
        check(f'normalized M and L retain all five-free-space points k{k}',M.degree()<=5 and L.degree()==6 and L.leading_coefficient()==8 and
              all(c.derivative(z)==0 for c in list(M)+list(L)))
        check(f'nonzero kappa is exactly nonzero Ms times z k{k}',normalize(phi(A(sec['forced_scaled_kappa'])))==2*Ms*z and Ms!=0)
        D=(Ms-M)//h
        quotient,E=(M*(W//h)).quo_rem(L)
        check(f'linear constant coefficient D has degree at most4 k{k}',D.degree()<=4)
        check(f'companion remainder E has degree at most5 k{k}',E.degree()<=5)
        check(f'all D,E coefficients are polynomials only in the four q variables k{k}',
              all(c.derivative(z)==0 for c in list(D)+list(E)))
        numerator=2*Ms*z-2*M*root
        residual=2*h*(D*z-E)
        check(f'full old divisibility remainder equals new six linear equations k{k}',
              numerator.quo_rem(h*L)[1]==residual)
        check(f'exact identity controls the quotient as well k{k}',numerator+2*h*L*quotient==residual)
        # Each nonzero D_i chart is certified algebraically without a search.
        chart_data=[]
        for i in range(5):
            Di,Ei=D[i],E[i]
            for j in range(5):
                check(f'chart cross-product identity k{k} i{i} j{j}',
                      Di*(D[j]*z-E[j])-D[j]*(Di*z-Ei)==D[j]*Ei-Di*E[j])
            chart_data.append(dict(index=i,previous_D_coefficients_zero=list(range(i)),pivot=f'D[{i}] != 0',
                recovered_root_constant=f'E[{i}]/D[{i}]',nonzero_condition=f'M(s)*E[{i}] != 0',
                equations=['E[5]=0']+[f'D[{i}]*E[{j}]-E[{i}]*D[{j}]=0' for j in range(5) if j!=i]))
        E_strings=polynomials(E)
        sectors.append(dict(multiplicities=sec['multiplicities'],P8=str(p8),x_as_polynomial_in_y=sec['x_as_polynomial_in_y'],
            old_variables=list(A.variable_names()),old_high_generators=polynomials(inputs),old_high_GB=polynomials(G),
            elimination_field_determinant=str(K(determinant)),elimination_field_inverse=str(1/K(determinant)),
            graph_eliminated_variable_maps={str(var):str(sol[var]) for var in sol},
            graph_generator_in_GB_coefficients={str(var):polynomials(vectors[var]) for var in sol},
            free_coordinates=['a3','c4','c5','c6','c7'],normalized_coordinates=list(B.variable_names()),
            root_translation_mu=str(K(mu)),root_translation_beta=str(beta),
            old_variable_maps_in_normalized_coordinates=polynomials([normalize(f) for f in images]),
            cusp_parameter=str(s),root_deviation_W=str(W),critical_remainder_M=str(M),Rprime_over_h_L=str(L),
            constant_remainder_Ms=str(Ms),linear_constant_D=str(D),
            E_arithmetic_circuit='remainder_v(M*(W/(v-s)), L)',
            E_coefficient_sha256=[hashlib.sha256(text.encode()).hexdigest() for text in E_strings],
            E_coefficient_degrees=[int(c.degree()) for c in E],E_coefficient_monomial_counts=[len(c.monomials()) for c in E],
            reduced_equations=['E[5]=0']+[f'D[{i}]*z-E[{i}]=0' for i in range(5)],
            nonzero_condition='M(s)*z != 0',saturated_equation='t*M(s)*z-1=0',
            nonzero_D_charts=chart_data,
            zero_D_stratum=dict(equations=['D[0]=...=D[4]=0','E[0]=...=E[5]=0'],
                nonzero_condition='M(s) != 0',remaining_root_constant='any z != 0',excluded=False),
            full_nonzero_kappa_locus_solved=False,high_locus_is_global_affine_five_space=True))
        print(f'PASS k{k}: global A5 high locus, six linear z equations, five pivot charts plus complete zero-D stratum',flush=True)
    cert=dict(as_of='2026-11-14',marked_points=16,sectors=sectors,
        high_locus_global_polynomial_isomorphism='A5',remaining_core_parameter_count=4,
        full_linearized_divisibility_equation_count=6,additional_root_constant_count=1,
        nonzero_D_chart_count=5,zero_D_stratum_retained=True,
        full_nonzero_kappa_degree8_ideal_solved=False,
        proof_method='exact number-field polynomial graph identities and monic univariate division; no full elimination')
    cp=HERE/'degree8_affine_reduction_certificates.json.gz'
    cp.write_bytes(gzip.compress((json.dumps(cert,indent=2)+'\n').encode(),mtime=0))
    inputs=[Path(__file__),HERE/'DEGREE_EIGHT_AFFINE_REDUCTION.md',cp,
            HERE/'square_boundary_gap_receipt.json',HERE/'square_boundary_degree8_checkpoint.json.gz',
            HERE/'SQUARE_BOUNDARY_DEGREE_GAP.md',HERE/'SECOND_NORMAL_DETERMINANT_OBSTRUCTION.md',
            HERE/'WRIGHT_MARKED_MODEL.md',HERE/'sources/degree8_constant_remainder_unit_attempt.py.txt']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,marked_points=16,
        high_locus_global_polynomial_isomorphism='A5',remaining_core_parameter_count=4,
        full_linearized_divisibility_equation_count=6,nonzero_D_chart_count=5,zero_D_stratum_retained=True,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(all_complex_embeddings_covered=True,global_parameterization_has_no_variable_denominator=True,
            high_ideal_isomorphism_replayed=True,old_divisibility_ideal_preserved=True,
            nonzero_kappa_condition_retained=True,all_zero_pivot_strata_retained=True,
            zero_D_stratum_excluded=False,full_degree8_nonzero_kappa_locus_solved=False,
            general_determinant_locus_solved=False,actual_global_root_hypothesis_removed=False,
            whole_marked_candidate_excluded=False,source_Keller_pair_constructed=False,
            field_degree_four_realization_constructed=False,global_surface_finiteness_proved=False,
            actual_ruled_compatibility_proved=False,contact_cases_classified=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'degree8_affine_reduction_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; exact degree8 reduction retains all special strata; full goal active',flush=True)

if __name__=='__main__':main()
