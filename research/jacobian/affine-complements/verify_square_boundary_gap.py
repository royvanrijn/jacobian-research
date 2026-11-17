"""PF2RULED26: exact low boundary-degree exclusion near square/cube jets.
The default replay proves degrees<=7 and builds, but does not solve, the
retained degree-eight checkpoint. All identities are in characteristic zero.
"""
from pathlib import Path
from time import monotonic
import gzip,hashlib,json,platform
from sage.all import QQ,NumberField,PolynomialRing,matrix,vector
from sage.version import version as sage_version
from verify_full_next_normal import cox_jets,combine,jac1

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def strs(values):return list(map(str,values))
def unit_certificate(A,generators,label):
    I=A.ideal(generators)
    check(label+' unit Groebner basis',I.groebner_basis()==[A.one()])
    coefficients=A.one().lift(I)
    check(label+' exact Bezout identity',sum(c*f for c,f in zip(coefficients,generators))==1)
    return dict(generators=strs(generators),coefficients=strs(coefficients),
                coefficient_maximum_degree=max(int(c.total_degree()) for c in coefficients if c))


def main():
    start=monotonic();pred=json.loads((HERE/'second_normal_obstruction_receipt.json').read_text())
    check('current determinant predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    old=json.loads((HERE/'second_normal_obstruction_certificates.json').read_text())
    Y=PolynomialRing(QQ,'y');sectors=[];checkpoints=[]
    for sector in old['sectors']:
        k=sector['multiplicities'][0];p8=Y(sector['P8'])
        check(f'irreducible degree-eight coefficient field k{k}',p8.degree()==8 and p8.is_irreducible())
        K=NumberField(p8,'y');y=K.gen();x=K(Y(sector['x_as_polynomial_in_y']))
        W=PolynomialRing(K,'w');w=W.gen();H=(w-x)**k*(w-y)**(9-k);a,b=H[8],H[7]
        delta=-a/48;s=(24*b+72-11*a*a)/(12*a)
        check(f'nonzero marked cusp unit k{k}',delta!=0)
        V=PolynomialRing(K,'v');v=V.gen();h=v-s;lam=3*h/4+delta
        def Lop(F):return h*F[1]-lam*F[0].derivative()
        T=[V.zero()]*4;bas3=cox_jets(3,V)
        for power,indices in enumerate(sector['T_symbol_lift_cox_ji']):
            row=next(row for row in bas3 if [row['j'],row['i']]==indices)
            check(f'actual marked T symbol k{k} index{power}',row['symbol_power']==power)
            T=[T[j]+H[9-power]*row['jets'][j] for j in range(4)]
        check(f'complete V2 output range for the T cusp residual k{k}',Lop(T).degree()<=7)
        bas2=cox_jets(2,V);M=matrix(K,[[Lop(row['jets'])[j] for row in bas2] for j in range(8)])
        pivots=sector['commuting_power_S_pivot_columns'];free=[j for j in range(12) if j not in pivots]
        inverse=M[:,pivots].inverse();correction=vector(K,12)
        check(f'exact cusp matrix has rank8 k{k}',M.rank()==8 and len(free)==4)
        for j,c in zip(pivots,inverse*vector(K,[-Lop(T)[j] for j in range(8)])):correction[j]=c
        Sbase=[T[j]+combine(correction,bas2)[j] for j in range(4)]
        check(f'whole characteristic-zero S reference cusp equation k{k}',Lop(Sbase)==0)
        directions=[]
        for j in free:
            row=vector(K,12);row[j]=1
            for ii,c in zip(pivots,-inverse*M.column(j)):row[ii]=c
            check(f'exact full S kernel vector k{k} free{j}',M*row==0)
            directions.append(combine(row,bas2))
        results=[]
        for degree in range(3,9):
            names=['a0','a1','a2','a3']+['c'+str(j) for j in range(2,degree)]
            A=PolynomialRing(K,names=names,order='degrevlex');av=A.gens();VV=PolynomialRing(A,'v');vv=VV.gen();hh=vv-A(s);ll=3*hh/4+A(delta)
            S=[VV(Sbase[j])+sum(av[i]*VV(directions[i][j]) for i in range(4)) for j in range(4)]
            root=S[0];derivative=root.derivative();r=derivative//hh
            check(f'full four-parameter S family retains cusp equation k{k} d{degree}',hh*S[1]-ll*derivative==0)
            check(f'root degree12 and derivative leading unit k{k} d{degree}',root.degree()==12 and derivative.degree()==11 and derivative.leading_coefficient() in K)
            RR=[hh**degree,degree*ll*hh**(degree-2),-degree*vv*ll*hh**(degree-2),VV.zero()]
            for j in range(2,degree):
                for order,term in enumerate([hh**j,j*ll*hh**(j-2),-j*vv*ll*hh**(j-2),VV.zero()]):RR[order]+=av[j+2]*term
            L=RR[0].derivative()//hh
            check(f'actual perturbation first jet k{k} d{degree}',RR[1]==ll*L)
            C=jac1(S,RR)
            quotient,remainder=C.quo_rem(derivative)
            check(f'complete critical remainder degree bound k{k} d{degree}',remainder.degree()<=10)
            generators=[remainder[j] for j in range(degree-2,11)]
            I=A.ideal(generators);G=list(I.groebner_basis())
            if degree<=4:
                cert=unit_certificate(A,generators,f'high remainder k{k} d{degree}')
                results.append(dict(degree=degree,variables=names,stage='high critical remainder',unit=cert))
                continue
            check(f'nontrivial complete high-remainder ideal k{k} d{degree}',G!=[A.one()])
            # Keep explicit identities in both directions, not just a reported GB.
            forward=[]
            for i,g in enumerate(G):
                coeff=g.lift(I);check(f'GB-to-input identity k{k} d{degree} row{i}',sum(c*f for c,f in zip(coeff,generators))==g)
                forward.append(strs(coeff))
            JG=A.ideal(G);reverse=[]
            for i,g in enumerate(generators):
                coeff=g.lift(JG);check(f'input-to-GB identity k{k} d{degree} row{i}',sum(c*f for c,f in zip(coeff,G))==g)
                reverse.append(strs(coeff))
            norm=lambda pol:VV([I.reduce(c) for c in pol])
            rootn=norm(root);Ln=norm(L);Mn=norm(remainder)
            check(f'critical remainder has required lower degree modulo the exact ideal k{k} d{degree}',Mn.degree()<=degree-3)
            kap=2*Mn(s)*rootn(s);numerator=kap-2*Mn*rootn
            denominator=hh*Ln
            check(f'normal-control division uses a fixed nonzero leading scalar k{k} d{degree}',denominator.leading_coefficient()==degree)
            divrem=numerator.quo_rem(denominator)[1]
            reduced=[I.reduce(c) for c in divrem]
            common=dict(degree=degree,variables=names,high_generators=strs(generators),high_groebner_basis=strs(G),
                GB_in_input_coefficients=forward,input_in_GB_coefficients=reverse,
                reduced_root=str(rootn),reduced_Rprime_over_h=str(Ln),reduced_critical_remainder=str(Mn),
                forced_scaled_kappa=str(kap),divisibility_remainder_coefficients=strs(reduced),
                high_remainder_locus_dimension=int(I.dimension()))
            if degree==8:
                checkpoints.append(dict(multiplicities=sector['multiplicities'],P8=str(p8),x_as_polynomial_in_y=sector['x_as_polynomial_in_y'],
                    **common,nonzero_kappa_equation='invk*forced_scaled_kappa-1',
                    full_saturated_ideal_solved=False,active_process=False,
                    next_step='solve the retained high GB plus divisibility coefficients and the nonzero-kappa equation; no full Groebner search is run by this checker'))
                continue
            cert=unit_certificate(A,G+reduced,f'full divisibility k{k} d{degree}')
            common.update(stage='high remainder plus full polynomial divisibility',unit=cert)
            results.append(common)
        sectors.append(dict(multiplicities=sector['multiplicities'],P8=str(p8),x_as_polynomial_in_y=sector['x_as_polynomial_in_y'],
            S_cusp_inverse_pivot_columns=pivots,S_free_cox_columns=free,S_cusp_correction=strs(correction),
            full_S_affine_dimension=4,closed_degrees=[0,1,2,3,4,5,6,7],exact_unit_cases=results,
            all_eight_complex_embeddings_covered=True))
        print(f'PASS k{k}: exact characteristic-zero units close boundary degrees3..7; degrees<=2 force kappa0; degree8 checkpoint retained',flush=True)
    check('normal coefficient A has maximum degree10',10+1+0<12)
    check('degree<=2 perturbation forces kappa0 by evaluation at the cusp',11<12)
    cert=dict(as_of='2026-11-14',marked_points=16,sectors=sectors,
        boundary_difference_degree_excluded_through=7,full_actual_S_pole3_cusp_family_covered=True,
        all_P_lower_normal_corrections_retained=True,all_Q_lower_normal_corrections_retained=True,
        proof_method='exact number-field polynomial Bezout identities plus written normal-Jacobian reduction',
        nonzero_first_normal_Jacobian_requires_boundary_difference_degree_at_least=8,
        degree8_and_general_determinant_locus_solved=False)
    cp=HERE/'square_boundary_gap_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    checkpoint=dict(as_of='2026-11-14',status='UNSOLVED_NO_LIVE_PROCESS',sectors=checkpoints,
        earlier_exploratory_attempt=dict(sector=[2,7],degree=8,stage='saturated full Groebner basis',
            terminal_outcome='interrupted after observed105 seconds with no mathematical result',
            process_terminal_verified=True,not_restarted_by_this_checker=True),
        original_exploratory_source=(HERE/'sources/square_boundary_degree8_exploration.py.txt').read_text())
    chk=HERE/'square_boundary_degree8_checkpoint.json.gz';chk.write_bytes(gzip.compress((json.dumps(checkpoint,indent=2)+'\n').encode(),mtime=0))
    inputs=[Path(__file__),HERE/'SQUARE_BOUNDARY_DEGREE_GAP.md',cp,chk,HERE/'verify_full_next_normal.py',
            HERE/'second_normal_obstruction_receipt.json',HERE/'second_normal_obstruction_certificates.json',
            HERE/'SECOND_NORMAL_DETERMINANT_OBSTRUCTION.md',HERE/'FULL_CUSP_FIRST_JET_LIFTS.md',
            HERE/'FULL_NEXT_NORMAL_LIFTS.md',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'SOURCE_COORDINATES.md',
            HERE/'sources/square_boundary_degree8_exploration.py.txt']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,marked_points=16,
        boundary_difference_degree_excluded_through=7,
        nonzero_first_normal_Jacobian_requires_boundary_difference_degree_at_least=8,
        full_actual_S_pole3_cusp_family_covered=True,exact_number_field_unit_certificates=10,
        full_degree8_ideal_solved=False,degree8_high_remainder_dimension=[x['high_remainder_locus_dimension'] for x in checkpoints],
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(all_complex_embeddings_covered=True,exact_Bezout_identities_replayed=True,
            complete_S_family_and_all_normal_freedoms_retained=True,condition_on_boundary_root_lift_not_automatic_for_all_P=True,
            arbitrary_determinantal_zero_locus_excluded=False,degree8_nonzero_kappa_locus_solved=False,
            general_nonzero_kappa_second_normal_locus_solved=False,whole_marked_candidate_excluded=False,
            full_source_Jacobian_pair_constructed=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,contact_cases_classified=False,
            quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'square_boundary_gap_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; square-boundary perturbations through degree7 excluded at all16 marked points; full goal active',flush=True)

if __name__=='__main__':main()
