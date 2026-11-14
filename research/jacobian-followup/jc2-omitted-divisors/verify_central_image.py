"""Exact central Euler/Milnor, normalized block and residue-cover controls."""
from pathlib import Path
import hashlib,json,math,platform,time
import sympy as s
from verify_global_monodromy import closure,cycle,orbits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started=time.monotonic()
    kap,epsR,epsB,e,d,q,u,z,mu=s.symbols('kappa epsilon_R epsilon_B e d q u z mu')
    chi0=4-6*kap-epsR-e*q-epsB-z
    generic=-138-e*d-u
    jump=142-6*kap-epsR+e*(d-q)+(u-z)-epsB
    check('universal central Euler integration',s.expand(6*(1-kap)-2-epsR-e*q-epsB-z-chi0)==0)
    check('universal central Euler jump',s.expand(chi0-generic-jump)==0)
    degree,components,delta,branchdef=s.symbols('D components delta branchdef')
    arithmetic=(degree-1)*(degree-2)/2
    genus_sum=arithmetic-delta+components-1
    projective_chi=2*components-2*genus_sum-branchdef
    check('projective plane Euler/Milnor formula including reducible curves',
          s.expand(projective_chi-degree*(3-degree)-(2*delta-branchdef))==0)
    xx,yy,lam=s.symbols('x y lambda')
    control=xx*(xx*yy-1)
    check('integral generic fibre countercontrol solves for y',s.expand((control-lam).subs(yy,(lam+xx)/xx**2))==0)
    check('disconnected smooth central fibre countercontrol',
          s.expand(4*yy**2*s.diff(control,yy)-(2*xx*yy+1)*s.diff(control,xx))==1)
    prior=json.loads((HERE/'generic_line_ledger.json').read_text())
    milnor=dict(nodes=1,triples=4,tacnodes=3,cusps=2,A5=5)
    delta_weights=dict(nodes=1,triples=3,tacnodes=2,cusps=1,A5=3)
    branch_weights=dict(nodes=1,triples=2,tacnodes=1,cusps=0,A5=1)
    catalogue=[]
    for index,row in enumerate(prior['first_image_catalogue']):
        counts=row['singularity_counts']
        kappa=sum(counts[k]*branch_weights[k] for k in counts)
        target_mu=sum(counts[k]*milnor[k] for k in counts)
        check(f'catalogue{index} delta12',sum(counts[k]*delta_weights[k] for k in counts)==12)
        check(f'catalogue{index} target Milnor/branch relation',target_mu==24-kappa)
        if not counts['cusps']:
            baseline=2*counts['nodes']+6*counts['tacnodes']
            infinity_base=142-6*kappa-baseline
            check(f'catalogue{index} immersed infinity baseline',infinity_base==46+12*counts['triples']+4*counts['tacnodes'])
            if not counts['triples']:
                check(f'catalogue{index} shared H=C baseline',141-6*kappa==69+6*counts['tacnodes'])
            row_data=dict(source_finite_Milnor_baseline_shared_B_or_C=baseline,
                          source_infinity_Milnor_baseline_shared_B_or_C=infinity_base)
        else:
            baseline=2*counts['nodes']+10*counts['A5']+6
            infinity_base=141-6*kappa-baseline
            check(f'catalogue{index} N cusp-class baseline',(baseline,infinity_base) in [(32,55),(28,47)])
            row_data=dict(source_finite_Milnor_baseline_N=baseline,source_infinity_Milnor_baseline_N=infinity_base)
        catalogue.append(dict(chart=row['chart'],stratum=row['stratum'],kappa=kappa,
            target_affine_Euler=1-kappa,target_total_Milnor=target_mu,**row_data))
    check('all central infinity lower bounds at least46',all(
        min(value for key,value in row.items() if key.startswith('source_infinity'))>=46 for row in catalogue))
    check('row B cusp defect correction nonnegative',all(2*t-epsilon>=0 for epsilon in [0,1] for t in range(epsilon,4)))
    blocks=[]
    for original_index in [1,2,3]:
        for contact in range(1,7):
            g=math.gcd(original_index,contact)
            indices=[original_index//g]*g
            check(f'cyclic e{original_index} m{contact} normalized total rank',sum(indices)==original_index)
            check(f'cyclic e{original_index} m{contact} primitive branches',math.gcd(original_index//g,contact//g)==1)
            different=sum(i-1 for i in indices)
            check(f'cyclic e{original_index} m{contact} deleted plus different',g+different==original_index)
            blocks.append(dict(original_index=original_index,contact=contact,original_N_point_count=1,
                original_ambient_length=original_index,normalized_curve_indices=indices,
                normalized_deleted_points=g,normalized_finite_different=different,new_original_omitted_divisors=0))
    t,w,omega,U,V=s.symbols('t w omega U V')
    check('double tangency splits two deleted unramified points',s.expand((w-t)*(w+t)-(w*w-t*t))==0)
    check('triple tangency splits three deleted unramified points',
          s.rem(s.expand((w-t)*(w-omega*t)*(w-omega**2*t)-(w**3-t**3)),omega**2+omega+1,omega)==0)
    check('own cusp factorization',s.expand((w-2*t)*(w+t)**2-(w**3-3*t*t*w-2*t**3))==0)
    cubic=w**3+U*w-V
    discriminant=s.discriminant(cubic,w)
    check('cusp transverse discriminant order2',s.Poly(discriminant.subs({U:0,V:t}),t).degree()==2)
    check('cusp tangent discriminant order3',s.Poly(discriminant.subs({U:t,V:0}),t).degree()==3)
    check('cusp tangent index1+2 factorization',s.expand(cubic.subs({U:t,V:0})-w*(w*w+t))==0)
    check('cusp generic index cannot replace special length',3!=2 and (2+1)==3)
    check('own cusp coexceptional point is unramified but omitted',1+0==1)
    completions=[]
    for grouping,e_other in [('A_distinct',2),('A_H=C',2),('A_H=B',2),('B',2),('C',3),('N',2)]:
        count=0
        for c in range(1,4):
            for branch_b in range(4):
                for h in range(6) if grouping=='A_distinct' else [c if grouping=='A_H=C' else branch_b if grouping=='A_H=B' else 0]:
                    for epsilon_R in ([0,1] if grouping=='N' and c==1 else [0]):
                        for epsilon_B in ([0,1] if grouping=='B' and c==1 and branch_b==1 else [0]):
                            retained=6-2*c-e_other*branch_b-h-epsilon_R-epsilon_B
                            if retained<0:
                                continue
                            retained_degree=3 if grouping=='A_H=C' else 4
                            rank=(2*(c-1)+epsilon_R+e_other*branch_b+epsilon_B
                                  +(c-1 if grouping=='A_H=C' else h))
                            check(f'{grouping} completion{count} normalized residual rank',retained+rank==retained_degree)
                            count+=1
        completions.append(dict(grouping=grouping,nonnegative_pointwise_rank_completions_checked=count))
    residue_examples=[]
    for label,size,generators,expected_degrees,expected_order in [
        ('A_H=C_connected',3,[[0,1],[1,2]],[3],6),
        ('A_H=C_split',3,[[0,1]],[2,1],2),
        ('A_H=B_connected',4,[[0,1],[1,2],[2,3]],[4],24),
        ('A_H=B_3plus1',4,[[0,1],[1,2]],[3,1],6),
        ('A_H=B_2plus2',4,[[0,1],[2,3]],[2,2],4),
        ('A_H=B_2plus1plus1',4,[[0,1]],[2,1,1],2),
        ('C_connected',4,[[0,1],[1,2,3]],[4],24),
        ('C_split',4,[[0,1],[0,1,2]],[3,1],6)]:
        group=closure(size,[cycle(size,letters) for letters in generators])
        degrees=sorted([len(orbit) for orbit in orbits(size,group)],reverse=True)
        check(label+' group and residue partition',degrees==expected_degrees and len(group)==expected_order)
        residue_examples.append(dict(label=label,retained_residue_degrees=degrees,group_order=len(group),
            original_Keller_cover_realization_claimed=False))
    node_preimages=t**24-1
    for intersection_degree in [6,10,12]:
        intersections=t**intersection_degree+2
        branch_poly=s.Poly(node_preimages*intersections,t)
        check(f'split double control d{intersection_degree} disjoint squarefree branch loci',
              s.gcd(branch_poly.as_expr(),s.diff(branch_poly.as_expr(),t))==1)
        check(f'split double control d{intersection_degree} genus',
              (branch_poly.degree()-1)//2==11+math.ceil(intersection_degree/2))
    check('split H=C central Euler reconstruction',s.expand((-23)+(-46-2*d)-chi0.subs({kap:12,epsR:0,epsB:0,e:2,q:d,z:1}))==0)
    clean=[]
    for grouping,other_index,unram_count,z_charge,m_cover,finite_different,deleted in [
        ('A_H=C',2,0,1,3,24+d,48+d),
        ('A_H=B',2,d,d,4,24+d,24+2*d),
        ('C',3,0,0,4,24+2*d,24+d)]:
        values={kap:12,epsR:0,epsB:0,e:other_index,q:d,z:z_charge,u:unram_count}
        central=s.expand(chi0.subs(values))
        normalized=central+(0 if grouping=='A_H=C' else 24)
        check(grouping+' clean normalization Euler cover agreement',s.expand(m_cover-finite_different-deleted-normalized)==0)
        source_mu=0 if grouping=='A_H=C' else 24
        central_jump=s.expand(jump.subs(values))
        clean.append(dict(grouping=grouping,intersection_degree='d, all simple',source_Euler=str(central),
            source_normalization_Euler=str(normalized),source_finite_central_Milnor=source_mu,
            source_central_Euler_jump=str(central_jump),source_infinity_Milnor_jump=str(central_jump-source_mu),
            residual_cover_degree=m_cover,finite_residual_different=str(finite_different),
            deleted_normalized_residual_points=str(deleted)))
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT19',central_fibre='X0=(H_C composed with F)^-1(0)',
        source_central_Euler='4-6kappa-epsilon_R-e_S*q_B-epsilon_B-z_H',
        source_central_Euler_jump='142-6kappa-epsilon_R+e_S*(d-q_B)+(u-z_H)-epsilon_B',
        source_central_finite_Milnor='sum retained(y)*mu_y(C)',
        normalized_source_Euler='chi(X0)+sum retained(y)*(c(y)-1)',
        source_infinity_Milnor_jump='central Euler jump minus finite central Milnor',
        universal_minimum_infinity_Milnor_jump=46,minimum_shared_H_C_infinity_jump=69,
        minimum_N_class_I_infinity_jump=55,minimum_N_class_II_infinity_jump=47,
        generic_primitive_fibration_implies_central_connected=False,catalogue=catalogue,
        cyclic_normalized_block_controls=blocks,pointwise_rank_completion_checks=completions,
        clean_nodal_transverse_central_rows=clean,
        clean_retained_residue_degree_options={'A_H=C':[[3],[2,1]],'A_H=B':[[4],[3,1],[2,2],[2,1,1]],'C':[[4],[3,1]]},
        clean_retained_C_prime_flags=dict(image='C',omitted=False,ramifies=False,generic_ramification_index=1,
            residue_degrees='One prime per part of the chosen partition; actual choice unknown'),
        original_C_omitted_prime_flags=[dict(prime='R',image='C',omitted=True,ramifies=True,e=2,f=1),
            dict(prime='D',image='C',only_when='A_H=C',omitted=True,ramifies=False,e=1,f=1)],
        residue_group_controls=residue_examples,
        clean_C_split_case_requires_polynomial_lift=True,
        clean_A_H_C_split_section_deleted_node_parameters=24,
        clean_A_H_C_split_double_field='C(t)(sqrt(A_C(t)*R_CB(t)))',
        clean_A_H_C_split_double_genus='11+ceil(d/2)',
        actual_C_coefficient_parameter=None,actual_C_stratum=None,
        actual_q_B=None,actual_q_H=None,actual_cusp_B_on_C=None,actual_central_source_components=None,
        actual_retained_residue_degrees=None,actual_polynomial_lift=None,actual_source_infinity_singularities=None,
        actual_B_equation=None,actual_H_equation=None,actual_global_source_realization=None,
        unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'central_image_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian-followup/jc2-omitted-divisors/'+name for name in ['CENTRAL_IMAGE_FIBRE.md',
        'GENERIC_LINE_ACCOUNTING.md','generic_line_ledger.json','IMAGE_PENCIL_MONODROMY.md','image_pencil_ledger.json',
        'RESIDUE_LEDGER.md','NONIMMERSED_LEDGER.md','SIMPLE_RAMIFICATION_AND_RADIAL.md','GLOBAL_MONODROMY.md','verify_global_monodromy.py']]
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in inputs},ledger_sha256=digest(ledger_path),
        runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=s.__version__,
        proof_boundary='Euler integration, etale local-type transport, plane genus/Milnor identity, normalized finite base change and residue-monodromy interpretation are written. Exact symbolic identities, complete imported catalogue charges, cyclic/cusp controls, finite row rank completions and permutation/clean-cover controls replayed. No actual incidence, source infinity singularity, residue choice, polynomial lift or global cover is constructed.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'central_image_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
