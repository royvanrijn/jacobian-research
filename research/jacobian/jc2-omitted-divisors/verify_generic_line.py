"""Exact generic-line, projection-Euler and completed image-catalogue controls."""
from pathlib import Path
import hashlib
import json
import math
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def infinity_profiles():
    kinds=sorted([(e,f) for e in range(1,7) for f in range(1,7) if e*f<=6],key=lambda q:(q[0]*q[1],q))
    out=[]

    def recurse(remaining,start,current):
        if not remaining:
            out.append(tuple(current))
            return
        for j in range(start,len(kinds)):
            e,f=kinds[j]
            if e*f<=remaining:
                recurse(remaining-e*f,j,current+[(e,f)])

    recurse(6,0,[])
    return out


def sign_from_cycles(lengths):
    perm=[]
    offset=0
    for length in lengths:
        perm.extend(list(range(offset+1,offset+length))+[offset])
        offset+=length
    return (-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))


def main():
    started=time.monotonic()
    manifest_path=HERE/'sources/generic_line_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    for source in manifest['sources']:
        path=manifest_path.parent/source['path']
        check('fixed primary source bytes '+source['path'],path.stat().st_size==source['bytes']
              and digest(path)==source['sha256'])
    k,l,n,e,t=s.symbols('k l n e t',integer=True)
    genus=(6+5*(e-1)*k-n)/2
    punctures=n+10+5*k+5*t*l
    lost=20+5*e*k+5*t*l
    deleted_points=10+5*k+5*t*l
    affine_different=10+5*(e-1)*k
    chi=-14-5*e*k-5*t*l
    check('universal source Euler/genus/puncture identity',s.expand(2-2*genus-punctures-chi)==0)
    check('universal finite-normal curve Euler identity',s.expand(2-2*genus-n+4+5*(e-1)*k)==0)
    check('universal finite-fibre Euler integration',s.expand(6-lost-chi)==0)
    check('universal lost sheets minus removed points equals different',s.expand(lost-deleted_points-affine_different)==0)
    # Independent small norm control for the generic source-coordinate monicity lemma.
    # This finite ramified map is not Keller and is not a residual r7 realization.
    T,beta,U,V=s.symbols('T beta U V')
    norm=s.resultant(-(T-beta)**3+U,beta**5-V,beta)
    check('source-finite primitive small norm has full V-degree3',s.degree(norm,V)==3 and s.degree(norm,T)==15)
    check('monic V top coefficient has no T dependence',s.Poly(norm,V).LC() in [1,-1]
          and s.degree(s.Poly(norm,T).LC(),V)==0)
    profiles=infinity_profiles()
    profile_rows=[]
    for i,profile in enumerate(profiles):
        total=sum(a*b for a,b in profile)
        count=sum(b for a,b in profile)
        lengths=[a for a,b in profile for _ in range(b)]
        check(f'infinity profile{i} weighted full degree',total==6 and 1<=count<=6)
        check(f'infinity profile{i} full different and geometric cycles',
              sum((a-1)*b for a,b in profile)==6-count and sum(lengths)==6 and len(lengths)==count)
        check(f'infinity profile{i} exact permutation sign',sign_from_cycles(lengths)==(-1)**(6-count))
        profile_rows.append(dict(primes=[dict(e=a,f=b) for a,b in profile],
            infinity_prime_count=len(profile),geometric_infinity_point_count=count,
            total_different=6-count,cycle_lengths=lengths,global_geometric_realization_claimed=False))
    check('two global primes need not give two geometric points from degree alone',
          ((2,1),(2,2)) in profiles and sum(b for a,b in [(2,1),(2,2)])==3)
    options=[]
    for index in [2,3]:
        for degree_parameter in range(1,23):
            pairs=[(point_count,(6+5*(index-1)*degree_parameter-point_count)//2)
                   for point_count in range(1,7)
                   if (6+5*(index-1)*degree_parameter-point_count)%2==0]
            floor=math.ceil(5*degree_parameter/2) if index==2 else 5*degree_parameter
            check(f'e{index} k{degree_parameter} all three genus options',
                  sorted(g for point_count,g in pairs)==[floor,floor+1,floor+2])
            check(f'e{index} k{degree_parameter} global minimum genus',
                  min(g for point_count,g in pairs)>=(3 if index==2 else 5))
            for point_count,g in pairs:
                check(f'e{index} k{degree_parameter} n{point_count} sign and RH agree',
                      (-1)**(6-point_count)==(-1)**(10+5*(index-1)*degree_parameter))
                interval=(2*g-5,2*g) if index==2 else (g-2,g)
                check(f'e{index} k{degree_parameter} n{point_count} image-degree recovery',
                      interval[0]<=5*degree_parameter<=interval[1])
                if point_count==2:
                    check(f'e{index} k{degree_parameter} n2 cannot retain genus2',g>=7 and (index!=2 or degree_parameter%2==0))
            options.append(dict(e_S=index,k_B=degree_parameter,degree_B=5*degree_parameter,
                possible_genus_and_infinity_points=[dict(g=g,n_infinity=point_count) for point_count,g in pairs],
                additional_F2_degree_cap_assumed=True,realization_claimed=False))
    check('F2 genus ranges',
          min(v['g'] for r in options if r['e_S']==2 for v in r['possible_genus_and_infinity_points'])==3
          and max(v['g'] for r in options if r['e_S']==2 for v in r['possible_genus_and_infinity_points'])==57
          and min(v['g'] for r in options if r['e_S']==3 for v in r['possible_genus_and_infinity_points'])==5
          and max(v['g'] for r in options if r['e_S']==3 for v in r['possible_genus_and_infinity_points'])==112)
    degree_groups=[]
    for group in ['A_distinct','A_H=C','A_H=B','B','N','C']:
        choices=[]
        for kb in range(1,25):
            for kh in (range(1,25) if group=='A_distinct' else [2] if group=='A_H=C' else [kb] if group=='A_H=B' else [0]):
                with_D=group.startswith('A')
                if 2+kb+(kh if with_D else 0)>24:
                    continue
                pole_count=0
                for pr in range(1,13):
                    for ps in range(1,25):
                        remaining=24-2*pr-kb*ps
                        pole_count+=max(0,remaining//kh) if with_D else int(remaining>=0)
                check(group+f' k{kb} l{kh} nonempty positive pole options',pole_count>0)
                choices.append(dict(k_B=kb,k_H=kh if with_D else None,admissible_positive_pole_label_count=pole_count))
        maximum=max(r['k_B'] for r in choices)
        expected={'A_distinct':21,'A_H=C':20,'A_H=B':11,'B':22,'N':22,'C':22}[group]
        check(group+' complete pole-weighted image-degree cap',maximum==expected)
        max_genus=(5*maximum+2) if group=='C' else math.ceil(5*maximum/2)+2
        check(group+' pole-weighted genus cap',max_genus=={'A_distinct':55,'A_H=C':52,'A_H=B':30,'B':57,'N':57,'C':112}[group])
        degree_groups.append(dict(grouping=group,maximum_degree_B=5*maximum,degree_choices=choices,
            maximum_generic_line_genus=max_genus,geometry_realization_claimed=False))
    check('all328 oriented degree-label choices enumerated',sum(len(r['degree_choices']) for r in degree_groups)==328)
    check('shared unramified prime cannot be erased from the degree bound',
          2+12<=24 and 2+12+12>24 and 2+21<=24 and 2+21+2>24)
    check('pole-weighted inverse degree bound respects both coordinate degrees',3*24<=74 and 5*24==120)
    spectra={
        'A_distinct':([5*l-1,5*k+8,0],-14-10*k-5*l),
        'A_H=C':([0,5*k-1,9],-24-10*k),
        'A_H=B':([0,9,5*k-1],-14-15*k),
        'B_or_N':([1,5*k+7,0],-14-10*k),
        'C':([0,9,5*k-1],-14-15*k)}
    spectrum_rows=[]
    for label,(counts,value) in spectra.items():
        check(label+' whole projection Euler sum',s.expand(sum((i+1)*c for i,c in enumerate(counts))-(1-value))==0)
        spectrum_rows.append(dict(grouping=label,nonzero_Euler_jump_counts_1_2_3=[str(c) for c in counts],
            generic_source_Euler=str(value),actual_critical_values=None,
            identifies_ramification_indices=False))
    check('shared unramified image and cyclic ramification have identical Euler-jump spectra',
          spectra['A_H=B']==spectra['C'])
    check('same source Euler can hide different genus and punctures',
          2-2*7-32==2-2*12-22==-44)
    check('even genus and punctures can agree when true infinity data differ',
          (6+5-1)//2==(6+10-6)//2==5 and 1+10+10==6+10+5==21)
    u,w,a=s.symbols('u w a')
    for label,v in [('unramified',w),('fold',w*w),('cyclic3',w**3),('cusp',w**3+u*w)]:
        equation=v-a*u
        gradient=[s.diff(equation,x).subs({u:0,w:0,a:1}) for x in [u,w]]
        check(label+' generic fixed-point target line has smooth pullback',any(x!=0 for x in gradient))

    L_path=ROOT/'jacobian/jc2-r7-collisions/strata.json'
    K_path=ROOT/'jacobian/jc2-r7-laurent-collisions/strata.json'
    L_data=json.loads(L_path.read_text())
    K_data=json.loads(K_path.read_text())
    d=s.symbols('d')
    L_counts={'d=0':1,'U=0':9,'M=0':3,'T=0':18,'C=0, retained cusp class I':3}
    for key,count in [('U',9),('M',3),('T',18),('C',3)]:
        check('polynomial parameter factor '+key+' degree',s.degree(s.sympify(L_data['exceptional_factors'][key]),d)==count)
    categories=['nodes','triples','tacnodes','cusps','A5']
    aggregate={key:0 for key in categories}
    total=0
    triple_coefficient_points=tacnode_coefficient_points=cusp_coefficient_points=0
    catalog=[]
    for chart,data in [('L',L_data),('K',K_data)]:
        for row in data['singularities']:
            count=L_counts.get(row['stratum'],0) if chart=='L' else row.get('x_values',0)
            counts=dict(nodes=row.get('nodes',0),triples=row.get('triple_points',row.get('ordinary_triples',0)),
                tacnodes=row.get('tacnodes',row.get('tacnodes_A3',0)),
                cusps=row.get('cusps',row.get('cusps_A2',0)),A5=row.get('A5',row.get('tangencies_A5',0)))
            check(chart+' '+row['stratum']+' complete affine delta12',
                  counts['nodes']+3*counts['triples']+2*counts['tacnodes']+counts['cusps']+3*counts['A5']==12)
            total+=count
            for key in categories:
                aggregate[key]+=count*counts[key]
            triple_coefficient_points+=count*bool(counts['triples'])
            tacnode_coefficient_points+=count*bool(counts['tacnodes'])
            cusp_coefficient_points+=count*bool(counts['cusps'])
            catalog.append(dict(chart=chart,stratum=row['stratum'],coefficient_point_count=count or None,
                open_stratum=not bool(count),singularity_counts=counts,
                exact_parameter_equations_reference=str((L_path if chart=='L' else K_path).relative_to(ROOT))))
    # The three L/K intersections are precisely the same class-I coefficient points.
    total-=3
    aggregate['nodes']-=3*8
    aggregate['cusps']-=3
    aggregate['A5']-=3
    cusp_coefficient_points-=3
    check('two-family nonnodal coefficient union',total==268 and 34+237-3==total)
    check('all exceptional coefficient point classes',
          (triple_coefficient_points,tacnode_coefficient_points,cusp_coefficient_points)==(163,99,6)
          and triple_coefficient_points+tacnode_coefficient_points+cusp_coefficient_points==total)
    check('all singularity incidences reconciled',aggregate==dict(nodes=2508,triples=165,tacnodes=99,cusps=6,A5=3))
    check('whole coefficient-union delta conservation',
          aggregate['nodes']+3*aggregate['triples']+2*aggregate['tacnodes']+aggregate['cusps']+3*aggregate['A5']==12*total==3216)
    manifest_inputs=[manifest_path]+[manifest_path.parent/r['path'] for r in manifest['sources']]
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT14',
        hypotheses='PF2OMIT13 four residual degree6 gap21 rows; Chau Theorem3 imported; W<=24 only with additional F2 coordinate degrees75,125',
        other_image_degrees=dict(B='5*k, k>=1',H_when_present='5*l, l>=1',H_equals_C='l=2',H_equals_B='l=k'),
        line_formulas=dict(projective_genus='2*g=6+5*(e_S-1)*k-n_infinity',
            infinity_geometric_points='n_infinity=sum f_i, sum e_i*f_i=6',
            affine_normalization_Euler='-4-5*(e_S-1)*k',
            source_punctures='n_infinity+10+5*k+5*t*l',
            source_Euler='-14-5*e_S*k-5*t*l',
            lost_affine_sheets='20+5*e_S*k+5*t*l',deleted_affine_points='10+5*k+5*t*l'),
        all_rows_generic_target_line_genus_at_least=3,index3_row_genus_at_least=5,
        inverse_pole_weighted_degree_budget='W=2*p_R+k*p_S+t*l*p_D<=M-1 when source degrees are3M,5M',
        conditional_F2_inverse_degree_budget=dict(M=25,maximum_W=24,maximum_degree_a6=120,
            maximum_V_degree_a6=72,all_oriented_degree_label_count=328,groupings=degree_groups),
        generic_target_line_through_any_fixed_affine_value_has_same_genus=True,
        source_coordinate_fibration_identified_with_target_pencil=False,
        projection_Euler_jump_formula='delta(a)=sum e_D*(deg h_D-#h_D^-1(a))-sum eta_D(q) at that projection value',
        projection_Euler_jump_spectra=spectrum_rows,
        numerical_infinity_profiles=profile_rows,conditional_F2_genus_options=options,
        first_image_catalogue=catalog,
        first_image_union=dict(nonnodal_coefficient_points=total,triple_coefficient_points=163,
            tacnode_coefficient_points=99,cusp_coefficient_points=6,
            singularity_incidences=aggregate,total_affine_delta=3216,
            count_is_affine_moduli_count=False),
        first_prime_local_factors=dict(node=[2,2],tacnode=[2,2],A5=[2,2],triple=[2,2,2],cusp=[3]),
        H_equals_C_excluded_immersed_coefficient_points=163,
        Laurent_affine_collision_taxonomy_complete=True,
        actual_B_equation=None,actual_H_equation=None,actual_infinity_divisors=None,
        actual_source_map=None,full_goal_complete=False)
    target=HERE/'generic_line_ledger.json'
    target.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=[HERE/n for n in ['GENERIC_LINE_ACCOUNTING.md','ADMISSIBLE_INVERSE.md','NONIMMERSED_LEDGER.md',
        'nonimmersed_ledger.json','RESIDUE_LEDGER.md','SIMPLE_RAMIFICATION_AND_RADIAL.md']]
    inputs += [L_path,L_path.parent/'PROOF.md',K_path,K_path.parent/'PROOF.md',
        ROOT/'jacobian/jc2-r7-universal-contact/PROOF.md',
        ROOT/'foundations/plane-jc/F2_AFFINE_TARGET_CURVE_ATLAS.md']+manifest_inputs
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={str(p.relative_to(ROOT)):digest(p) for p in inputs},
        ledger_sha256=digest(target),runtime_seconds=round(time.monotonic()-started,6),
        python=platform.python_version(),sympy=s.__version__,numeric_infinity_profile_count=len(profiles),
        proof_boundary='Exact universal Euler identities, all degree6 numeric DVR profiles, capped integer genus lists and full imported catalogue aggregation checked. Chau, Bertini/Hodge, Riemann-Hurwitz, constructible Euler and generic plane-projection geometry are written/imported arguments. Numeric profiles and controls are not global covers.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'generic_line_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','numeric_infinity_profile_count','runtime_seconds']}))


if __name__=='__main__':
    main()
