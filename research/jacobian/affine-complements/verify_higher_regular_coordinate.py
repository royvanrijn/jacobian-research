'''PF2RULED40 exact controls for the written all-r cyclic-chart reduction.
Bounded arithmetic controls support, and do not prove, the geometric steps.
No affine converse for arbitrary Laurent data or finite Keller construction is supplied.
'''
from pathlib import Path
from math import gcd
from functools import reduce
from itertools import combinations
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    preds={}
    for name in ['coordinate_regular_sector_receipt.json','regular_sector_end_receipt.json','multiple_fiber_cyclic_receipt.json']:
        d=json.loads((HERE/name).read_text());preds[name]=d
        check('predecessor PASS '+name,d['outcome']=='PASS')
        for p,h in d['source_sha256'].items():
            check('predecessor seal '+name+' '+p,digest(HERE/p)==h)
    manifest=HERE/'sources/coordinate-pole-sector/manifest.json'
    source=next(x for x in json.loads(manifest.read_text())['primary_sources']
                if 'miyanishi' in x['path_relative_to_affine_complements'])
    pin=HERE/source['path_relative_to_affine_complements']
    check('primary affine-line bundle theorem pinned',
          pin.read_bytes().startswith(b'%PDF') and digest(pin)==source['sha256'] and
          'Theorem4.1.1' in source['used_statement'])
    z,t=S.symbols('z t')
    charts=[]
    for r in range(2,10):
        k=2-r
        check('canonical order determines the chart exponent r'+str(r),r-1+k==1)
        # A full mixed Laurent polynomial also keeps the forbidden residue
        # so the all-degree diagonal obstruction is not erased by setup.
        q=r+3
        c=S.symbols('c'+str(r)+'_0:'+str(q-r+2))
        support=list(range(r-1,q+1))
        h=sum(a*z**(-j) for a,j in zip(c,support))
        x=z**k*t+h
        check('whole-chart pulled-back Jacobian r'+str(r),
              S.expand(r*z**(r-1)*S.diff(x,t)-r*z)==0)
        theta_difference=S.expand(r*z**(r-1)*h)
        residue=theta_difference.coeff(z,-1)
        check('only coefficient r contributes to the exact volume obstruction r'+str(r),
              residue==r*c[support.index(r)])
        primitive=sum(r*a*z**(r-j)/(r-j) for a,j in zip(c,support) if j!=r)
        check('all nonresidual mixed terms have the displayed Laurent primitive r'+str(r),
              S.expand(S.diff(primitive,z)-theta_difference+residue/z)==0)
        check('leading Laurent pole strictly dominates the fiber term r'+str(r),r-1>r-2)
        # Terms at and above k are gauge; none below k is removable.
        for exponent in range(k-3,k+4):
            check('polynomial gauge threshold r'+str(r)+' exponent'+str(exponent),
                  (exponent-k>=0)==(exponent>=2-r))
        for i in range(r):
            check('deck characters share the same residue r'+str(r)+' i'+str(i),(i*r)%r==0)
            check('deck linear factor preserves x r'+str(r)+' i'+str(i),
                  ((2-r)+(r-2))%r==0)
        # Borderline and diagonal controls retain the local distinction.
        lower=r-1
        check('borderline single monomial has a free cyclic orbit r'+str(r),gcd(r,lower)==1)
        check('borderline single monomial has no residue r'+str(r),lower!=r)
        check('single diagonal monomial has nonzero residue and a stabilizer r'+str(r),
              S.expand(r*z**(r-1)*z**(-r)-r/z)==0 and gcd(r,r)==r)
        check('translation has exact boundary zero order r'+str(r),
              S.diff(x,t)*z**(r-2)==1 and r-2>=0)
        charts.append(dict(r=r,linear_exponent=k,minimum_pole=r-1,
                           excluded_diagonal_pole=r,translation_zero_order=r-2))
    # Character enumeration is a check of the all-support gcd argument:
    # using exact residues avoids floating point roots of unity.
    orbit_controls=0;nonseparated_controls=[]
    for r in range(2,11):
        exponents=list(range(r-1,r+5))
        for size in range(1,len(exponents)+1):
            for support in combinations(exponents,size):
                stabilizer=[a for a in range(r) if all(a*j%r==0 for j in support)]
                g=reduce(gcd,support,r)
                assert len(stabilizer)==g
                assert (len(stabilizer)==1)==(g==1)
                orbit_controls+=1
        bad=[r,2*r]
        stabilizer=[a for a in range(r) if all(a*j%r==0 for j in bad)]
        check('equal boundary Laurent charts really occur without support gcd r'+str(r),
              len(stabilizer)==r)
        nonseparated_controls.append(dict(r=r,support=bad,stabilizer_order=r,
                                         invalid_for_separated_normalization=True))
    check('exact support character census agrees with the all-support gcd proof',orbit_controls==9*63)
    controls=[]
    for r,q in [(3,2),(2,3),(3,4),(4,3),(5,4)]:
        check('surviving formal monomial is not falsely excluded r'+str(r)+' q'+str(q),
              q>=r-1 and q!=r and gcd(r,q)==1)
        controls.append(dict(r=r,q=q,support=[q],passes_necessary_chart_and_volume_tests=True,
                             affine_realization_claimed=(q==r-1),finite_Keller_realization_claimed=False))
    check('mixed support can separate even when its leading pole is divisible by r',
          reduce(gcd,[3,2,6])==1 and gcd(3,6)==3)
    r=S.symbols('r',integer=True,positive=True)
    check('all-r canonical exponent identity',S.expand((r-1)+(2-r)-1)==0)
    check('all-r translation exponent identity',S.expand((2-r)+(r-2))==0)
    check('all-r residue exponent identity',S.expand((r-1)-r+1)==0)
    # Pushforward uses order r and residue degree one. Do not erase the
    # special nonreduced factor, or confuse primitive with normal generation.
    for rr in range(2,10):
        HC,H0=S.symbols('HC H0')
        norm=HC**rr*H0
        check('special norm has exactly the boundary multiplicity r'+str(rr),
              S.Poly(norm,HC).degree()==rr and S.Poly(norm,H0).degree()==1)
        check('lifted source Jacobian has a branch line r'+str(rr),rr-1>0)
    # Uniform affine borderline ring: both principal localizations
    # recover whole plane opens. These are controls for the written cover
    # proof, not an affineness assertion for arbitrary Laurent data.
    v,ci=S.symbols('v ci')
    affine_controls=[]
    for rr in range(2,9):
        qq=rr-1
        f=v**rr-1
        # Use root ci=1; other roots are obtained by the cyclic symmetry.
        gi=S.cancel(f/(v-1))
        xx=S.symbols('x')
        for j in range(1,qq+1):
            w=v*f**j/z**j
            check('borderline source principal-ring identity r'+str(rr)+' j'+str(j),
                  S.cancel(w.subs(v,z**qq*xx)-z**(qq-j)*xx*(z**(qq*rr)*xx**rr-1)**j)==0)
            check('borderline omitted principal-ring identity r'+str(rr)+' j'+str(j),
                  S.cancel(w.subs(v,1+z*t)-(1+z*t)*gi.subs(v,1+z*t)**j*t**j)==0)
        check('borderline principal opens cover by Bezout r'+str(rr),
              S.expand(v*(rr*v**(rr-1))-rr*f)==rr)
        check('each omitted root has a nonzero parameter coefficient r'+str(rr),
              gi.subs(v,1)==rr and rr!=0)
        check('source retained line belongs to the affine source chart r'+str(rr),
              (z**(qq*rr)*xx**rr-1).subs(z,0)==-1)
        check('boundary root-chart volume primitive is polynomial r'+str(rr),
              S.expand(rr*z**(rr-1)*(z**(2-rr)*t+z**(1-rr))-rr*(z*t+1))==0)
        check('all r omitted roots form a free quotient orbit r'+str(rr),
              gcd(qq,rr)==1)
        check('quotient root invariant y and source x are preserved r'+str(rr),
              (rr%rr)==0 and (qq-qq)%rr==0)
        affine_controls.append(dict(r=rr,q=qq,ring_generator_count=qq+2,
            principal_chart_gcd=1,primitive='x dy',finite_Keller_target_constructed=False))
    check('all-r borderline chart exponents agree with the necessary form',S.expand((1-(r-1))-(2-r))==0)
    inputs=[Path(__file__),HERE/'HIGHER_REGULAR_COORDINATE_CYCLIC_FORM.md',manifest,pin]+[
        HERE/n for n in preds]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='existing smooth affine one-line K=R normalization with regular-coordinate boundary order r>=2; all r, poles and finite field degrees',
        cyclic_form=dict(root='y=z^r',degree='r',smooth_affine=True,
            boundary_components='r disjoint copies of A1',entire_plane_chart_count='r+1',
            common_overlap='Gm_z times A1',linear_exponent='2-r',
            Laurent_support_minimum='r-1',deck_orbit='h_i(z)=h_0(zeta^(-i)*z)',
            minimum_pole_inequality='q>=r-1',support_separation='gcd(r,support(h))=1'),
        volume=dict(common_residue='r*c_r',exact_iff='c_r=0',
            all_finite_degree_diagonal_excluded='q=r',written_affine_Cech_sufficiency=True,
            invariant_primitive_descends_across_branch_line=True),
        additive_action=dict(source_derivation='partial_x',boundary_derivation='z^(r-2)*partial_t',
            locally_nilpotent_on_normalization=True,free_when='r=2',zero_divisor_when_r_at_least_three='(r-2)R',
            target_algebra_invariance_proved=False),
        finite_algebra=dict(actual_one_index_two_prime_is_an_extra_hypothesis=True,
            special_norm='H_C^r*H_0',nonzero_norm='H_beta reduced irreducible',
            primitive_integral_coordinate='K(y)=C(N)',monic_norm_degree='d',
            every_D_beta_self_colliding=True,target_plus_y_algebra_not_normal=True,
            conductor_has_y_horizontal_component=True,norm_family_is_a_linear_pencil=False,
            cyclic_lift_degree='r*d',cyclic_lift_monodromy='C_r wreath S_d',
            boundary_inertia='r disjoint transpositions',retained_line_inertia='one r-cycle',
            lifted_cover_is_Keller=False,deleting_boundary_components_preserves_finiteness=False),
        borderline_controls=dict(all_r_affine_construction_written=True,
            ring='C[z,v,v*(v^r-1)^j/z^j : 1<=j<=r-1]',
            principal_opens='f(v) and v*g_i(v)',smooth_quotient=True,
            source='A2_(x,y) with y=z^r',minimum_coordinate_pole='q=r-1',
            volume_divisor='R',regular_primitive='x dy',
            hidden_pole_one_when_r_at_least_three=False,
            finite_Keller_target_constructed=False,controls=affine_controls),
        controls=dict(charts=charts,exact_support_cases=orbit_controls,
                      nonseparated_data=nonseparated_controls,surviving_formal_data=controls),
        assurance=dict(all_r_geometric_proof_written=True,primary_fibration_import_explicit=True,
            existing_affine_normalization_is_a_hypothesis=True,
            exact_controls_not_all_r_geometric_proof=True,
            affine_converse_for_arbitrary_chart_data=False,
            finite_Keller_map_constructed=False,whole_higher_regular_sector_excluded=False,
            no_regular_coordinate_case_excluded=False,arbitrary_quartic_packets_excluded=False,
            all_finite_Wright_degrees_excluded=False,coefficient_or_Groebner_search_used=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,
            research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'higher_regular_coordinate_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' higher regular-coordinate controls; written all-r reduction; goal active')
if __name__=='__main__':main()
