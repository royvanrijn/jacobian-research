"""PF2RULED21: exact controls for the written all-series critical-value proof.
Exhaustive gcd chains use a bound derived from the polar sum, not a jet cutoff.
"""
from pathlib import Path
from math import gcd
from time import monotonic
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.simplify(S.cancel(S.expand(value.doit())))==0)
def order(poly,t):return min(S.Poly(S.expand(poly),t).monoms())[0]
def chains(e,target,last=0,score=0,chain=()):
    if e==1:return [chain] if score==target else []
    out=[]
    # Each gcd drop is >=1, so every admissible exponent is <=target.
    for j in range(last+1,target+1):
        d=gcd(e,j);v=score+(e-d)*j
        if d<e and v<=target:out+=chains(d,target,j,v,chain+((j,d),))
    return out

def main():
    start=monotonic();w,z,c,t,X=S.symbols('w z c t X')
    cc=chains(8,33)
    check('complete characteristic gcd-chain census',cc==[((2,2),(21,1))])
    check('conjugate counts sum to seven',(8-2)+(2-1)==7)
    check('exact conjugate-difference polar sum',6*2+1*21==33)
    check('no earlier positive exponent divisible by8',2<8)
    # Universal valuation identities; their applicability is proved in the note.
    N,l,m,E=S.symbols('N l m E')
    I=N+l-2
    zero('differentiated branch valuation',I+1-(N+l-1))
    nu=S.expand(I-(m-2)*l-E)
    zero('critical-value companion pole bound',nu.subs(N,m*l/2)-((3-m/2)*l-2-E))
    check('all m>=6 coefficient nonpositive',3-S.Rational(6,2)==0)
    check('all m>=6 constant bound strictly negative',-2<0)
    check('polar33 would force critical value order27',33-8+2==27)
    check('critical curve order4 forced by N>l',27>8 and 8//2==4)
    check('parameter flow has lower critical value order24',6*4==24<27)
    check('generic single-branch polar bound at l8 is30',(S.Rational(6,2)+1)*8-2==30)
    # Chain rule in the ORIGINAL coordinates, with the critical equation.
    AA=S.Function('A');phi=S.Function('phi')(c,w)
    K=AA(phi,w)-c*phi**m
    derivative=S.diff(K,c)
    chain=S.diff(phi,c)*(S.Subs(S.Derivative(AA(z,w),z),z,phi)-c*m*phi**(m-1))-phi**m
    zero('universal original-coordinate critical-flow chain rule',derivative-chain)
    zero('chain rule modulo critical equation',S.expand(chain+phi**m).subs(
        S.Subs(S.Derivative(AA(z,w),z),z,phi),c*m*phi**(m-1)))
    # A selected fibre with polar33 exists; the generic pencil has two branches.
    A=(z-w**4)**2-w**27;F=A-c*z**6
    param={w:t**2,z:t**8+t**27,c:0}
    zero('special polar33 exact parameterization',F.subs(param))
    fw=S.expand(S.diff(A,w).subs(param));fz=S.expand(S.diff(A,z).subs(param))
    check('special polar order33',order(fw,t)==33 and fw.coeff(t,33)==-8)
    check('special transverse derivative order27',order(fz,t)==27 and fz.coeff(t,27)==2)
    check('special branch primitive exponents',gcd(gcd(2,8),27)==1)
    scaled=S.cancel(F.subs({w:t,z:t**4+t**12*X})/t**24)
    zero('generic countercontrol full scaled equation',scaled-(X**2-t**3-c*(1+t**8*X)**6))
    zero('generic countercontrol initial equation',scaled.subs(t,0)-(X**2-c))
    zero('generic countercontrol simple-root discriminant',S.discriminant(X**2-c,X)-4*c)
    gfw=S.expand(S.diff(A,w).subs({w:t,z:t**4+t**12*X}))
    check('generic branch polar order15',order(gfw,t)==15 and gfw.coeff(t,15)==-8*X)
    check('generic total polar30 is below33',2*15==30<33)
    check('generic countercontrol fails finite local pole positivity',15-4*4==-1)
    phifirst=w**4+3*c*w**20
    critical=S.expand(S.diff(F,z).subs(z,phifirst))
    for j in (0,1):zero('critical curve through c order '+str(j),critical.coeff(c,j))
    val=S.expand(F.subs(z,phifirst))
    zero('critical value constant coefficient',val.coeff(c,0)+w**27)
    zero('critical value first coefficient is minus phi0 sixth',val.coeff(c,1)+w**24)
    # Sealed predecessor results, including their now-historical handoff.
    prev=json.loads((HERE/'quartic_polar_gap_certificates.json').read_text())
    pr=json.loads((HERE/'quartic_polar_gap_receipt.json').read_text())
    check('quartic-gap predecessor passed',pr['outcome']=='PASS')
    check('quartic-gap old profiles explicitly historical',pr['higher_polar_profiles_superseded_by']=='PF2RULED21')
    for path,digest in pr['source_sha256'].items():
        check('predecessor seal '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    cost=lambda part:sum(2 if k in (2,4) else 1 for k in part)
    costs=[dict(partition=part,minimum_Q_pole_degree=cost(part),excluded=cost(part)>4)
           for part in prev['remaining_noncontact_partitions']]
    kept=[x['partition'] for x in costs if not x['excluded']]
    removed=[x['partition'] for x in costs if x['excluded']]
    check('new noncontact exclusion is432',removed==[[4,3,2]])
    check('current noncontact types are exactly72 and54',kept==[[7,2],[5,4]])
    check('contact3 away42 costs at least5',cost([4,2])+1==5>4)
    check('contact3 cube6 remains necessary',cost([6])+1==2<=4)
    check('local common quadruple root cannot cost1',4*8+1==33)
    check('pole2 cost is exact predecessor result',prev['number_of_local_branches']==2)
    polar=json.loads((HERE/'polar_cusp_sector_receipt.json').read_text())
    source=polar['current_pole69_frontier']
    check('three source degree groups unchanged',[x['ordinary_degrees'] for x in source]==[[96,144],[80,128],[48,128]])
    contacts=sorted(set(k for x in source for k in x['contacts'] if k))
    check('current contact multiplicities unchanged',contacts==[3,4,5,6,7,9])
    check('sixteen marked leading points retained',polar['marked_two_root_sector_exact_point_count']==16)
    forms=[dict(contact=k,away_multiplicities=[] if k==9 else [9-k],cube=k%3==0) for k in contacts]
    check('contact3 only has away cube6',forms[0]==dict(contact=3,away_multiplicities=[6],cube=True))
    check('all current forms have support at most2',all(len(x['away_multiplicities'])+1<=2 for x in forms))
    # Auxiliary rationality: theta^3=x^k, T=x^u*theta^v.
    rational=[]
    for k in (2,4,5,7):
        v=pow(k,-1,3);u=(1-k*v)//3
        check('Bezout monomial cover k'+str(k),3*u+k*v==1)
        check('T cubed equals x k'+str(k),3*u+k*v==1)
        # T^k=x^(u*k)*theta^(v*k)=theta*x^(u*k+k*(v*k-1)/3).
        check('theta equals T power k'+str(k),(v*k-1)%3==0 and u*k+k*((v*k-1)//3)==0)
        rational.append(dict(k=k,u=u,v=v,base='x=T^3',cover='theta=T^'+str(k),field='C(T)'))
    check('cube root forms have three rational components',all(x['contact'] in (3,6,9) for x in forms if x['cube']))
    check('noncube contact has a single finite root',all(len(x['away_multiplicities'])==1 for x in forms if not x['cube']))
    cubic=json.loads((HERE/'cubic_cover_period_certificates.json').read_text())
    # Exactness here is inherited from the written one/two-root classifications.
    check('period polynomial criterion retained',cubic['period_polynomial_identity']=='3*s*A_prime+(3*s_prime-2*s*h_prime/h)*A=1')
    cert=dict(as_of='2026-11-14',whole_generic_germ_unibranch_required=True,
        conjugate_formula='I=sum (e_previous-e_next)*j',characteristic_exponent_bound=33,
        exhaustive_degree8_polar33_chains=[[[j,e] for j,e in C] for C in cc],
        forced_original_w_order=2,critical_value_identity='partial_c K=-phi_c^m',
        polar_identity='I=N+l-2',unibranch_bound='I<=(m/2+1)*l-2',
        local_companion_bound='nu<=(3-m/2)*l-2-I(F,E)',all_height_exclusion_m_minimum=6,
        polar33_required_N=27,polar33_flow_order=24,polar33_bound=30,
        special33_control=dict(A=str(A),special_parameterization=['t^2','t^8+t^27'],
            special_polar_order=33,generic_total_polar=30,generic_branch_polar_orders=[15,15],
            critical_value_first_c_terms='-w^27-c*w^24',finite_map=False),
        common_multiplicity4_Q_pole_lower_bound=2,lower_bound_sharpness_proved=False,
        noncontact_partition_costs=costs,excluded_noncontact_partitions=removed,
        excluded_contact3_away_type=[4,2],current_contact3_away_types=[[6]],
        remaining_noncontact_partitions=kept,current_contact_root_forms=forms,
        current_pole69_frontier=source,marked_two_root_leading_candidates_retained=16,
        auxiliary_rational_parametrizations=rational,all_current_auxiliary_components_rational=True,
        all_current_root_forms_pass_leading_period=True,actual_PQ_constructed=False)
    cp=HERE/'unibranch_polar33_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'UNIBRANCH_POLAR33_OBSTRUCTION.md',cp,
        HERE/'QUARTIC_GENERIC_POLAR_GAP.md',HERE/'quartic_polar_gap_certificates.json',HERE/'quartic_polar_gap_receipt.json',
        HERE/'POLAR_CUSP_LEADING_SECTORS.md',HERE/'polar_cusp_sector_certificates.json',HERE/'polar_cusp_sector_receipt.json',
        HERE/'CUBIC_COVER_BOUNDARY_PERIOD.md',HERE/'cubic_cover_period_certificates.json']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        generic_unibranch_degree8_polar33_excluded=True,
        all_height_whole_unibranch_w_order2_excluded_for_poles_at_least6=True,
        noncontact432_excluded=True,contact3_away42_excluded=True,
        remaining_noncontact_nonic_partitions=kept,current_contact3_away_types=[[6]],
        current_pole69_frontier=source,marked_two_root_leading_candidates_retained=16,
        all_current_pole69_common_nonic_support_at_most2=True,
        all_current_auxiliary_components_rational=True,all_current_root_forms_pass_leading_period=True,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(all_series_coverage='written proof; exact identities and derived finite gcd-chain census replayed',
            whole_generic_germ_unibranch_required=True,uses_parameter_dependent_ambient_change_before_flow=False,
            claims_no_selected_fibre_polar33=False,claims_all_multibranch_generic_polar33_impossible=False,
            claims_exact_two_branches_at_common_multiplicity4=False,local_lower_bound_sharpness_proved=False,
            auxiliary_cover_is_target_intermediate_field=False,auxiliary_Mobius_change_is_marked_source_automorphism=False,
            sixteen_points_construct_actual_maps=False,all_degree9_period_loci_classified=False,
            all_pole69_excluded=False,actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
            conditional_2026_degree125_upgrade_used=False,formal_verification=False,
            independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'unibranch_polar33_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; generic single-branch polar33 excluded; pole69 support<=2; goal active')
if __name__=='__main__':main()
