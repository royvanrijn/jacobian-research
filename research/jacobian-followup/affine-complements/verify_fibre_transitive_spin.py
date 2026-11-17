"""PF104 exact controls; universal kernel/lattice/Euler proofs are written."""
from pathlib import Path
import hashlib,json,math,platform
import sympy as S
import full_peripheral_linking as N
import fibre_transitive_cabling as T
import multiple_boundary_contact as C
from verify_full_peripheral_linking import natural_seed

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    predecessors=['full_peripheral_linking_receipt.json','multiple_boundary_contact_receipt.json']
    for name in predecessors:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/fibre-transitive-spin/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/source pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('complete primary lattice paper read visually',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[1,2,3])
    r,s,e,a,F=S.symbols('r s e a F',nonzero=True)
    corrected=lambda R,G:R*R*G+a*(1-R*R)/e
    ident('UNBOUNDED full-shear two-stage telescoping',corrected(s,corrected(r,F))-corrected(r*s,F))
    ident('UNBOUNDED weighted balance retained',r*r*a+a*(1-r*r)-a)
    ident('UNBOUNDED pushed longitude correction',r*(a*r)+e*a*(1-r*r)/e-a)
    q=S.symbols('q',integer=True,positive=True)
    ident('UNBOUNDED spin Euler degree threshold',8*q+3-11-8*(q-1))
    degree,count,knowncost,knowncount,eta=S.symbols('d n C b eta')
    ident('UNBOUNDED all-prime budget upper bound',
          (degree-1-(count-knowncount+knowncost+eta)).subs(count,degree-1-knowncost+knowncount-eta))
    n,A=S.symbols('N A',integer=True,positive=True)
    ident('UNBOUNDED primitive kernel relation Nt-rA=1',n*((1+r*A)/n)-r*A-1)
    # Exact pattern matrices test the independent lifted-path implementation.
    pattern_controls=[]
    for pars in [(1,2,0,7,31),(3,2,0,7,43),(4,1,2,5,9),(4,1,2,9,13),
                 (5,2,1,11,21),(3,2,0,37,31)]:
        packet=T.pattern_packet(*pars)
        check('full cycle/order criterion '+str(pars),packet['fibre_transitive'] and packet['fibre_order']==pars[0]*pars[1])
        check('exact torsion-free pattern H1 and integral port basis '+str(pars),
              packet['H1']==dict(free=2,torsion=[]) and packet['quotient_by_proposed_basis']==dict(free=0,torsion=[]))
        check('both full-shear port identities INTEGRAL '+str(pars),packet['port_difference_orders']==[1,1])
        pattern_controls.append({k:v for k,v in packet.items() if k not in ['rows','ports']})
    failure=T.pattern_packet(3,1,1,4,7)
    check('without fibre transitivity torsion remains despite vanishing port differences',
          not failure['fibre_transitive'] and failure['fibre_order']==1 and
          failure['H1']==dict(free=2,torsion=[4,4]) and
          failure['quotient_by_proposed_basis']==dict(free=0,torsion=[4,4]) and
          failure['port_difference_orders']==[1,1])
    family=[]
    for initialpq,rho,lists,label in [
        ((2,3),natural_seed(),[[],[(7,31)],[(7,31),(13,79)],[(7,43)]*6,
             [(6*10**30+1,6*(6*10**30+1)*10**30+1)]],'natural8'),
        ((2,3),natural_seed(regular=True),[[],[(7,43)]],'regular24'),
        ((2,3),[(1,0,2,3),(0,2,3,1)],[[],[(5,9)],[(5,9),(9,13)]],'S4_nonzero_shear'),
        ((2,9),natural_seed(q=9),[[(7,31)]],'closed_Z3_retained')]:
        p0,q0=initialpq;before=N.peripheral_cover(p0,q0,[],rho);seed=N.rational_core_linking(before)
        mu=before['meridian_monodromy'];lam=before['longitude_monodromy']
        z=tuple(mu[lam[i]] for i in range(len(mu)));M=math.lcm(T.order(mu),T.order(lam))
        profiles=[(b['normal_index'],b['longitudinal_degree'],b['shear']) for b in before['boundaries']]
        check('seed satisfies full fibre-transitivity hypotheses '+label,all(T.fibre_transitive(*b) for b in profiles))
        for cables in lists:
            R=math.prod(r for r,k in cables)
            check('independent cable parameters preserve both permutation orders '+str((label,cables)),
                  all(r>1 and k>1 and math.gcd(r,k)==1 and r%M==k%M==1 for r,k in cables))
            for reverse in ([False,True] if len(cables)<=1 and all(r<1000 for r,k in cables) else [False]):
                after=N.peripheral_cover(p0,q0,cables,rho+[z]*len(cables),reverse)
                result=N.rational_core_linking(after)
                check('actual full peripheral permutations retained '+str((label,cables,reverse)),
                      after['meridian_monodromy']==mu and after['longitude_monodromy']==lam)
                check('full-shear ALL-LENGTH predicted matrix '+str((label,cables,reverse)),
                      result['matrix']==T.transfer_matrix(seed['matrix'],profiles,R))
                check('filled H1 and exact individual core orders retained '+str((label,cables,reverse)),
                      result['closed_H1']==seed['closed_H1'] and result['core_orders']==[
                          order//math.gcd(order,R) for order in seed['core_orders']])
                check('full profiles including shear retained '+str((label,cables,reverse)),
                      list(zip(result['normal_indices'],result['longitudinal_degrees'],result['shears']))==profiles)
                if not reverse:family.append(dict(label=label,cables=cables,winding_product=str(R),result=result))
        if label=='S4_nonzero_shear':
            check('naive r^2F fails at nonzero shear',T.transfer_matrix(seed['matrix'],profiles,5)==[['1/2']] and
                  str(25*S.Rational(seed['matrix'][0][0]))!='1/2')
    # E8 at rank8 is retained; no false impossibility at all ranks.
    B=2*S.eye(8)
    for i,j in [(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7)]:B[i,j]=B[j,i]=-1
    e8=C.canonical_packet(B,[1]*7+[3])
    check('E8 positive even unimodular countercontrol',B.det()==1 and e8['rank']==8 and e8['canonical_energy']=='8')
    check('formal rank8 d11 budget remains possible, NOT a realization',
          T.spin_budget(11,[(3,1)])['necessary_rank_interval_nonempty'])
    check('all d<=10 smooth all-odd homology-sphere source packets fail',
          all(not T.spin_budget(d)['necessary_rank_interval_nonempty'] for d in range(2,11)))
    eight=T.spin_budget(8,[(3,2)])
    check('eight-sheet actual cost6 leaves at most2 primes',eight['maximum_omitted_primes']==2 and not eight['necessary_rank_interval_nonempty'])
    check('regular24 cost24 exceeds d-1 even without spin',4*(3*2)>24-1)
    check('rank-one even determinant2 is NOT unimodular',S.Matrix([[2]]).det()!=1)
    check('rank-one odd unimodular is NOT even',S.Matrix([[1]])[0,0]%2==1)
    U=S.Matrix([[0,1],[1,0]])
    check('indefinite even unimodular rank2 shows definiteness is essential',U.det()==-1 and set(U.eigenvals())=={-1,1})
    # Genuine finite degree8 polynomial map; its complement is C* x C.
    x,y,t,p,q=S.symbols('x y t p q')
    P0=3*x*x+12*x*y**3+8*y**6;Q0=x*y+y**4
    ident('finite polynomial control Jacobian',S.diff(P0,x)*S.diff(Q0,y)-S.diff(P0,y)*S.diff(Q0,x)-6*x*x)
    ident('complete branch pullback multiplicities',P0**2-64*Q0**3-x**3*(9*x+8*y**3))
    ident('degree8 monic integrality identity',y**8-6*Q0*y**4+P0*y*y-3*Q0**2)
    check('linear-P irreducibility primitive coefficient gcd',S.gcd(y*y,y**8-6*q*y**4-3*q*q)==1)
    check('generic y equation monic degree8',S.Poly(y**8-6*q*y**4+p*y*y-3*q*q,y).degree()==8)
    G=S.groebner([P0,Q0],x,y)
    check('special fibre monomials prove local length8',
          list(G)==[3*x*x-4*y**6,x*y+y**4,y**7])
    ident('weighted homogeneous first target coordinate',3*x*S.diff(P0,x)+y*S.diff(P0,y)-6*P0)
    ident('weighted homogeneous second target coordinate',3*x*S.diff(Q0,x)+y*S.diff(Q0,y)-4*Q0)
    ident('index3 boundary has nonprimitive residue2',P0.subs(x,0)-8*(y*y)**3)
    ident('same normalized cusp residue parameter',Q0.subs(x,0)-(y*y)**2)
    proof=(P/'FIBRE_TRANSITIVE_SHEAR_RECURSION_AND_SPIN_BUDGET_EXCLUSION.md').read_text()
    normalized=' '.join(proof.lower().split())
    for word in ['COLORED','ANY extra unramified','realization is inferred','No algebraic','NOT the marked Keller source',
                 'not a surface','singular eight-sheet','No Rokhlin','NO torsion','a_i*(1-r^2)/e_i']:
        check('proof/scope guard '+word,word.lower() in normalized)
    files=['FIBRE_TRANSITIVE_SHEAR_RECURSION_AND_SPIN_BUDGET_EXCLUSION.md','fibre_transitive_cabling.py',
           'verify_fibre_transitive_spin.py','sources/fibre-transitive-spin/manifest.json']+predecessors
    files+= [x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED104',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        universal=dict(arbitrary_degree=True,arbitrary_finite_length=True,independent_coprime_parameters=True,
                       stable_nonzero_shear=True,F_rule='R^2 F0+diag((1-R^2)*a/e)',
                       filled_H1_preserved_integrally=True,core_class_multiplied_by_R=True,
                       finite_samples_used_as_proof=False),
        spin_source=dict(all_odd_indices_and_integral_homology_sphere_end=True,
                         smooth_actual_source_rank_positive_multiple8=True,field_degree_at_least11=True,
                         known_pairs_budget='r<=d-1-sum(ef-1)',natural_eight_budget=eight),
        actual_exclusions=dict(natural_eight_all_lengths_all_parameters_smooth_source=True,
                               arbitrary_extra_unramified_omitted_primes_allowed=True,
                               regular_twenty_four_all_lengths_normal_source=True,
                               specific_colored_covers_only=True),
        pattern_controls=pattern_controls,nontransitive_torsion_countercontrol={
             k:v for k,v in failure.items() if k not in ['rows','ports']},
        exact_cover_controls=family,E8_formal_countercontrol=e8,
        finite_degree8_polynomial_countercontrol=dict(genuine=True,Jacobian='6x^2',
             complement='C* times C, not A2',finite_Keller_realization=False,precise_seed_label_identification_claimed=False),
        assurance=dict(written_universal_proofs=True,primary_lattice_theorem_imported=True,
                       actual_higher_degree_genuine_iterated_colored_family_excluded=True,
                       all_representations_of_underlying_knots_excluded=False,
                       singular_natural_eight_sources_excluded=False,all_degree_eight_Keller_maps_excluded=False,
                       all_higher_degree_Keller_maps_excluded=False,global_Laurent_finite_algebra_closed=False,
                       topological_covers_declared_nonexistent=False,formal_verification=False,
                       external_review=False,external_review_deferred=True,research_goal_completed=False),
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'fibre_transitive_spin_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',claim='PF2RULED104',check_count=len(checks))))


if __name__=='__main__':main()
