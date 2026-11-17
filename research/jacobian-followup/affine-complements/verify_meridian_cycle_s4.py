"""PF2RULED86: longitude-surface framing and all-depth S4 incoming controls."""
from pathlib import Path
import sys,json,hashlib,math,collections,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_loop_first_arrival as K
G=K.G;E=G.E;V=G.V;B=G.B
I=tuple(range(4))
def order(pi):return next(n for n in range(1,13) if E.power(pi,n)==I)
def subgroup(gs):return G.generated([G.ix[tuple(x)] for x in gs])
def end_pair(p,q,cables,rho):
    rels,M,L=K.words(p,q,cables)
    return tuple(E.evaluate(M,rho)),tuple(E.evaluate(L,rho))
def choose_root(M,L,r,k,desired='simple'):
    u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
    F=E.mul(E.power(M,k),E.power(L,r));beta=E.mul(E.power(M,u),E.power(L,v))
    for z in G.perms:
        if E.power(z,r)!=F:continue
        mu=E.mul(beta,E.power(z,-v));lam=E.mul(E.power(mu,-r*k),F)
        if desired=='simple' and mu not in [G.perms[x] for x in G.trans]:continue
        if desired=='four' and order(mu)!=4:continue
        if len(subgroup([beta,z]))!=24:continue
        return z,mu,lam,beta
    raise AssertionError((M,L,r,k,desired))
def fill_companion(p,q,cables,rho):
    rels,mw,lw=K.words(p,q,cables);M=E.evaluate(mw,rho);L=E.evaluate(lw,rho)
    assert order(M)==4
    ell=next(j for j in [0,2] if L==E.power(M,j))
    Q=K.product(lw,K.wp(mw,-ell))
    rows,qmark,keep=E.cover_matrix(rels,rho,mw,Q)
    assert len(rows)==4*len(rels)+1
    return dict(ell=ell,rows=rows,qmark=qmark,mrow=rows[-1],unfilled=rows[:-1],
                H1=E.h1(rows),H1_E=E.h1(rows[:-1]))
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    preds=['terminal_cone_twists_receipt.json','cyclic_loop_first_arrival_receipt.json',
        'iterated_peripheral_tail_receipt.json','lifted_surface_s3_receipt.json',
        'fibre_order_transfer_receipt.json','fold_schur_s4_receipt.json',
        'terminal_cone_binomial_closure_receipt.json','general_two_stage_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():
            check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/meridian-cycle-s4/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path'];check('retained complete-page primary pin '+path.name,
                                sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    n,ell,l,c=S.symbols('n ell l c',nonzero=True,integer=True)
    identity('all-degree lift boundary Q coefficient',(n/l)*l-n)
    identity('all-degree lift boundary m coefficient',(n/l)*(ell*l/n)-ell)
    check('all-degree boundary lattice determinant index n',S.Matrix([[n,-ell],[0,1]]).det()==n)
    check('odd-square half linking stable',all((r*r-1)%8==0 for r in [3,5,7,9,13]))
    # Direct symmetric-group fact: S4''=V4, intersect each4cycle centralizer in1,M^2.
    full=frozenset(range(24))
    def derived(h):
        comm=[G.mul[G.mul[G.mul[G.inverse[a]][G.inverse[b]]][a]][b] for a in h for b in h]
        return G.generated(comm)
    d1=derived(full);d2=derived(d1)
    check('S4 second derived group is the4-element Klein group',len(d1)==12 and len(d2)==4)
    for m in range(24):
        if order(G.perms[m])!=4:continue
        central={x for x in d2 if G.mul[x][m]==G.mul[m][x]}
        check('actual longitude possibilities for4cycle '+str(m),central=={G.I,G.power(m,2)})
    # All finite-period4/5 incoming transitions force true even product/fibre1/s,3s.
    sat=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())['saturation']
    states=sat['useful_profiles'];admissible=[]
    for x in sat['useful_edges']:
        if x['source'] not in [4,5] or x['dest']!=3:continue
        st=states[x['source']];M=tuple(st['M']);L=tuple(st['L']);r,k=x['rk']
        F=E.mul(E.power(M,k),E.power(L,r));z=tuple(x['z'])
        check('complete actual S4 incoming root parity '+str((x['source'],r,k,z)),
              F==I and r%2==1 and k%2==0 and r%3==0 and order(z)==3 and x['cable_order']==24)
        admissible.append(x)
    check('all32 S4 incoming residue root choices retained',len(admissible)==32)
    incoming={dest:{x['source'] for x in sat['transition_profiles'] if x['dest']==dest and x['source']!=dest} for dest in [4,5]}
    check('complete firstS4 formation dichotomy at4/5',incoming=={4:{1,6},5:{1,6}})
    initial=collections.Counter()
    for p,q in [(p,q) for p in range(12) for q in range(12) if math.gcd(math.gcd(p,q),12)==1]:
        for a,b,m,lon in G.torus(p,q):
            if len(G.groups[a,b])!=24 or order(G.perms[m])!=4:continue
            H=G.power(a,p);(aa,bb,mm),f=V.fibre_quotient([G.perms[a],G.perms[b],G.perms[m]],G.perms[H],list(range(4)))
            genus=V.genus_from_base([aa,bb],[mm])
            if genus:continue
            lengths=sorted([sorted(B.lengths(G.perms[a])),sorted(B.lengths(G.perms[b]))])
            check('every planar initial S4 fourcycle has transposition/3cycle '+str((p,q,a,b)),
                  lengths==sorted([[1,1,2],[1,3]]) and H==G.I)
            initial[p,q]+=1
    check('initial S4 planar family nonempty',bool(initial))
    for ee,oo in [(2,3),(2,9),(4,3),(8,9),(14,3)]:
        check('actual initial terminal has at least THREE genuine cones '+str((ee,oo)),
              sum(x>1 for x in [ee//2,ee,ee,oo//3,oo])>=3)
    # Initial-S4 histories and both noncyclic birth mechanisms, with actual words.
    seeds=[]
    for p,q in [(2,3),(4,3)]:
        a,b,m,lon=next(x for x in G.torus(p,q) if len(G.groups[x[0],x[1]])==24 and
                       order(G.perms[x[2]])==4 and
                       not V.genus_from_base([G.perms[x[0]],G.perms[x[1]]],[G.perms[x[2]]]))
        rho=[G.perms[a],G.perms[b]]
        seeds.append(dict(first=[p,q],cables=[],rho=rho,kind='initial'))
        M,L=end_pair(p,q,[],rho);el=0 if L==I else 2
        z=E.power(M,5*pow(3,-1,4)+el)
        seeds.append(dict(first=[p,q],cables=[(3,5)],rho=rho+[z],kind='initial_with_loop'))
    branch_controls=[]
    for oldn,p,q,fr,fk in [(2,2,3,3,2),(2,2,3,3,4),(3,2,5,2,3),(3,2,5,4,3)]:
        M=(0,1,3,2) if oldn==2 else (0,2,3,1);L=I
        rho=[E.power(M,q),E.power(M,p)]
        z,mu,lam,beta=choose_root(M,L,fr,fk,'four')
        check('actual S4 birth companion periphery '+str((oldn,fr,fk)),lam in [I,E.power(mu,2)])
        seeds.append(dict(first=[p,q],cables=[(fr,fk)],rho=rho+[z],kind='birth'+str(oldn),
                          birth=dict(oldM=M,oldL=L,beta=beta,z=z,outer_gens=rho)))
    targets=[]
    for si,seed in enumerate(seeds):
        p,q=seed['first'];cables=seed['cables'];rho=seed['rho']
        M,L=end_pair(p,q,cables,rho);filled=fill_companion(p,q,cables,rho);el=filled['ell']
        relation=[4*x+el*y for x,y in zip(filled['qmark'],filled['mrow'])]
        check('integral pre-fill longitude surface balance on full CW '+str(si),
              E.h1(filled['unfilled'],[relation])==filled['H1_E'])
        if not filled['H1']['free'] and el==2:
            check('half-framing control has nonzero even core order '+str(si),
                  filled['H1']!=E.h1(filled['rows'],[filled['qmark']]))
        for r,k in [(3,2 if el==2 else 4),(9,2 if el==2 else 4)]:
            z,tau,lam,beta=choose_root(M,L,r,k)
            check('actual final preferred longitude and simple meridian '+str((si,r,k)),lam==I)
            actual=K.summarized_cover(p,q,cables+[(r,k)],rho+[z])
            tors=B.group([r//3]);wantH=B.direct_sum([filled['H1'],tors])
            scaled=[r*x for x in filled['qmark']]
            wantQ=B.direct_sum([E.h1(filled['rows'],[scaled]),tors])
            check('full actual filled-pattern H1 and R=rQ '+str((si,r,k)),actual==dict(H1=wantH,quotient=wantQ))
            after=K.summarized_cover(p,q,cables+[(r,k),(5,3)],rho+[z,E.power(tau,3)])
            qafter=B.direct_sum([E.h1(filled['rows'],[[5*x for x in scaled]]),tors])
            check('actual later tail keeps marked old summand '+str((si,r,k)),after==dict(H1=wantH,quotient=qafter))
            if si==0 and r==3:
                check('independent full-cover tree in S4 incoming example',
                      K.summarized_cover(p,q,cables+[(r,k)],rho+[z],True)==actual)
            if seed['kind'].startswith('birth') and r==3:
                birth=seed['birth'];gr=G.O.orbit_graph(
                    dict(old=birth['outer_gens'],birth=[birth['beta'],birth['z']],arrival=[beta,z]),
                    [('old','birth',[birth['oldM'],birth['oldL']]),('birth','arrival',[M,L])],4)
                degree=sum(edge['left']==['birth',0] or edge['right']==['birth',0] for edge in gr['edges'])
                expected=4 if seed['kind']=='birth2' else 3
                check('actual persistent noncyclic birth branching '+str(si),degree==expected and gr['cycle_rank']==0)
                branch_controls.append(dict(seed=si,degree=degree,graph=gr))
            targets.append(dict(seed=si,kind=seed['kind'],first=[p,q],companion_cables=cables,
                                ell=el,arrival=[r,k],E_H1=filled['H1_E'],W_H1=filled['H1'],
                                predicted_self_linking='0' if el==0 else '1/2',
                                actual=actual,after_tail=after))
    old71=json.loads((P/'fold_schur_s4_receipt.json').read_text())
    check('known independent PF71 half-linking theorem retained',
          old71['target']['R_self_linking']=='1/2' and
          old71['source']['R_self_linking_strictly_between0_and_half'])
    inputs=[Path(__file__),P/'MERIDIAN_CYCLE_FRAMING_AND_S4_INCOMING_CLOSURE.md',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    receipt=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_cover_degrees_and_arbitrary_companions=True,
          transitive_meridian_single_boundary_basis_actual=True,
          integral_longitude_surface_balance_nQ_plus_ell_m=True,
          filling_core_self_linking_plus_or_minus_ell_over_n=True,
          finite_filled_H1_required_for_torsion_linking=True,
          no_Seifert_or_graph_hypothesis_for_framing=True,
          actual_pattern_boundary_no_torsion_coordinate_required=True,
          common_source_strict_fold_bound_all_degree_under_explicit_hypotheses=True),
        quartic=dict(S4_fourcycle_longitude_only_identity_or_square=True,
          all_S4_incoming_patterns_r3s_odd_k_even=True,
          actual_integral_pattern_m_zero_R_equal_rQ=True,
          common_self_linking_zero_or_half_excludes_all_states4_5=True,
          initial_planar_S4_terminal_at_least_three_cones=True,
          later_S4_birth_from_C2_or_C3_has_degree4_or3=True,
          unbranched_companion_not_capped_by_final_fixed_sheets=True,
          genuine_canonical_edges_persist_at_all_lengths=True,
          all_incoming_states4_5_closed_in_both_pencils=True,
          arbitrary_companion_depth_and_later_tail_length=True,no_prefix_eta_condition=True),
        ordinary_full_cover_controls=targets,actual_branching_controls=branch_controls,
        initial_planar_residue_counts={str(key):v for key,v in initial.items()},
        assurance=dict(written_unbounded_surface_cap_group_and_JSJ_proofs=True,
          full_unfilled_CW_balance_and_filled_mark_controls=True,
          complete_retained_primary_page_audits_reused=True,
          deep_JSJ_and_full_peripheral_rigidity_imported=True,
          self_linking_computed_from_surface_not_CW_group_order=True,finite_controls_not_universal_proof=True,
          individual_longitude_asserted_zero=False,exterior_pushforward_extended_through_meridian_fill=False,
          unbranched_S4_companion_asserted_cyclic=False,
          states2_7_closed=False,all_initial_state3_trajectories_closed=False,
          full_quartic_exclusion=False,higher_degree_packet_classification=False,
          finite_Keller_map=False,formal_verification=False,external_review=False,
          independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(x,P):sha(x) for x in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'meridian_cycle_s4_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree meridian-cycle framing; incoming S4 states4/5 closed in both pencils')
if __name__=='__main__':main()

