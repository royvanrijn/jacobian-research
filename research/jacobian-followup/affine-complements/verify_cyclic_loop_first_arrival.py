"""PF2RULED84: cyclic cable models, all-length DAG bounds and true-cover controls."""
from pathlib import Path
import sys,math,json,hashlib,collections,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_general_two_stage as G
import verify_persistent_jsj_terminal_modulus as J
E=G.E;V=G.V;B=G.B
product=J.product;wp=J.wp
def words(p,q,cables):
    u,v,d=S.gcdex(q,p);assert d==1
    M=product(wp([1],int(u)),wp([2],int(v)))
    L=product(wp([1],p),wp(M,-p*q))
    rels=[product(wp([1],p),wp([2],-q))]
    for j,(r,k) in enumerate(cables,3):
        u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
        F=product(wp(M,k),wp(L,r))
        rels.append(product(wp([j],r),E.wi(F)))
        M=product(wp(M,u),wp(L,v),wp([j],-v));L=product(wp(M,-r*k),F)
    return rels,M,L
def full(data,rho,reverse=False):
    rels,M,L=data;identity=tuple(range(len(rho[0])))
    assert all(E.evaluate(w,rho)==identity for w in rels) and E.evaluate(L,rho)==identity
    rows,mark,_=E.cover_matrix(rels,rho,M,L,reverse_tree=reverse)
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[mark]))
def summarized_cover(p,q,cables,rho,reverse=False,raw=False):
    """Exact lifted cellular paths; composition avoids expanding peripheral words."""
    d=len(rho[0]);n=len(rho);width=d*n;ident=tuple(range(d))
    unit=(ident,[[0]*width for s in range(d)])
    def compose(A,B):
        da,fa=A;db,fb=B
        return (tuple(db[da[s]] for s in range(d)),
                [[x+y for x,y in zip(fa[s],fb[da[s]])] for s in range(d)])
    def invert(A):
        dest,flows=A;back=E.inv(dest)
        return back,[[-x for x in flows[back[s]]] for s in range(d)]
    def power(A,k):
        if k<0:return power(invert(A),-k)
        out=unit
        while k:
            if k%2:out=compose(out,A)
            A=compose(A,A);k//=2
        return out
    def prod(*args):
        out=unit
        for A in args:out=compose(out,A)
        return out
    generators=[]
    for g,pi in enumerate(rho):
        flows=[[0]*width for s in range(d)]
        for s in range(d):flows[s][g*d+s]=1
        generators.append((E.inv(pi),flows))
    u,v,gg=S.gcdex(q,p);assert gg==1
    M=prod(power(generators[0],int(u)),power(generators[1],int(v)))
    L=prod(power(generators[0],p),power(M,-p*q))
    rels=[prod(power(generators[0],p),power(generators[1],-q))]
    for j,(r,k) in enumerate(cables,2):
        u,v,gg=S.gcdex(r,k);assert gg==1;u=int(u);v=-int(v)
        F=prod(power(M,k),power(L,r))
        rels.append(prod(power(generators[j],r),invert(F)))
        M=prod(power(M,u),power(L,v),power(generators[j],-v))
        L=prod(power(M,-r*k),F)
    assert all(A[0]==ident for A in rels) and L[0]==ident
    edges=[(g,s) for g in range(n) for s in range(d)];parent=list(range(d));tree=set()
    def root(s):
        while parent[s]!=s:s=parent[s]
        return s
    for g,s in (list(reversed(edges)) if reverse else edges):
        t=E.inv(rho[g])[s]
        if root(s)!=root(t):parent[root(s)]=root(t);tree.add((g,s))
    assert len(tree)==d-1
    keep=[e for e in edges if e not in tree];indices=[g*d+s for g,s in keep]
    project=lambda flow:[flow[i] for i in indices]
    rows=[project(A[1][s]) for A in rels for s in range(d)]
    seen=set();ram=None
    for s in range(d):
        if s in seen:continue
        cyc=[];t=s
        while t not in cyc:cyc.append(t);seen.add(t);t=M[0][t]
        closed=power(M,len(cyc));assert closed[0][s]==s
        rows.append(project(closed[1][s]))
        if len(cyc)==2:ram=s
    mark=project(L[1][0 if ram is None else ram])
    if raw:return rows,mark,keep
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[mark]))
def cyclic_rho(p,q,cables,M):
    wind=math.prod(r for r,k in cables)
    out=[E.power(M,q*wind),E.power(M,p*wind)]
    for i,(r,k) in enumerate(cables):
        out.append(E.power(M,k*math.prod(rr for rr,kk in cables[i+1:])))
    return out
def pattern(r,k,M,z):
    u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
    F=wp([1],k);meridian=product(wp([1],u),wp([2],-v))
    longitude=product(wp(meridian,-r*k),F)
    return full(([product(wp([2],r),E.wi(F))],meridian,longitude),[M,z])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    preds=['persistent_jsj_terminal_modulus_receipt.json','iterated_peripheral_tail_receipt.json',
        'meridian_balance_two_stage_receipt.json','general_two_stage_receipt.json',
        'cyclic_companion_receipt.json','cyclic_parity_receipt.json','fibre_order_transfer_receipt.json',
        'terminal_cone_binomial_closure_receipt.json','seifert_base_genus_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():
            check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/cyclic-loop-first-arrival/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path'];check('retained complete-page primary pin '+path.name,
                                 sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    cyclic_controls=[]
    for n in [2,3,4,5,7,8,11,12]:
        M=tuple((i+1)%n for i in range(n))
        for r,k in [(2,3),(3,2),(3,5),(5,3),(5,4),(7,2)]:
            if math.gcd(r,n)!=1:continue
            for ell in ([0,2] if n==4 else [0,1]):
                u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
                z=E.power(M,(k*pow(r,-1,n)+ell)%n)
                L=E.power(M,ell);F=E.mul(E.power(M,k),E.power(L,r))
                beta=E.mul(E.power(M,u),E.power(L,v))
                mu=E.mul(beta,E.power(z,-v));lam=E.mul(E.power(mu,-r*k),F)
                (zz,bb,tt),f=V.fibre_quotient([z,beta,mu],F,list(range(n)))
                g=math.gcd(n,k+ell*r)
                tag=str((n,r,k,ell))
                check('exact cyclic true peripheral maps '+tag,
                      mu==E.power(M,pow(r,-1,n)) and lam==E.power(mu,ell*r*r))
                check('actual cyclic base fibre degrees and repeated cones '+tag,
                      f==n//g and len(zz)==g and zz==tuple(range(g)) and
                      V.cycles(bb)==V.cycles(tt)==1 and V.genus_from_base([zz],[bb,tt])==0)
                cyclic_controls.append(dict(n=n,r=r,k=k,ell=ell,base_degree=g,fibre_degree=f,cone_orders=[r]*g))
    old82=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())
    sat=old82['saturation'];profiles=sat['transition_profiles']
    adjacency=collections.defaultdict(list)
    for x in profiles:
        if x['source']!=x['dest']:adjacency[x['source']].append((x['dest'],x['cable_order'] not in [1,2,3,4]))
    visiting=set();memo={}
    def depth(i):
        assert i not in visiting,'nontrivial state cycle'
        if i in memo:return memo[i]
        visiting.add(i)
        ans=max([0]+[1+depth(j) for j,cost in adjacency[i]])
        visiting.remove(i);memo[i]=ans;return ans
    check('all nonself transitions form DAG with at most three changes',max(depth(i) for i in range(10))==3)
    memo2={}
    def cost(i):
        if i in memo2:return memo2[i]
        memo2[i]=max([0]+[int(w)+cost(j) for j,w in adjacency[i]])
        return memo2[i]
    check('at most TWO noncyclic cable levels at arbitrary length',max(cost(i) for i in range(10))==2)
    states=sat['useful_profiles'];loop_count=0
    check('sharp bound TWO noncyclic levels INCLUDING initial torus',
          max(int(states[i]['order'] not in [1,2,3,4])+cost(i) for i in range(10))==2)
    for row in sat['useful_edges']:
        if row['source']!=row['dest']:continue
        source=states[row['source']];M=tuple(source['M']);L=tuple(source['L'])
        n=next(j for j in range(1,13) if E.power(M,j)==tuple(range(4)))
        r,k=row['rk'];ell=2 if row['source']==5 else 0
        z=tuple(row['z']);u,v=next((u,v) for u in range(12) for v in range(12) if (r*u-k*v)%12==1)
        F=E.mul(E.power(M,k),E.power(L,r));beta=E.mul(E.power(M,u),E.power(L,v))
        mu=E.mul(beta,E.power(z,-v));lam=E.mul(E.power(mu,-r*k),F)
        check('every quartic loop is actual cyclic model '+str((row['source'],r,k,z)),
              math.gcd(r,n)==1 and row['cable_order']==n and L==E.power(M,ell) and
              z in [E.power(M,j) for j in range(n)] and mu==E.power(M,pow(r,-1,n)) and
              lam==E.power(mu,ell*r*r))
        loop_count+=1
    check('complete self-loop root choices retained',loop_count==688)
    a,kap,pc,sc,T,nr,q=S.symbols('a kappa pc sc T nr q',integer=True)
    p=a*pc;ss=a*sc;last=ss-1+kap
    nd=S.expand(a*a*T+(a-1)*(2*p-last)+2*p-ss+1)
    identity('general terminal modulus without multiple-tail assumption',
             nd-(a*a*(T+2*pc-sc)+a-(a-1)*kap))
    identity('general full terminal determinant modulus',
             nr*nd-q*q*p*p-(a*a*(nr*(T+2*pc-sc)-q*q*pc*pc)+nr*(a-(a-1)*kap)))
    check('even first pair forces only terminal5 or25',
          [d for d in [1,5,25] if d>=3 and d%2]==[5,25])
    residues={x*x%13 for x in range(13)}
    check('terminal5 zero mark impossible modulo13',(-5)%13 not in residues)
    check('terminal25 forces winding-product square3 modulo13',(-pow(9*25,-1,13))%13==3)
    terminal_candidates=[]
    for aa in range(5,26,2):
        for kk in range(3,aa,2):
            D=13*(aa-(aa-1)*kk)-1
            if D%(aa*aa)==0 and math.gcd(aa,kk)==math.gcd(aa,kk-1)==1:
                terminal_candidates.append([aa,kk,D//(aa*aa)])
    check('all forced larger-terminal candidates excluded after divisor bound',terminal_candidates==[])
    m=S.symbols('m',integer=True)
    identity('forced divisor a-minus1 certificate',
             (m*a*a+13*a-1)-(m+12)-(a-1)*(m*(a+1)+13))
    for degree in [2,3,4]:
        shift=tuple((i+1)%degree for i in range(degree))
        rho=cyclic_rho(2,3,[],shift);rels,meridian,longitude=words(2,3,[])
        check('exact compressed CW comparison for cyclic meridian degree'+str(degree),
              summarized_cover(2,3,[],rho,raw=True)==E.cover_matrix(rels,rho,meridian,longitude))
    # True full branched-cover controls, not realizations of polynomial curves.
    targets=[]
    cases=[(3,2,5,[(3,2)],4,3),(3,2,5,[(5,2),(3,2)],4,9),
           (4,2,3,[(2,3)],3,4),(4,3,5,[(3,5)],3,4),
           (4,3,5,[(3,5),(5,3)],9,4),(4,2,25,[(9,5)],3,4)]
    for n0,p0,q0,prefix,r,k in cases:
        co=n0
        rows=[x for x in G.cell(p0,q0,r,k)['orbits'] if x['companion_order']==co and x['cable_order']==24]
        check('necessary cyclic arrival root exists '+str((n0,p0,q0,r,k)),bool(rows))
        row=rows[0];M=tuple(row['M']);z=tuple(row['rho'][2])
        rho=cyclic_rho(p0,q0,prefix,M);data=words(p0,q0,prefix)
        check('deep actual cyclic companion periphery '+str((n0,p0,q0,prefix)),
              E.evaluate(data[1],rho)==M and E.evaluate(data[2],rho)==tuple(range(4)))
        # Compute the branched cyclic companion homology in its own connected degree.
        shift=tuple((i+1)%n0 for i in range(n0))
        comp=summarized_cover(p0,q0,prefix,cyclic_rho(p0,q0,prefix,shift))['H1']
        pat=pattern(r,k,M,z)
        check('actual unknotted pattern has zero R '+str((n0,r,k)),pat['H1']==pat['quotient'])
        expected=B.direct_sum([pat['H1'],comp])
        actual=summarized_cover(p0,q0,prefix+[(r,k)],rho+[z])
        check('deep ordinary all-cycle-filled replacement and R0 '+str((n0,p0,q0,prefix,r,k)),
              actual==dict(H1=expected,quotient=expected))
        check('independent full arrival spanning tree '+str((n0,p0,q0,prefix,r,k)),
              summarized_cover(p0,q0,prefix+[(r,k)],rho+[z],True)==actual)
        if len(targets)==0:
            rels,meridian,longitude=words(p0,q0,prefix+[(r,k)])
            check('compressed paths equal every expanded CW row and actual R mark',
                  summarized_cover(p0,q0,prefix+[(r,k)],rho+[z],raw=True)==
                  E.cover_matrix(rels,rho+[z],meridian,longitude))
        tau=tuple(row['mu']);tail=(5,3)
        after=summarized_cover(p0,q0,prefix+[(r,k),tail],rho+[z,E.power(tau,3)])
        check('actual further C2 tail preserves cyclic arrival R0 '+str((n0,p0,q0,prefix,r,k)),
              after==actual)
        targets.append(dict(n=n0,first=[p0,q0],companion_cables=prefix,arrival=[r,k],
                            pattern=pat,branched_companion=comp,arrival_groups=actual,after_tail=after))
    inputs=[Path(__file__),P/'CYCLIC_LOOP_GEOMETRY_AND_FIRST_ARRIVAL_REDUCTION.md',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    receipt=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_degrees_and_cyclic_level_counts=True,
          actual_cyclic_root_and_preferred_longitude_shift=True,gcd_r_n_one_required_for_single_boundary_model=True,
          fibre_degree_n_over_gcd_n_k_plus_ell_r=True,annulus_g_identical_order_r_cones=True,
          nonidentity_fibres_and_fixed_sheet_orbits_retained=True,
          arbitrary_companion_cyclic_integral_replacement_retained=True,
          exact_composed_cellular_path_matrices_without_word_expansion=True,
          all_depth_terminal_modulus_arbitrary_free_tail=True),
        quartic=dict(at_most_two_noncyclic_cable_levels=True,at_most_two_with_initial_torus=True,
          no_cyclic_level_geometric_contraction=True,
          cyclic3_arrival_any_companion_depth_any_tail_both_pencils_excluded=True,
          cyclic4_arrival_any_nonempty_later_tail_both_pencils_excluded=True,
          cyclic4_arrival_winding_at_least9_both_pencils_excluded=True,
          no_prefix_eta_condition_for_new_arrival_closures=True,
          remaining_cyclic4_arrival_final_winding3=True,
          residual_companion_all_windings_and_coefficients_odd=True,
          residual_even_first_pair_only_2_25=True,
          residual_odd_terminal_winding_smaller_and_coprime6_13=True),
        cyclic_level_controls=cyclic_controls,nonself_state_change_bound=3,noncyclic_cable_bound=2,
        quartic_self_loop_choices=loop_count,larger_terminal_forced_candidates=terminal_candidates,
        ordinary_deep_companion_controls=targets,
        assurance=dict(written_unbounded_cyclic_geometry_and_first_arrival_proofs=True,
          finite_DAG_certificate_after_exact_period_proof=True,
          exact_compressed_cellular_paths_cross_checked_against_expanded_CW=True,
          all_actual_meridian_cycle_fillings_retained=True,retained_primary_audits_reused=True,
          deep_JSJ_and_peripheral_rigidity_imported=True,finite_controls_not_universal_proof=True,
          loops_geometrically_deleted=False,prefix_source_packet_assumed=False,
          necessary_action_declared_affine_realization=False,polynomial_curve_realization_of_controls_asserted=False,
          residual_2_25_closed=False,all_final_cyclic4_winding3_closed=False,
          all_noncyclic_arrivals_closed=False,full_quartic_exclusion=False,
          higher_degree_packet_classification=False,finite_Keller_map=False,formal_verification=False,
          external_review=False,independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(x,P):sha(x) for x in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'cyclic_loop_first_arrival_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree cyclic levels and all-length cyclic first-arrival reductions')
if __name__=='__main__':main()

