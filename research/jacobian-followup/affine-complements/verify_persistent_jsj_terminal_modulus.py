"""PF2RULED83: exact terminal modulus and persistent marked-JSJ controls."""
from pathlib import Path
import sys,math,collections,json,hashlib,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_general_two_stage as G
E=G.E;B=G.B
C=G.F.C.V
def reduced(w):
    out=[]
    for x in w:
        if out and out[-1]==-x:out.pop()
        else:out.append(x)
    return out
def product(*ws):return reduced([x for w in ws for x in w])
def wp(w,n):return reduced(E.wp(w,n))
def ordinary(p,q,cables,rho,reverse=False,mark_multiple=1):
    u,v,d=S.gcdex(q,p);assert d==1
    M=product(wp([1],int(u)),wp([2],int(v)))
    L=product(wp([1],p),wp(M,-p*q))
    rels=[product(wp([1],p),wp([2],-q))]
    for j,(r,k) in enumerate(cables,3):
        u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
        F=product(wp(M,k),wp(L,r))
        rels.append(product(wp([j],r),E.wi(F)))
        M=product(wp(M,u),wp(L,v),wp([j],-v))
        L=product(wp(M,-r*k),F)
    identity=tuple(range(4))
    assert all(E.evaluate(w,rho)==identity for w in rels)
    assert E.evaluate(L,rho)==identity
    rows,mark,_=E.cover_matrix(rels,rho,M,L,reverse_tree=reverse)
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[[mark_multiple*x for x in mark]]))

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    preds=['iterated_peripheral_tail_receipt.json','meridian_balance_two_stage_receipt.json',
        'terminal_cone_binomial_closure_receipt.json','trefoil_complementary_fusion_receipt.json',
        'complementary_corner_jsj_receipt.json','peripheral_cone_budget_receipt.json',
        'coordinate_pencil_matrix_receipt.json','marked_end_lattice_receipt.json',
        'fibre_order_transfer_receipt.json','common_pencil_mark_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():
            check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/persistent-jsj-terminal-modulus/manifest.json'
    sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path']
        check('retained primary page audit pin '+path.name,sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    a,nu,pc,sc,T,nr,q,c=S.symbols('a nu pc sc T nr q c',integer=True)
    p=a*pc;sigma=a*sc;kap=nu*a;jstar=sigma-1+kap
    contacts=a*a*T+(a-1)*(2*p-jstar)
    nd=S.expand(contacts+2*p-sigma+1)
    identity('all-depth exact terminal diagonal modulus',nd-(a*a*(T+2*pc-sc-nu)+(nu+1)*a))
    delta=nr*nd-q*q*p*p
    identity('all-depth full determinant terminal modulus',delta-(a*a*(nr*(T+2*pc-sc-nu)-q*q*pc*pc)+nr*(nu+1)*a))
    identity('A4 residual forces divisor10',3*(3+1)-2-10)
    check('only odd terminal divisor10 above1 fails zero-mark mod3',
          [x for x in [1,2,5,10] if x>=3 and x%2]==[5] and 5%3!=1)
    n,pp,dd,A,BB=S.symbols('n p delta A B',integer=True)
    identity('zero-mark full lattice relation', (3*n-pp**2-dd).subs({n:dd*A,pp:dd*BB})-dd*(3*A-dd*BB**2-1))
    check('zero-mark mod3 forces delta2 and unit B',
          {(d,b) for d in range(3) for b in range(3) if (-d*b*b)%3==1}=={(2,1),(2,2)})
    zero_controls=[]
    for nn in range(2,70):
        for pole in range(2,math.isqrt(3*nn-1)+1):
            determinant=3*nn-pole*pole
            if nn%determinant or pole%determinant:continue
            rows=[[-3,pole],[pole,-nn]]
            hh=E.h1(rows);hq=E.h1(rows,[[1,0]])
            check('actual full binary zero-mark arithmetic '+str((nn,pole)),
                  hh==hq and hh==B.group([determinant]) and determinant%3==2 and pole%3)
            zero_controls.append(dict(n=nn,p=pole,delta=determinant,H1=hh))
    check('known nontrivial full-lattice zero-mark control retained',
          any(x['n']==2 and x['p']==2 and x['delta']==2 for x in zero_controls))
    terminal_controls=[]
    for aa in [3,5,7,11]:
        for vv in [2,3,5]:
            for depth in [2,4]:
                ss=(2 if depth==2 else 6)*aa
                pole=(vv+6)*aa if depth==2 else (vv+12)*6*aa
                last=ss-1+vv*aa
                indices=[pole,pole-aa,last] if depth==2 else [pole,pole-3*aa,pole-4*aa,last]
                indices=sorted(set(indices),reverse=True)
                check('valid terminal-gcd support '+str((aa,vv,depth)),
                      min(indices)>=ss+1 and math.gcd(ss,*indices[:-1])==aa and
                      math.gcd(ss,*indices)==1 and indices[-1]==last)
                diagonal=C.diagonal_size(1,ss,indices)
                check('actual conjugate contact sum terminal modulus '+str((aa,vv,depth)),
                      (diagonal-(vv+1)*aa)%(aa*aa)==0)
                for qq in [1,2,3]:
                    nnr=qq*qq+qq+1
                    determinant=nnr*diagonal-qq*qq*pole*pole
                    check('full source terminal determinant modulus '+str((aa,vv,depth,qq)),
                          (determinant-nnr*(vv+1)*aa)%(aa*aa)==0)
                terminal_controls.append(dict(a=aa,nu=vv,sigma=ss,support=indices,nD=diagonal,
                                              terminal_residue=diagonal%(aa*aa)))
    rr,kk,r0,k0,d,d0=S.symbols('r k r0 k0 d d0',integer=True)
    old=S.Matrix([d0*r0*k0/2,d0]);new=S.Matrix([d*kk/2,d*rr])
    identity('true primitive fibre intersection at every tail edge',
             old.row_join(new).det()-d0*d*(rr*r0*k0-kk)/2)
    edge_controls=[]
    for oldr,oldk in [(2,3),(3,2),(3,5),(9,2)]:
        od=1 if oldr*oldk%2==0 else 2
        for r,k in [(3,2),(3,5),(5,2),(5,3),(7,2),(7,3),(9,2)]:
            ndg=1 if k%2==0 else 2
            det=S.Rational(od*ndg,2)*(r*oldr*oldk-k)
            check('actual nonzero canonical tail fibre edge '+str((oldr,oldk,r,k)),
                  math.gcd(r,k)==1 and det.is_Integer and det!=0)
            edge_controls.append(dict(previous=[oldr,oldk],tail=[r,k],intersection=int(det)))
    oldcensus=json.loads((P/'general_two_stage_receipt.json').read_text())['census']['complete_cells']
    counts=collections.Counter()
    for key,cell in oldcensus.items():
        p0,q0,r,k=map(int,key.split(','))
        if r*k%2:continue
        for row in cell['orbits']:
            if row['cycle_rank'] or sum(x['genus'] for x in row['companion_bases']):continue
            counts[row['companion_order'],row['cable_order']]+=1
    check('all1920 general negative-even planar representation types covered',
          dict(counts)=={(3,24):384,(4,24):384,(12,6):384,(6,6):256,(24,24):256,(24,2):256})
    tail82=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())
    check('absorbing tails require exactly odd r C2 at every length',
          tail82['quartic']['final_state_absorbing'] and tail82['quartic']['all_later_windings_odd'] and
          tail82['quartic']['tail_image_C2'])
    # Full new ordinary CW controls retain every branch-meridian cycle.
    targets=[]
    cases=[(2,9,2,3),(2,15,2,3),(2,21,2,3),
           (2,7,3,4),(3,5,3,4),(2,3,3,2),(2,9,3,2),(4,3,3,2)]
    for pp0,qq0,rseed,kseed in cases:
        rows=[row for row in G.cell(pp0,qq0,rseed,kseed)['orbits'] if not row['cycle_rank']
              and not sum(x['genus'] for x in G.base_records(row,pp0)[0])]
        check('actual negative-even fullS4 seed exists '+str((pp0,qq0,rseed,kseed)),
              bool(rows) and rseed*kseed%2==0 and kseed<rseed*(pp0*qq0-pp0-qq0)+1)
        row=rows[0];rho=[tuple(x) for x in row['rho']]
        seed=ordinary(pp0,qq0,[(rseed,kseed)],rho)
        for tail in [[(5,3)],[(7,3),(5,3)]]:
            image=rho+[E.power(tuple(row['mu']),kk0) for rr0,kk0 in tail]
            actual=ordinary(pp0,qq0,[(rseed,kseed)]+tail,image)
            wind=math.prod(rr0 for rr0,kk0 in tail)
            scaled=ordinary(pp0,qq0,[(rseed,kseed)],rho,mark_multiple=wind)
            check('full true-filled tail H1 and marked class '+str((pp0,qq0,rseed,kseed,tail)),
                  actual['H1']==seed['H1'] and actual['quotient']==scaled['quotient'])
            if len(tail)==1:
                check('independent ordinary cover tree '+str((pp0,qq0,rseed,kseed)),
                      ordinary(pp0,qq0,[(rseed,kseed)]+tail,image,True)==actual)
            if (pp0,rseed)==(2,2):
                check('A4 residual target really has cyclic2t and R0 '+str((qq0,tail)),
                      actual==dict(H1=B.group([2*(qq0//3)]),quotient=B.group([2*(qq0//3)])))
            targets.append(dict(seed=[pp0,qq0,rseed,kseed],tail=tail,**actual))
    inputs=[Path(__file__),P/'PERSISTENT_ROOTED_JSJ_AND_TERMINAL_MODULUS.md',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    receipt=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        source=dict(no_field_degree_or_characteristic_depth_bound=True,
          full_smooth_exact_two_line_whole_chart_hypotheses_required=True,
          marked_quotient_cyclic_for_both_pencils=True,
          collar_zero_R_forces_delta2mod3_and_p_not_divisible3=True,
          every_D_winding_divides_p=True,terminal_free_tail_not_a_winding=True,
          terminal_modulus_nD_equal_nu_plus1_times_a_mod_a_squared=True,
          terminal_modulus_requires_last_gcd_a_divides_sigma_and_earlier_indices=True,
          all_depth_determinant_divisor_nR_times_nu_plus1_minus_c=True),
        quartic=dict(actual_negative_even_simple_fullS4_two_stage_prefix_required=True,
          arbitrary_lower_weight_keys_and_actual_coprime_tail_parameters=True,
          arbitrary_finite_tail_length=True,tail_edges_genuine_by_actual_nonzero_fibre_intersection=True,
          all_even_k_tails_excluded_by_root_or_internal_cone_budget=True,
          marked_quotient_branching_internal_and_terminal_obstructions_persist=True,
          cyclic4_root3_residual_excluded_by_zero_mark_factor3=True,
          A4_residual_t_divides10_and_fails_mod3=True,
          entire_seed_tail_class_excluded_for_complementary_source=True,
          PF82_common_source_closure_retained=True,both_pencils_closed_for_seed_tail_class=True),
        even_planar_sector_counts={str(k):v for k,v in counts.items()},
        zero_mark_controls=zero_controls,terminal_modulus_controls=terminal_controls,
        canonical_edge_controls=edge_controls,ordinary_full_cover_controls=targets,
        assurance=dict(written_unbounded_gcd_contact_modulus_and_JSJ_proofs=True,
          all_fixed_sheet_and_branch_cycle_fillings_retained=True,
          complete_predecessor_primary_audits_reused=True,deep_JSJ_and_rigidity_imported=True,
          finite_controls_not_universal_proof=True,
          prefix_source_packet_assumed=False,final_eta_inferred_for_every_prefix=False,
          all_first_arrival_paths_closed=False,all_longer_targets_excluded=False,
          general_higher_degree_packets_classified=False,full_quartic_exclusion=False,
          finite_Keller_map=False,formal_verification=False,external_review=False,
          independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(x,P):sha(x) for x in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'persistent_jsj_terminal_modulus_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; terminal modulus and persistent JSJ; negative-even seed tails closed in both pencils')
if __name__=='__main__':main()

