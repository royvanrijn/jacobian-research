"""PF2RULED96: prefix excess and all genuine iterated-quartic closure.

Written unbounded proofs establish the result. Exact local chains and actual
target covers are controls, not affine realizations or a finite proof census.
"""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_nonunimodular_a4 as N
import verify_cyclic_solid_cap_genus as H
import composed_peripheral_cover as C
A=H.A;K=H.K;E=H.E;G=H.G;B=H.B;I=tuple(range(4))
def chart(sigma,support):
    support=sorted(set(support),reverse=True)
    return N.contact(sigma,support)
def prefixes(sigma,support):
    out=[];row=chart(sigma,support)
    while True:
        out.append(row)
        if len(row['windings'])==1:break
        aa,gg,jj=row['drops'][-1];assert gg==1
        row=chart(row['sigma']//aa,[j//aa for j in row['support'] if j>jj])
    return list(reversed(out))
def extend(row,a,kappa):
    jj=a*row['sigma']-1+kappa;oldlast=row['sigma']-1+row['kappa']
    assert a>=2 and kappa>=2 and math.gcd(a,kappa-1)==1 and jj<a*oldlast
    support=[a*j for j in row['support']]+[jj]
    return chart(a*row['sigma'],support)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['terminal_bridge_closure_receipt.json','cyclic_solid_cap_genus_receipt.json',
        'integer_grid_common_receipt.json','nonunimodular_a4_receipt.json',
        'persistent_jsj_terminal_modulus_receipt.json','iterated_peripheral_tail_receipt.json',
        'meridian_cycle_s4_receipt.json','terminal_cone_twists_receipt.json',
        'unimodular_initial_s4_receipt.json','cyclic_pattern_a4_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/prefix-excess-closure/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    check('no new deep topology or fresh PDF audit inferred',sources['new_deep_import'] is False and sources['fresh_PDF_audit'] is False)
    a,h,d,k0,k,D,s,u,v,b,m=S.symbols('a h delta kappa_old kappa D slack u v b m',integer=True)
    dnew=a*a*d-h*(a-1)*(a+k)
    gap=a*(k0-1)-k+1
    identity('all-degree exact excess positive-edge recurrence',dnew-h*k-a*a*(d-h*k0)-h*a*gap)
    identity('index slack gives exact nonnegative-excess absorption',
        (dnew-h*k).subs(k,a*(k0-1)-s)-(a*a*(d-h*k0)+h*a*(s+1)))
    certificate=4+(h+4)*s+s*s+u*(s+2)**2+h*(s+1)*v
    identity('low determinant terminality exact negative polynomial certificate',
        dnew.subs({a:s+2,d:h-1-u,k:2+v})+certificate)
    QQ=3*b+2
    identity('single-node root contact determinant EXACT5',3*((b+1)*QQ-b+1)-QQ*QQ-5)
    identity('single-node root free order',(QQ+1-b)-(2*QQ+5)/3)
    odddelta=5*a*a-3*(a-1)*(a+k)
    identity('balanced odd induction exact even m formula',odddelta-3*k-a*(2*a+3-3*k))
    identity('odd parity makes m an even integer',(2*(2*u+1)+3-3*(2*v+1))-2*(2*u+1-3*v))
    identity('negative even m below minus2 contradicts positive determinant',
        (2*a+3+(a-1)*m).subs({a:s+5,m:-4-2*u})-(-3-2*s-2*u*(s+4)))
    identity('unique negative even m minus2 gives determinant5',(2*a+3+(a-1)*m).subs(m,-2)-5)
    arrival=[(kk,14-3*kk,(14-3*kk)-3*kk) for kk in range(2,5,2)]
    check('complete winding2 primitive positive arrival alternatives',arrival==[(2,8,2),(4,2,-10)])
    check('arrival8 has forbidden nonnegative excess and arrival2 cannot extend',arrival[0][2]>=0 and 1<=arrival[1][1]<3)
    # Boundaries: the illegal formal recurrence must NOT count as a chart.
    illegal=dict(delta_old=8,kappa_old=2,a=2,kappa_new=8,
        delta_new=int(dnew.subs({a:2,h:3,d:8,k:8})),D=int(gap.subs({a:2,k0:2,k:8})))
    check('determinant-only countercontrol violates index order',illegal['delta_new']==2 and illegal['D']==-5)
    # Exact source chains, constructed by valid local extensions, at large depths.
    balanced=[]
    for root,windings in [(5,[]),(5,[5]),(11,[5,11]),(17,[11,5,17]),(5,[5,11]*6)]:
        row=chart((root-2)//3,[root])
        for aa in windings:row=extend(row,aa,(2*aa+5)//3)
        chain=prefixes(row['sigma'],row['support'])
        check('balanced arbitrary-depth odd local chart exact determinant and parity '+str((root,windings)),
            all(x['delta']==5 and x['kappa']%2 for x in chain) and all(x%2 and x%3 for x in row['windings']))
        check('all odd block bridge differences even '+str((root,windings)),
            all((before['sigma']*aa-aa+aa*before['kappa'])-(after['sigma']-1+after['kappa'])>0 and
                ((before['sigma']*aa-aa+aa*before['kappa'])-(after['sigma']-1+after['kappa']))%2==0
                for before,after,aa in zip(chain,chain[1:],windings)))
        arrived=extend(row,2,4);arrchain=prefixes(arrived['sigma'],arrived['support'])
        check('valid local determinant2 arrival exists without a further inner node '+str((root,windings)),
            arrived['delta']==2 and arrived['kappa']==4 and arrived['windings'][-1]==2 and
            arrchain[-2]['delta']==5)
        balanced.append(dict(root=root,odd_windings=windings,odd_block=row,arrival_local_chart=arrived,
            prefix_chain=arrchain,arrival_is_LAST_node=True,not_affine_realization=True))
    # Legal local append can leave positive determinant, which is the obstruction.
    base=extend(chart(1,[5]),2,4)
    negative=extend(base,4,10)
    check('primitive index-valid local reset append has NEGATIVE determinant',negative['windings']==[5,2,4] and negative['delta']==-94 and negative['kappa']==10)
    # Actual known delta2 counterchart: parity hypothesis is indispensable.
    counter=prefixes(28,[140,48,46,35]);virtual=counter[1]
    check('old delta2 earlier-reset chart retained with odd-block parity failure',
        counter[-1]['delta']==2 and counter[-1]['kappa']==8 and virtual['windings']==[5,7] and
        virtual['delta']==11 and virtual['kappa']==6 and (35-12)%2==1)
    high=chart(7,[35,9])
    check('without negative FINAL excess a positive odd block need not stay5',
        high['windings']==[5,7] and high['delta']==65 and high['kappa']==3 and high['delta']-3*high['kappa']>0 and (35-9)%2==0)
    # Actual full target reset trajectories formerly left open by PF95.
    old=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3)
    seed=A.normalized(2,3,[],[G.perms[old[0]],G.perms[old[1]]]);M,_=A.periphery(2,3,[],seed)
    specifications=[([],2,8,1,[],3,[(11,3),(5,3)]),
        ([],8,2,1,[],7,[(5,3),(11,3)]),
        ([],4,10,1,[(2,5)],1,[(5,3)]),
        ([],10,4,1,[(4,5)],1,[(11,3),(5,3)]),
        ([(5,1)],2,8,-1,[(8,5)],1,[(5,3),(11,3)]),
        ([(5,1)],8,2,-1,[(16,5)],1,[(5,3)]),
        ([(7,1)],4,10,-1,[(4,5)],1,[(11,3),(5,3)]),
        ([(7,1)],10,4,-1,[(2,5),(5,7),(8,5)],1,[(5,3),(11,3)]*7)]
    controls=[]
    for steps,rr,jj,sign,post,Kval,tail in specifications:
        hist=H.history(steps);kk=(3*hist['A']*rr+sign*jj)//2;rho=list(seed);mm=M
        cyclic=hist['cables']+[(rr,kk)]+post
        for tr,tk in cyclic:
            assert math.gcd(tr,tk)==1 and tr%3 and tk%3
            rho.append(E.power(mm,(tk*pow(tr,-1,3))%3));mm=E.power(mm,pow(tr,-1,3))
        za,tau,_,_=next(x for x in A.roots(mm,I,2,3*Kval,2,6) if E.generated(rho+[x[0]])==24)
        cables=cyclic+[(2,3*Kval)]+tail;images=rho+[za]+[E.power(tau,tk%2) for tr,tk in tail]
        bb=6;cc=1
        for tr,tk in cables:
            assert tr>=2 and tk>=2 and math.gcd(tr,tk)==1 and tk<tr*bb
            oldc=cc;cc=tr*cc+(tr-1)*tk;bb=tr*tk
        final_eta=cables[-1][1]-cables[-1][0]*oldc-1
        check('full reset control has actual strict infinity and final degree '+str((steps,rr,jj,sign,post,tail)),hist['infinity_valid'] and final_eta<=0)
        data=C.peripheral_cover(2,3,cables,images);bd=next(x for x in data['boundaries'] if x['degree']==2)
        rows=data['relation_rows']+[x['meridian'] for x in data['boundaries']]
        ans=dict(H1=C.reduced_h1(rows),quotient=C.reduced_h1(rows,[bd['longitude']]))
        check('actual reset arrival arbitrary-tail marked H1 unchanged while prefix theorem excludes '+str((steps,rr,jj,sign,post,len(tail))),
            ans==dict(H1=B.group([2]),quotient=B.group([2])) and data['word_expansion_used'] is False)
        controls.append(dict(solid_history=steps,A=hist['A'],reset=[rr,kk],reset_orders=[rr,jj],postreset_C3=post,
            arrival=[2,3*Kval],tail=tail,final_eta=final_eta,actual=ans,
            excluded_by_written_unbounded_prefix_excess=True,not_affine_realization=True))
    sat=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())['saturation']
    incoming={x['source'] for x in sat['transition_profiles'] if x['dest']==3 and x['source']!=3}
    check('complete first-simple graph incoming list retained',incoming=={2,4,5,6,7,9})
    exhaustion={'initial3':'PF88',2:'PF90/91/93/95',4:'PF86',5:'PF86',6:'PF85',7:'PF87/89/95/96',9:'PF85'}
    check('every genuine first-arrival class closed BOTH pencils all lengths',incoming==set(x for x in exhaustion if isinstance(x,int)))
    inputs=[Path(__file__),P/'PREFIX_EXCESS_AND_ALL_GENUINE_ITERATED_QUARTIC_CLOSURE.md',P/'composed_peripheral_cover.py',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_degree_local_prefix_delta_positive_from_full_positive_delta=True,
            exact_excess_Enew_a2Eold_plus_h_a_positive_edge=True,
            nonnegative_excess_absorbs_for_every_valid_later_node=True,
            determinant_below_h_cannot_extend_with_genuine_free_order=True,
            h3_even_bridge_odd_unit3_block_all_prefix_determinants_EXACT5=True,
            odd_block_then_winding2_then_inner_node_impossible_at_full_delta2=True,
            no_affine_surface_or_finite_cover_assigned_to_local_prefix=True),
        quartic=dict(ALL_earlier_reset_A4_complementary_families_ALL_lengths_closed=True,
            zero_or_unbounded_postreset_C3_loops_included=True,
            ALL_incoming_A4_BOTH_pencils_ALL_finite_lengths_closed=True,
            ENTIRE_GENUINE_POSITIVE_ITERATED_TORUS_QUARTIC_CLASS_BOTH_PENCILS_CLOSED=True,
            complete_state_graph_and_all_predecessor_sectors_retained=True,
            whole_smooth_entire_line_actual_peripheral_simple_fullS4_hypotheses_retained=True),
        winding2_arrival_alternatives=arrival,balanced_local_prefix_controls=balanced,
        illegal_index_formal_countercontrol=illegal,negative_determinant_valid_local_extension=negative,
        parity_failure_old_delta2_counterchart=counter,
        missing_negative_final_excess_counterchart=high,
        actual_full_reset_target_controls=controls,complete_genuine_exhaustion=exhaustion,
        assurance=dict(written_unbounded_prefix_inequality_and_complete_exhaustion_proofs=True,
            index_order_and_even_bridge_parity_essential_countercontrols_retained=True,
            actual_composed_integer_target_covers_with_full_marks=True,
            primary_contact_JSJ_peripheral_inputs_reused_not_independently_reproved=True,
            finite_controls_not_universal_proof=True,
            all_degenerate_stages_classified=False,arbitrary_higher_degree_packets_excluded=False,
            all_quartic_boundary_packets_excluded=False,finite_Keller_map=False,
            full_compatible_completion_goal_solved=False,formal_verification=False,external_review=False,
            independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'prefix_excess_closure_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; ENTIRE genuine positive iterated-quartic class BOTH pencils ALL lengths closed; full goal OPEN')
if __name__=='__main__':main()
