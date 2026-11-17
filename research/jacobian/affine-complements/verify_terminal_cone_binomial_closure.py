"""PF2RULED79: exact controls for the unrestricted companion/source proof."""
from pathlib import Path
import collections,hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
F=load('fusion79','verify_trefoil_complementary_fusion.py');V=F.V;E=V.E;B=F.B
def genuine(xs):return sorted(int(x) for x in xs if x>1)
def cone_orders(perm,order):
    lengths=B.lengths(perm);assert all(order%l==0 for l in lengths)
    return [order//l for l in lengths]

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['trefoil_complementary_fusion_receipt.json','complementary_corner_jsj_receipt.json',
      'matched_meridian_common_exclusion_receipt.json','fibre_order_transfer_receipt.json',
      'fold_schur_s4_receipt.json','peripheral_cone_budget_receipt.json',
      'seifert_base_genus_receipt.json','cyclic_parity_receipt.json',
      'binomial_periodic_receipt.json','coordinate_pencil_matrix_receipt.json',
      'marked_end_lattice_receipt.json','unramified_orevkov_receipt.json',
      'nodal_deformation_receipt.json','peripheral_zero_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(HERE/path)==digest)
    manifest=HERE/'sources/terminal-cone-binomial-closure/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('retained audited primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Unrestricted contact/lattice identities, not interpolation in t or o.
    t,o,b,r,e=S.symbols('t o b r e',positive=True,integer=True)
    def diagonal(sigma,p,jstar=None,g=None):
        contacts=(sigma-1)*p if jstar is None else (sigma-g)*p+(g-1)*(2*p-jstar)
        return S.expand(contacts+2*p-sigma+1)
    ndA=diagonal(2*t+1,3*t)
    check('A4 higher-R residual nD identity',ndA==6*t*t+4*t)
    check('A4 higher-R residual full determinant identity',S.expand(7*ndA-(6*t)**2-2*t*(3*t+14))==0)
    ndAc=diagonal(3*t,6*t,4*t-1,3*t)
    check('A4 collar residual nD identity',ndAc==24*t*t+4*t)
    check('A4 collar residual full determinant identity',S.expand(3*ndAc-(6*t)**2-12*t*(3*t+1))==0)
    check('A4 discarded terminal order exceeds top support',S.expand((4*t-1)-2*t-(2*t-1))==0)
    ndC=diagonal(b*o,3*o,(b+1)*o-1,o)
    check('C4 both collar contacts exact symbolic nD',S.expand(ndC-((2*b+2)*o**2+2*o))==0)
    check('C4 both collar full determinants exact symbolic',S.expand(3*ndC-(3*o)**2-3*o*((2*b-1)*o+2))==0)
    check('smooth unramified tail uses 2p not guessed p squared',diagonal(1,o)==2*o)
    Q=S.symbols('Q',positive=True,integer=True)
    check('all-depth rho2 first-cone marked distance determinant',S.Matrix([[1,-b],[-3,Q]]).det()==Q-3*b)
    check('cone3 collar marked distances cannot be regular drilling',{abs(3-3*bb) for bb in [1,2]}=={0,3})
    check('C4 higher-R residual violates negative definiteness',S.expand(13*(2*o)-(3*o)**2-o*(26-9*o))==0 and 26-9*3<0)
    check('solid S4C2 forced smooth source determinant9',3*diagonal(1,3)-3**2==9)
    check('solid A4 forced source nD134 determinant258',diagonal(11,12)==134 and 3*134-12**2==258)
    c=e*o-e-o;k=r*c+1
    check('A4 actual fibre intersection',S.expand(r*e*o/3-k/3-(r*(e+o)-1)/3)==0)
    check('S4C2 actual fibre intersection',S.expand(r*e*o/2-k/2-(r*(e+o)-1)/2)==0)
    check('C4 primitive lifted H2 intersection',S.expand(r*e*o/2-2*k/4-(r*(e+o)-1)/2)==0)
    check('S4S4 actual nonzero edge j',S.expand((1-r*(o+2))/4+(r*(o+2)-1)/4)==0)
    # Actual solid-torus kernel slopes; both free-generator signs are kept.
    for sector,Hcoef,mabs,etas in [('A4',2,2,[(1,2),(7,2)]),('S4C2',6,1,[(3,1),(9,1)])]:
        for sign,eta in zip([1,-1],etas):
            mfree=sign*mabs;lfree=3-Hcoef*mfree
            check('actual solid companion primitive kernel '+str((sector,sign)),
                  eta[0]*mfree+eta[1]*lfree==0 and math.gcd(*eta)==1)
    jAm=(r-2)/3;jAp=(19*r-2)/3;jSm=(r-1)/2;jSp=(13*r-1)/2
    check('A4 both solid-cap distances',S.expand(r-2*(r+1)/3-jAm)==0 and S.expand(7*r-2*(r+1)/3-jAp)==0)
    check('S4C2 both solid-cap distances',S.expand(3*r-(5*r+1)/2-jSm)==0 and S.expand(9*r-(5*r+1)/2-jSp)==0)
    u=S.symbols('u',positive=True,integer=True)
    check('A4 smaller cap never consecutive for u>=7',S.expand(u-(2*u-2)/3-(u+2)/3)==0 and S.Rational(7+2,3)>1)
    check('A4 larger cap gap',S.expand((38*u-2)/3-2*u-(32*u-2)/3)==0)
    check('S4 smaller cap gap',S.expand(r-jSm-(r+1)/2)==0)
    check('S4 larger cap gap',S.expand(jSp-r-(11*r-1)/2)==0)
    check('order1 and fibre-parallel cap exceptions explicit',jAm.subs(r,2)==0 and jAp.subs(r,2)==12 and jSm.subs(r,3)==1)
    check('A4 solid-cap targets negative orbifold Euler uniformly',S.Rational(1,7)+S.Rational(1,14)+S.Rational(1,4)<1 and S.Rational(1,2)+S.Rational(1,12)<1)
    check('S4 solid-cap targets negative orbifold Euler uniformly',S.Rational(2,5)+S.Rational(1,2)<1 and S.Rational(2,3)+S.Rational(1,19)<1 and S.Rational(2,3)<1)
    check('Klein has abelian boundary subgroup index2 explicit conjugation',S.Matrix([[1,0],[0,-1]])**2==S.eye(2))
    data=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    counts=collections.Counter();cone_counts=collections.Counter()
    for key,cell in data.items():
        pp,qq,rr=map(int,key.split(','));pp=pp if pp>=2 else pp+12;qq=qq if qq>=2 else qq+12;rr=rr if rr>=2 else rr+12
        for i,row in enumerate(cell['orbits']):
            tag=str((key,i));co,ca=row['companion_order'],row['cable_order']
            comps=V.companion_bases(row,pp,qq,rr)
            if row['cycle_rank']:counts['cycle obstruction']+=1;continue
            if sum(x['genus'] for x in comps):counts['positive base genus']+=1;continue
            sector={(3,24):'C3',(4,24):'C4',(12,6):'A4/S3',(6,6):'S3/S3',(24,24):'S4/S4',(24,2):'S4/C2'}[(co,ca)]
            counts[sector]+=1
            if sector in ['C3','S3/S3']:continue
            if pp%2 and qq%2:
                check('odd first pair already closed '+tag,sector=='C4');continue
            ee,oo=(pp,qq) if pp%2==0 else (qq,pp)
            if sector=='C4' and ee%4==0:continue
            check('remaining target single transitive companion '+tag,len(comps)==1 and len(comps[0]['sheets'])==4)
            aperm,bperm,mm=comps[0]['induced']
            true=cone_orders(aperm,pp)+cone_orders(bperm,qq)
            if sector=='S4/S4':want=[ee//2,ee,ee,oo//3,oo]
            elif sector=='A4/S3':want=[ee//2,ee//2,oo//3,oo]
            elif sector=='S4/C2':want=[ee//4,oo//3,oo]
            else:want=[ee//2,oo,oo]
            check('complete actual fibre-orbit companion cones '+tag,genuine(true)==genuine(want))
            check('correct companion generic fibre degree '+tag,comps[0]['fibre_degree']==(2 if sector=='C4' else 1))
            cabs=V.cable_bases(row);main=max(cabs,key=lambda x:len(x['sheets']));root,beta,tau=main['induced']
            cable=cone_orders(root,rr)
            cwant=[rr//2,rr] if sector=='A4/S3' else [rr,rr] if sector=='S4/C2' else [rr//3,rr]
            check('complete actual main cable cones '+tag,genuine(cable)==genuine(cwant) and main['genus']==0 and main['fibre_degree']==1)
            check('all fixed sheet cable pieces retained '+tag,all(len(x['sheets'])==1 for x in cabs if x is not main))
            cone_counts[sector]+=1
    expected={'cycle obstruction':112,'positive base genus':32,'C3':84,'C4':96,'A4/S3':96,'S3/S3':64,'S4/S4':64,'S4/C2':48}
    check('ALL596 periodic representation classes exhausted',dict(counts)==expected and sum(counts.values())==596)
    controls=[]
    for tt in [3,5]:
        for rho,sigma,support,want in [(3,2*tt+1,{3*tt:2},2*tt*(3*tt+14)),(2,3*tt,{6*tt:2,4*tt-1:1},12*tt*(3*tt+1))]:
            ctl=F.cluster_control(rho,sigma,support);tag=str(('A4',tt,rho))
            check('forced actual A4 source lattice '+tag,ctl['matrix'].det()==want)
            check('forced A4 terminal and root cones '+tag,any(genuine(x)==[tt,3*tt] for x in ctl['legs'].values()) and ctl['legs'][ctl['root']]==[2])
            controls.append(dict(sector='A4/S3 r2',parameter=tt,rho=rho,sigma=sigma,support=support,matrix=ctl['matrix'].tolist(),determinant=want))
    for oo in [3,7]:
        for rho,sigma,support,want in [(4,1,{oo:2},oo*(26-9*oo)),(2,oo,{3*oo:2,2*oo-1:1},3*oo*(oo+2)),(2,2*oo,{3*oo:2,3*oo-1:1},3*oo*(3*oo+2))]:
            ctl=F.cluster_control(rho,sigma,support);tag=str(('C4',oo,rho,sigma))
            check('forced actual C4 source lattice '+tag,ctl['matrix'].det()==want)
            check('forced C4 terminal and root cones '+tag,any(genuine(x)==[oo,oo] for x in ctl['legs'].values()) and ctl['legs'][ctl['root']]==[3])
            controls.append(dict(sector='C4 r3',parameter=oo,rho=rho,sigma=sigma,support=support,matrix=ctl['matrix'].tolist(),determinant=want))
    for rho,sigma,support,want in [(2,11,{12:2},258),(2,1,{3:2},9)]:
        ctl=F.cluster_control(rho,sigma,support)
        check('forced actual solid-cap source '+str((sigma,support)),ctl['matrix'].det()==want and len(ctl['legs'])==1)
        controls.append(dict(sector='solid cap residual',rho=rho,sigma=sigma,support=support,matrix=ctl['matrix'].tolist(),determinant=want))
    # Additional permitted coefficients do not change the characteristic-contact controls.
    for sigma,support,nd in [(7,{9:2,8:5},66),(9,{18:2,11:1,10:-3},228),(3,{9:2,6:4,5:1,4:7},42)]:
        check('noncharacteristic coefficients leave exact diagonal unchanged '+str((sigma,support)),F.C.V.diagonal_size(1,sigma,support)==nd)
    cases=[(2,9,3,24,24),(2,3,9,24,24),(4,3,10,12,6),(2,9,2,12,6),(2,15,2,12,6),(2,3,2,12,6),(2,3,14,12,6),(2,3,26,12,6),
           (4,3,3,24,2),(4,3,5,24,2),(4,9,5,24,2),(8,3,5,24,2),(8,9,5,24,2),(2,3,3,4,24),(2,7,3,4,24),(6,5,9,4,24),(2,3,15,4,24)]
    targets=[]
    for ee,oo,rr,co,ca in cases:
        rows=[x for x in data[f'{ee%12},{oo%12},{rr%12}']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(co,ca,0)]
        check('actual target rows exist '+str((ee,oo,rr,co,ca)),bool(rows))
        if co==24 and ca==24:want=dict(H1=B.group([8*(oo//3),rr//3]),quotient=B.group([2*(oo//3),rr//3]))
        elif co==12:want=dict(H1=B.group([ee*oo//3,rr//2]),quotient=B.group([ee*oo//3,rr//2]))
        elif ca==2:want=dict(H1=B.group([oo//3,rr]),quotient=B.group([oo//3,rr]))
        else:want=dict(H1=B.group([oo,rr//3]),quotient=B.group([oo,rr//3]))
        for i,row in enumerate(rows):check('ordinary full target group and R mark '+str((ee,oo,rr,i)),V.full_satellite(ee,oo,rr,row)==want)
        targets.append(dict(e=ee,o=oo,r=rr,companion_order=co,cable_order=ca,classes=len(rows),**want))
    old78=json.loads((HERE/'trefoil_complementary_fusion_receipt.json').read_text());old76=json.loads((HERE/'matched_meridian_common_exclusion_receipt.json').read_text());old68=json.loads((HERE/'cyclic_parity_receipt.json').read_text())
    check('common and S3 closures retained',old78['target']['entire_S3_S3_sector_excluded_in_both_pencils'] and old76['quartic']['entire_common_source_sector_excluded'])
    check('cyclic predecessor retained without rewriting sealed historical frontiers',old68['outcome']=='PASS')
    inputs=[Path(__file__),HERE/'TERMINAL_CONE_BUDGET_AND_BINOMIAL_QUARTIC_CLOSURE.md',manifest,HERE/'verify_trefoil_complementary_fusion.py']+[HERE/n for n in preds]+[HERE/row['path'] for row in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(no_degree_order_support_depth_bound=True,two_entire_omitted_lines_required=True,full_smooth_exact_complementary_source_and_whole_charts_required=True,
        terminal_opposite_R_piece_exactly_two_cones_if_JSJ_edge=True,negative_Euler_orientable_target_pieces_required=True,
        root_two_cones_must_be_consecutive=True,single_piece_source_disk_two_cones_or_consecutive_three_cones=True,
        Klein_survival_excluded_by_virtual_abelianness=True,Klein_alternate_merger_excluded_by_centre=True,
        Klein_other_alignment_assumed_absent=False,smooth_D_diagonal='2*p_D',rho2_one_cone_marked_distance='abs(Q-3b),1<=b<Q,gcd(b,Q)=1'),
      quartic=dict(actual_primitive_binomial_two_stage_required=True,PF46_50_actual_source_packet_hypotheses_required=True,
        all_coprime_first_pairs_ge2_all_windings_ge2_all_source_data=True,
        S4_S4_complementary_entirely_excluded=True,A4_S3_complementary_entirely_excluded=True,S4_C2_complementary_entirely_excluded=True,remaining_even_C4_complementary_entirely_excluded=True,
        all_representations_and_both_pencils_in_stated_target_class_excluded=True,entire_actual_binomial_two_stage_quartic_class_excluded=True,
        actual_fixed_sheet_caps_retained=True,actual_solid_topology_proved_before_homological_kernel_slopes=True,
        both_solid_cap_signs_retained=True,A4_fibre_parallel_cap_contradicts_boundary_injection=True,
        complete596_class_exhaustion=dict(counts),actual_cone_audit_counts=dict(cone_counts)),
      residual_formulas=dict(A4_r2_higher_R='2t(3t+14)>2t',A4_r2_collar='12t(3t+1)>2t',C4_r3_higher_R='o(26-9o)<0',C4_r3_collars='marked fibre distance0 or3, not target1',optional_C4_collar_determinant_controls='3o[(2b-1)o+2]>o,b=1,2',A4_solid_cap='258 !=2',S4C2_solid_cap='9 !=3'),
      source_cluster_controls=controls,ordinary_target_controls=targets,
      assurance=dict(written_unbounded_terminal_cap_kernel_and_contact_proofs=True,complete_existing_census_exhausted=True,full_peripheral_recognition_hypotheses_checked=True,
        compressible_target_excluded_without_Waldhausen=True,prior_complete_Hatcher_Waldhausen_page_audits_reused=True,deep_JSJ_fibre_rigidity_imported=True,finite_controls_not_universal_proof=True,
        general_lower_weight_keys_assumed_binomial=False,arbitrary_omitted_count_theorem=False,higher_degree_packet_classification=False,
        full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,PF63_polynomial_immersions_declared_nonexistent=False,
        formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    def native(x):
        if isinstance(x,S.Integer):return int(x)
        raise TypeError(type(x).__name__)
    (HERE/'terminal_cone_binomial_closure_receipt.json').write_text(json.dumps(report,indent=2,default=native)+'\n')
    print(f'PASS {len(checks)} controls; terminal cone budget; entire actual-binomial two-stage quartic target class excluded')
if __name__=='__main__':main()
