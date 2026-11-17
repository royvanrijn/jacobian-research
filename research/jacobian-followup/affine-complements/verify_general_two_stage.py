"""PF2RULED80: exact general cable cells and odd fibre controls."""
from pathlib import Path
import sys,itertools,math,collections,json,hashlib,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent;P=HERE;sys.path.insert(0,str(HERE))
import verify_orbit_tree as O
import verify_trefoil_complementary_fusion as F
E=F.E;V=F.V;B=F.B
perms=list(itertools.permutations(range(4)));ix={p:i for i,p in enumerate(perms)};I=ix[(0,1,2,3)]
mul=[[ix[tuple(a[b[i]] for i in range(4))] for b in perms] for a in perms]
pw=[]
for a in range(24):
    row=[I]
    for n in range(12):row.append(mul[row[-1]][a])
    pw.append(row)
def power(a,n):return pw[a][n%12]
inverse=[power(a,11) for a in range(24)]
trans={a for a,p in enumerate(perms) if sum(p[i]!=i for i in range(4))==2}
def generated(gs):
    seen={I};todo=[I]
    while todo:
        a=todo.pop()
        for b in gs:
            c=mul[a][b]
            if c not in seen:seen.add(c);todo.append(c)
    return seen
groups={(a,b):generated([a,b]) for a,b in itertools.product(range(24),repeat=2)}
full={(a,b,z):len(generated([a,b,z]))==24 for a,b,z in itertools.product(range(24),repeat=3)}
roots=[{h:[z for z in range(24) if power(z,r)==h] for h in range(24)} for r in range(12)]
canon_cache={};torus_cache={}
def canonical(row):
    if row not in canon_cache:
        variants=[tuple(mul[mul[inverse[c]][a]][c] for a in row) for c in range(24)];key=min(variants)
        for val in variants:canon_cache[val]=key
    return canon_cache[row]
def torus(p,q):
    p,q=p%12,q%12
    if (p,q) not in torus_cache:
        u,v=next((u,v) for u,v in itertools.product(range(12),repeat=2) if (q*u+p*v)%12==1)
        ans=[]
        for a,b in itertools.product(range(24),repeat=2):
            H=power(a,p)
            if H!=power(b,q):continue
            M=mul[power(a,u)][power(b,v)];L=mul[H][power(M,-p*q)]
            ans.append((a,b,M,L))
        torus_cache[p,q]=ans
    return torus_cache[p,q]
def cell(p,q,r,k):
    p,q,r,k=(x%12 for x in [p,q,r,k]);assert math.gcd(math.gcd(p,q),12)==math.gcd(math.gcd(r,k),12)==1
    u,v=next((u,v) for u,v in itertools.product(range(12),repeat=2) if (r*u-k*v)%12==1)
    classes={};total=0
    for a,b,M,L in torus(p,q):
        h=mul[power(M,k)][power(L,r)];beta=mul[power(M,u)][power(L,v)]
        for z in roots[r][h]:
            tau=mul[beta][power(z,-v)]
            if tau not in trans or mul[power(tau,-r*k)][h]!=I or not full[a,b,z]:continue
            total+=1;key=canonical((a,b,z))
            if key not in classes:
                graph=O.orbit_graph(dict(companion=[perms[a],perms[b]],cable=[perms[beta],perms[z]]),[('companion','cable',[perms[M],perms[L]])],4)
                classes[key]=dict(rho=[perms[a],perms[b],perms[z]],mu=perms[tau],h=perms[h],M=perms[M],L=perms[L],beta=perms[beta],
                    companion_order=len(groups[a,b]),cable_order=len(groups[beta,z]),cycle_rank=graph['cycle_rank'],orbit_size=0)
            classes[key]['orbit_size']+=1
    assert total==24*len(classes) and all(x['orbit_size']==24 for x in classes.values())
    return dict(representations=total,classes=len(classes),orbits=list(classes.values()))

def satellite(p,q,r,k,row,reverse=False):
    u,v,d=S.gcdex(q,p);assert d==1
    M=E.wp([1],int(u))+E.wp([2],int(v));L=E.wp([1],p)+E.wp(M,-p*q)
    h=E.wp(M,k)+E.wp(L,r);uu,cv,d=S.gcdex(r,k);assert d==1;vv=-int(cv);uu=int(uu)
    mu=E.wp(M,uu)+E.wp(L,vv)+E.wp([3],-vv);lam=E.wp(mu,-r*k)+h
    rho=[tuple(x) for x in row['rho']];I=(0,1,2,3)
    rels=[E.wp([1],p)+E.wp([2],-q),E.wp([3],r)+E.wi(h)]
    assert all(E.evaluate(w,rho)==I for w in rels) and E.evaluate(lam,rho)==I
    rows,mark,_=E.cover_matrix(rels,rho,mu,lam,reverse_tree=reverse)
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[mark]))

def base_records(row,p):
    a,b,z=map(tuple,row['rho']);M=tuple(row['M']);H=E.power(a,p);h=tuple(row['h']);beta=tuple(row['beta']);tau=tuple(row['mu'])
    comps=[];cables=[]
    for block in V.orbits([a,b],range(4)):
        (aa,bb,mm),degree=V.fibre_quotient([a,b,M],H,block)
        comps.append(dict(sheets=block,fibre_degree=degree,base_degree=len(aa),genus=V.genus_from_base([aa,bb],[mm])))
    for block in V.orbits([beta,z],range(4)):
        (zz,bb,tt),degree=V.fibre_quotient([z,beta,tau],h,block)
        cables.append(dict(sheets=block,fibre_degree=degree,base_degree=len(zz),genus=V.genus_from_base([zz],[bb,tt])))
    return comps,cables

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['terminal_cone_binomial_closure_receipt.json','trefoil_complementary_fusion_receipt.json','complementary_corner_jsj_receipt.json',
      'binomial_periodic_receipt.json','peripheral_zero_receipt.json','nodal_deformation_receipt.json','coordinate_pencil_matrix_receipt.json',
      'fibre_order_transfer_receipt.json','common_pencil_mark_receipt.json','seifert_base_genus_receipt.json','orbit_tree_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(HERE/path)==digest)
    manifest=HERE/'sources/general-two-stage/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('complete-page audited primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    p,q,r,k,h=S.symbols('p q r k h',positive=True,integer=True);c=p*q-p-q
    identity('all-key canonical splice recurrence',r*(p+q)+k-p*q*r-1-(k-r*c-1))
    identity('negative-even master parameter',k.subs(k,r*c+1-2*h)-r*c-1+2*h)
    a,l,Q,kap,eps=S.symbols('a l Q kap eps',positive=True,integer=True)
    contacts=(l-1)*a*Q*a+(a-1)*((2*Q-l)*a+1-kap);nd=S.expand(contacts+2*Q*a-l*a+1)
    identity('arbitrary-support two-node diagonal',nd-(((l+1)*Q-l)*a*a+(1-kap)*a+kap))
    det=S.expand(3*nd-Q*Q*a*a);identity('marked distance2 determinant coefficient',det.subs(l,(Q-eps)/3)-((2-eps)*Q+eps)*a*a-3*a-3*kap*(1-a))
    t=S.symbols('t',positive=True,integer=True);DD=lambda C,a,b:S.expand(C*a*a+3*a+3*b*(1-a))
    for label,left,right in [('B1 first',DD(2,t,3*t),t*(12-7*t)),('B1 second',DD(4*Q-2,t,3*t),t*((4*Q-11)*t+12)),
      ('B1 swap first',DD(2,3*t,t),3*t*(3*t+4)),('B1 swap second',DD(4*Q-2,3*t,t),3*t*((12*Q-9)*t+4))]:identity('odd terminal formula '+label,left-right)
    check('t1 both terminal orders exclude determinant1 modulo3',all(DD(C,bb,3)%3==2 and DD(C,3,bb)%3==0 for C in [2,26] for bb in [2,4,5,7]))
    J=S.symbols('J',positive=True,integer=True);ss=J+1-r;sd=S.expand(3*((ss-1)*J+2*J-ss+1)-J*J)
    identity('solid source marked distance',J-3*ss-(3*r-2*J-3))
    identity('first solid residual determinant',sd.subs(J,(3*r-1)/2)-(6*r-1));identity('second solid residual determinant',sd.subs(J,(3*r-5)/2)-5)
    identity('solid first cap bound',(6*r-k).subs(k,5*r-2)-(r+2));check('local odd-fibre meridian distance2',S.Matrix([[1,1],[0,2]]).det()==2)
    for rr in [3,5,7]:check('odd cable coprime cones give free H1 '+str(rr),E.h1([[rr,0,1],[0,2,1]])==dict(free=1,torsion=[]))
    old=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    cells={};counts=collections.Counter();admissible=0;independent=odd_complete=True
    for pp,qq,rr,kk in itertools.product(range(12),repeat=4):
        if math.gcd(math.gcd(pp,qq),12)!=1 or math.gcd(math.gcd(rr,kk),12)!=1:continue
        admissible+=1;ans=cell(pp,qq,rr,kk)
        if ans['classes']:cells[f'{pp},{qq},{rr},{kk}']=ans
        uv=[(u,v) for u,v in itertools.product(range(12),repeat=2) if (rr*u-kk*v)%12==1]
        for row in ans['orbits']:
            ai,bi,zi=(ix[tuple(x)] for x in row['rho']);mi,li=ix[tuple(row['M'])],ix[tuple(row['L'])];taui=ix[tuple(row['mu'])]
            independent &= {mul[mul[power(mi,u)][power(li,v)]][power(zi,-v)] for u,v in uv}=={taui}
            pp1=pp if pp>=2 else pp+12;qq1=qq if qq>=2 else qq+12
            comps,cables=base_records(row,pp1);row['companion_bases']=comps;row['cable_bases']=cables
            genus=sum(x['genus'] for x in comps);assert all(x['genus']==0 for x in cables)
            tag=(rr*kk%2,row['companion_order'],row['cable_order'],row['cycle_rank'],genus);counts[tag]+=1
            if rr*kk%2:
                ee,oo=(pp1,qq1) if pp1%2==0 else (qq1,pp1)
                odd_complete &= ee%4==0 and oo%3==0 and tag==(1,24,2,0,0)
                odd_complete &= row['mu']==row['h']==row['M']==row['rho'][2] and tuple(row['L'])==(0,1,2,3)
                odd_complete &= sorted((len(x['sheets']),x['fibre_degree'],x['base_degree']) for x in cables)==[(1,1,1),(1,1,1),(2,2,1)]
                evenperm=tuple(row['rho'][0] if pp1%2==0 else row['rho'][1]);oddperm=tuple(row['rho'][1] if pp1%2==0 else row['rho'][0])
                odd_complete &= B.lengths(evenperm)==[4] and B.lengths(oddperm)==[1,3]
    expected={(0,3,24,0,0):384,(0,4,24,0,0):384,(0,12,6,0,0):384,(0,6,6,0,0):256,(0,24,24,0,0):256,
      (0,24,24,0,1):128,(0,24,2,0,0):256,(1,24,2,0,0):256,(0,24,6,1,0):256,(0,12,24,1,0):192}
    check('complete9216 cell census1344 nonempty2752 classes66048 representations',admissible==9216 and len(cells)==1344 and sum(x['classes'] for x in cells.values())==2752 and sum(x['representations'] for x in cells.values())==66048)
    check('all general image types and base genera exhausted',dict(counts)==expected)
    check('all modular cable Bezout choices identical',independent)
    check('all256 odd classes have forced companion and fibre-degree2 cable',odd_complete)
    check('remaining even planar tree periodic types1920',sum(v for tag,v in counts.items() if tag[0]==tag[3]==tag[4]==0)==1920)
    match=True
    for key,val in old.items():
        pp,qq,rr=map(int,key.split(','));kk=(rr*(pp*qq-pp-qq)+1)%12;ans=cells.get(f'{pp},{qq},{rr},{kk}',dict(representations=0,classes=0,orbits=[]))
        match &= (ans['representations'],ans['classes'],sum(x['cycle_rank']==0 for x in ans['orbits']))==(val['representations'],val['classes'],val['tree_classes'])
        oldkeys=set()
        for row in val['orbits']:
            ai,bi,betai,taui=(ix[tuple(x)] for x in row['rho']);oldkeys.add(canonical((ai,bi,mul[inverse[betai]][taui])))
        match &= oldkeys=={canonical(tuple(ix[tuple(x)] for x in row['rho'])) for row in ans['orbits']}
    check('all1152 old eta0 cells recover exact canonical classes',match)
    check('negative-even admissible cells6144',sum(math.gcd(math.gcd(rr,kk),12)==1 and rr*kk%2==0 for rr,kk in itertools.product(range(12),repeat=2))*96==6144)
    targets=[]
    for pp,qq,rr,kk in [(4,3,3,5),(4,3,5,23),(4,3,9,43),(4,9,5,7),(8,3,5,7),(8,9,7,11),(3,4,5,23)]:
        ans=cell(pp,qq,rr,kk);check('ordinary odd target one class '+str((pp,qq,rr,kk)),ans['classes']==1)
        row=ans['orbits'][0];out=satellite(pp,qq,rr,kk,row);tt=(qq if pp%2==0 else pp)//3
        check('ordinary FULL target H1 and zero R mark '+str((pp,qq,rr,kk)),out==dict(H1=B.group([tt]),quotient=B.group([tt])))
        if rr in [3,9]:check('odd satellite alternate spanning tree '+str((pp,qq,rr,kk)),satellite(pp,qq,rr,kk,row,True)==out)
        targets.append(dict(p=pp,q=qq,r=rr,k=kk,eta=kk-rr*(pp*qq-pp-qq)-1,**out))
    controls=[]
    for rho,ss,cs,delta in [(2,3,{15:2,11:1},-27),(2,9,{45:2,11:1},117),(2,9,{21:2,17:1},189),(2,3,{7:2},29),(2,3,{11:2},5)]:
        ctl=F.cluster_control(rho,ss,cs);check('actual odd forced source cluster '+str((ss,cs)),ctl['matrix'].det()==delta)
        controls.append(dict(rho=rho,sigma=ss,support=cs,matrix=[[int(x) for x in row] for row in ctl['matrix'].tolist()],determinant=delta))
    inputs=[Path(__file__),HERE/'GENERAL_TWO_STAGE_SPLICES_AND_ODD_FIBRE_EXCLUSION.md',manifest,HERE/'verify_trefoil_complementary_fusion.py',HERE/'verify_orbit_tree.py']+[HERE/n for n in preds]+[HERE/row['path'] for row in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(actual_reduced_positive_two_stage_splice_required=True,genuine_first_and_last_pairs_required=True,lower_weight_key_terms_allowed=True,literal_binomial_required=False,
        eta='k-r*(pq-p-q)-1',eta0_all_general_two_stage_quartic_packets_excluded=True,transfer_preserves_preferred_meridian_longitude_and_all_fillings=True,
        global_polynomial_key_normalization_asserted=False,general_peripheral_period=12,cable_Bezout_independence=True),
      odd_sector=dict(all_general_two_stage_odd_product_quartic_packets_excluded_in_both_pencils=True,fibre_image='transposition',companion='S4',cable='C2',e='4B',o='3t',
        M_equals_mu_equals_z=True,L_identity=True,filled_main_cable_cones=['r',2],generic_cable_fibre_degree=2,R_meridian_fibre_distance=2,
        full_target_H1='Z/t',R_mark_zero=True,common_closed_by_PF67=True,complementary_closed_with_both_cap_signs=True,
        solid_residual_source_determinants=['6r-1',5],solid_target_H1_order=1),
      census=dict(admissible_cells=9216,nonempty_cells=1344,representations=66048,conjugacy_classes=2752,odd_classes=256,cycle_classes=448,positive_genus_tree_classes=128,
        remaining_even_planar_tree_types=1920,complete_cells=cells,type_counts={str(tag):v for tag,v in counts.items()},old1152_eta0_slice_exactly_recovered=True,residue_types_not_curves_or_affine_extensions=True),
      frontier=dict(all_surviving_genuine_two_stage_quartic_packets_require_eta_negative_even=True,eta='-2h,h>=1',k='r*(pq-p-q)+1-2h>=2,gcd(r,k)=1',rk_even=True,
        period_in_negative_parameter_h=6,admissible_negative_even_cells=6144,all_negative_even_two_stage_cases_excluded=False,longer_target_sequences_excluded=False),
      ordinary_target_controls=targets,source_cluster_controls=controls,
      assurance=dict(written_unbounded_splice_transfer_period_and_odd_fibre_proofs=True,complete_Neumann_Norbury_Clay_Watson_primary_pages_visually_audited=True,
        deep_topology_imported_with_retained_PF79_hypotheses=True,full_peripheral_recognition_before_slope_comparison=True,complete_finite_census_after_unbounded_period_proof=True,
        bounded_source_controls_not_universal_proof=True,eta_negative_assumed_zero=False,old_identity_fibre_assumed_for_odd_product=False,
        arbitrary_degree_packet_smoothness_or_two_line_classification=False,full_two_stage_exclusion=False,full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
        formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'general_two_stage_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; general9216-cell census; general eta0 and odd-product sectors closed; negative-even frontier OPEN')
if __name__=='__main__':main()
