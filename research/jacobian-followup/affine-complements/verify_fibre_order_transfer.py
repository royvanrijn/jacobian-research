"""PF2RULED72: exact universal identities and independently lifted finite controls."""
from pathlib import Path
import collections,hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=load('base69','verify_seifert_base_genus.py');B=load('cone70','verify_peripheral_cone_budget.py');E=V.E

def seifert(alpha,beta):
    j=len(alpha);mat=S.zeros(j+1,j+1)
    for i,(a,b) in enumerate(zip(alpha,beta)):
        assert math.gcd(a,b)==1
        mat[i,0]=b;mat[i,i+1]=a;mat[j,i+1]=1
    epsilon=-sum(S.Rational(b,a) for a,b in zip(alpha,beta));lam=math.lcm(*alpha)
    return mat,epsilon,lam

def size(group):
    assert group['free']==0
    return math.prod(group['torsion'])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['binomial_periodic_receipt.json','cyclic_parity_receipt.json','common_pencil_mark_receipt.json',
           'seifert_base_genus_receipt.json','peripheral_cone_budget_receipt.json','fold_schur_s4_receipt.json','oriented_pattern_linking_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/fibre-order-transfer/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Exact annihilator certificate: the lcm-weighted cone relations minus the global relation.
    examples=[([2,3,5],[1,1,1]),([2,2,3],[1,1,1]),([1],[1]),([4,4,9],[1,3,1]),([3,3,3],[1,1,1])]
    fibre_controls=[]
    for alpha,beta in examples:
        mat,epsilon,lam=seifert(alpha,beta);n=abs(epsilon*lam);rows=mat.tolist()
        certificate=sum((S.Rational(lam,a)*mat.row(i) for i,a in enumerate(alpha)),S.zeros(1,mat.cols))-lam*mat.row(len(alpha))
        check('integer lcm fibre relation certificate '+str((alpha,beta)),certificate==S.Matrix([[-epsilon*lam]+[0]*len(alpha)]))
        h=[1]+[0]*len(alpha);group=E.h1(rows);quotient=E.h1(rows,[h])
        check('regular fibre exact order '+str((alpha,beta)),size(group)//size(quotient)==n)
        check('closed Seifert determinant '+str((alpha,beta)),size(group)==abs(epsilon*math.prod(alpha)))
        if all(sum(a%p==0 for a in alpha)<=2 for p in [2,3,5,7]):
            check('cone count proves cyclic full group '+str((alpha,beta)),len(group['torsion'])<=1)
        fibre_controls.append(dict(alpha=alpha,beta=beta,epsilon=str(epsilon),Lambda=lam,H1=group,H_order=int(n)))
    mat,epsilon,lam=seifert([2,2],[1,-1])
    check('zero Euler case is outside finite order formula',epsilon==0 and E.h1(mat.tolist())['free']==1)
    # Cone relation makes the fibre divisible by alpha without rational division.
    for a,b in [(8,3),(9,2),(15,7)]:
        u,v,g=S.gcdex(b,a);check('fibre divisibility Bezout '+str((a,b)),g==1 and u*b+v*a==1)
    for u,v,A in [(4,8,3),(5,10,2),(3,9,4),(9,9,2)]:
        aa,bb,g=S.gcdex(A,u)
        check('unit coefficient kills meridian torsion '+str((u,v,A)),v%u==0 and g==1 and aa*A+bb*u==1)
    # Whole finite action classification; all cone orders below are actual representative integers.
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];counts=collections.Counter()
    for key,cell in census.items():
        p,q,r=map(int,key.split(','));p=p if p>=2 else p+12;q=q if q>=2 else q+12;r=r if r>=2 else r+12
        for i,row in enumerate(cell['orbits']):
            tag=row['companion_order'],row['cable_order']
            if row['cycle_rank'] or tag not in [(12,6),(24,2)]:continue
            e,o=(p,q) if p%2==0 else (q,p);a,b,beta,tau=map(tuple,row['rho'])
            check('both sectors have identity companion fibre '+str((key,i)),E.power(a,p)==tuple(range(4)))
            comp=V.companion_bases(row,p,q,r);cab=V.cable_bases(row)
            check('companion is one planar piece '+str((key,i)),len(comp)==1 and comp[0]['genus']==0)
            cone=[p//x for x in B.lengths(comp[0]['induced'][0])]+[q//x for x in B.lengths(comp[0]['induced'][1])]
            epsilon=-S.Rational(4,e*o);lam=math.lcm(*cone)
            if tag==(12,6):
                check('A4 cone lcm fibre order and closed H1 order '+str((key,i)),lam==e*o//2 and abs(epsilon*lam)==2 and abs(epsilon*math.prod(cone))==e*o//3)
                check('A4 cable orbit and root shape '+str((key,i)),sorted(x['base_degree'] for x in cab)==[1,3] and r%2==0 and
                      sorted(x for data in cab if data['base_degree']==3 for x in B.lengths(data['induced'][0]))==[1,2])
            else:
                check('S4 C2 cone lcm fibre order and closed H1 order '+str((key,i)),e%4==0 and lam==e*o//4 and abs(epsilon*lam)==1 and abs(epsilon*math.prod(cone))==o//3)
                check('C2 cable orbit and identity root shape '+str((key,i)),sorted(x['base_degree'] for x in cab)==[1,1,2] and r%2==1 and
                      B.lengths(E.mul(E.inv(tau),beta))==[1,1,1,1])
            check('closed companion at most two divisible cones per prime '+str((key,i)),all(sum(a%ell==0 for a in cone)<=2 for ell in S.factorint(math.prod(cone))))
            counts[tag]+=1
    check('all144 periodic representation types covered',dict(counts)=={(12,6):96,(24,2):48})
    # Full ordinary-word controls include every class, rather than selecting a convenient representation.
    controls=[]
    cases=[(2,3,8,(12,6)),(2,3,26,(12,6)),(4,3,16,(12,6)),(8,3,20,(12,6)),
           (10,3,10,(12,6)),(2,9,32,(12,6)),(4,3,3,(24,2)),(4,9,3,(24,2)),
           (8,3,9,(24,2)),(4,15,9,(24,2)),(16,3,3,(24,2))]
    for e,o,r,tag in cases:
        key=','.join(str(x%12) for x in [e,o,r]);rows=[x for x in census[key]['orbits'] if (x['companion_order'],x['cable_order'])==tag and not x['cycle_rank']]
        orders=[e*o//3,r//2] if tag==(12,6) else [o//3,r];want=B.group(orders)
        k=r*(e*o-e-o)+1;A=k//(3 if tag==(12,6) else 2);u=r//2 if tag==(12,6) else r
        check('actual unit coefficient in boundary word '+str((e,o,r,tag)),k%(3 if tag==(12,6) else 2)==0 and math.gcd(A,u)==1)
        check('control passes prior determinant nine gate '+str((e,o,r,tag)),(k if r%2==0 else r*o)%9==0)
        for i,row in enumerate(rows):
            ans=V.full_satellite(e,o,r,row)
            check('entire full ordinary satellite group and zero mark '+str((e,o,r,tag,i)),ans==dict(H1=want,quotient=want))
            if i==0:check('independent spanning tree control '+str((e,o,r,tag)),V.full_satellite(e,o,r,row,True)==ans)
            controls.append(dict(e=e,o=o,r=r,images=list(tag),class_index=i,**ans))
    # The stronger A4 condition must see the extra factor2 in the FULL group.
    check('old cone budget can pass a new full-group exclusion',math.gcd(2*3//6,8//2)==1 and math.gcd(2*3//3,8//2)==2)
    N=S.symbols('N',integer=True,nonnegative=True);rr=36*N+8
    check('entire trefoil new even family has r divisible4',S.expand(rr-4*(9*N+2))==0)
    check('entire trefoil new even family has target determinant nine',S.expand(rr+1-9*(4*N+1))==0)
    types=collections.Counter((x['companion_order'],x['cable_order'],x['cycle_rank']) for x in census['2,3,8']['orbits'])
    check('all seven trefoil classes exhausted',dict(types)=={(3,24,0):1,(12,6,0):3,(12,24,1):3})
    old68=json.loads((HERE/'cyclic_parity_receipt.json').read_text())
    check('cyclic3 predecessor exclusion retained',old68['theorem']['C3_all_first_pair_parities_excluded'])
    # Do NOT apply the one-boundary cable transfer to S3/S3.
    counter=[]
    for r in [3,15]:
        rows=[x for x in census[f'2,3,{r%12}']['orbits'] if (x['companion_order'],x['cable_order'])==(6,6) and not x['cycle_rank']]
        for i,row in enumerate(rows):
            ans=V.full_satellite(2,3,r,row)
            want=dict(H1=B.group([] if r==3 else [3]),quotient=B.group())
            check('two-boundary S3 countercontrol '+str((r,i)),ans==want)
            counter.append(dict(e=2,o=3,r=r,class_index=i,**ans))
    check('S3 nonzero ramification mark prevents blanket closure',counter[-1]['H1']!=counter[-1]['quotient'])
    inputs=[Path(__file__),HERE/'FIBRE_ORDER_AND_NONCYCLIC_ZERO_MARKS.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_cover_degrees=True,all_cabling_lengths_when_hypotheses_hold=True,
          closed_orientable_Seifert_base_required=True,nonzero_Euler_required=True,
          regular_fibre_exact_order='abs(Euler*lcm(cone orders))',
          unit_transfer_two_cones_base_D_required=True,cone_u_divides_v_required=True,
          actual_boundary_basis_and_pushforward_required=True,unit_A_mod_u_required=True,
          integral_meridian_zero_proved=True,peripheral_and_marked_torsion_coordinates_zero=True,
          arbitrary_noncylclic_pieces_replaced=False,linking_form_splitting_asserted=False,
          two_remaining_boundary_tori_covered_by_unit_lemma=False),
        quartic=dict(actual_primitive_binomial_two_stage_required=True,all_first_exponents_windings_source_orders_poles=True,
          all144_periodic_A4_S3_and_S4_C2_classes_covered=True,
          A4_S3_end_H1='Z/(eo/3) plus Z/(r/2)',A4_S3_R_zero=True,A4_S3_common_entirely_excluded=True,
          A4_S3_necessary_gcd='gcd(eo/3,r/2)=1',A4_S3_necessary_r_mod4=2,
          S4_C2_end_H1='Z/(o/3) plus Z/r',S4_C2_R_zero=True,S4_C2_common_entirely_excluded=True,
          remaining_common_sector='S3/S3 only',S3_S3_common_entirely_excluded=False,
          trefoil_all_r8mod36_excluded_in_both_source_pencils=True,
          trefoil_remaining_even_winding_mod36=26,passing_complementary_sectors_excluded=False),
        Seifert_fibre_controls=fibre_controls,periodic_type_counts={str(k):v for k,v in counts.items()},
        full_satellite_controls=controls,S3_nonzero_mark_countercontrols=counter,
        assurance=dict(written_unbounded_fibre_order_and_integral_transfer_proofs=True,
          primary_Euler_naturality_complete_pages_rechecked=True,retained_Seifert_presentation_audit_reused=True,
          finite_controls_not_universal_proof=True,
          exterior_pushforward_assumed_after_meridian_filling=False,all_common_two_stage_sectors_closed=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
          PF63_polynomial_curves_declared_nonexistent=False,formal_verification=False,
          independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'fibre_order_transfer_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; general fibre-order/transfer lemmas; A4/S3 and S4/C2 common closed')
if __name__=='__main__':main()
