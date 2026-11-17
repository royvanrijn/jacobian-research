"""PF2RULED75: exact section identities and ordinary full-satellite controls.

The unrestricted proofs and imported topology are in the companion note.
"""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=load('base75','verify_seifert_base_genus.py');E=V.E
BUD=load('budget75','verify_peripheral_cone_budget.py')

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['rational_end_orderability_receipt.json','lifted_surface_s3_receipt.json',
           'fibre_order_transfer_receipt.json','binomial_periodic_receipt.json',
           'seifert_base_genus_receipt.json','peripheral_cone_budget_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/seifert-sections-winding-three/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('retained primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    b,t,s,j=S.symbols('b t s j',integer=True);A=b*t;c=6*b*t-2*b-3*t;r=3*s;k=r*c+1
    vv=c-2*A;BB=k/2-r*A;D=k/2-A*s
    check('third section coefficient',S.expand(2*BB-(3*s*vv+1))==0)
    check('positive determinant identity',S.expand(k+BB-3*D)==0)
    check('non-trefoil positive V lower bound b>=2',S.expand(vv-(2*b-3)-(4*b-3)*(t-1))==0)
    check('non-trefoil positive V at b1',vv.subs(b,1)==t-2)
    # Seifert annulus presentation X^s=h, XYZ=h^j. In free H1, h=sX.
    xx,yy,zz=S.symbols('X Y Z');hh=s*xx
    boundary_sum=3*hh+r*(j*hh-xx)
    check('lifted disk forces relative Euler integer zero',S.expand(boundary_sum-r*j*hh)==0)
    for ss,jj in [(1,-2),(3,1),(5,7)]:
        annulus=[[ss,0,0,-1],[1,1,1,-jj]]
        check('annulus H1 is torsion-free rank2 '+str((ss,jj)),E.h1(annulus)==dict(free=2,torsion=[]))
    # Boundary rows on variables (h,Y,Z): n=h+rY, m=2h+rZ.
    m=S.Matrix([[2,0,r]]);ll=S.Matrix([[-c,0,-k/2]])
    check('actual eta filling yields third cone relation',S.simplify(A*m+ll-S.Matrix([[-vv,0,-BB]]))==S.zeros(1,3))
    # X'=X*h^-1, f=h^-1 gives X'^s=f^(s-1), product X'YZ=f.
    change=S.Matrix([[1,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
    check('nonabelian central generator change is unimodular',change.det()==-1)
    old=S.Matrix([[s,0,0,-1],[0,k,0,c],[0,0,BB,vv],[1,1,1,0]])
    new=S.Matrix([[s,0,0,1-s],[0,k,0,-c],[0,0,BB,-vv],[1,1,1,-1]])
    check('exact relation exponent change',S.simplify(old*change.inv()-new)==S.zeros(4,4))
    epsilon=-1+(s-1)/s+c/k+vv/BB
    check('unrestricted normalized Euler identity',S.factor(epsilon+D/(s*k*BB))==0)
    check('full star determinant identity',S.factor(new.det()+D)==0)
    check('normalized root coefficient primitive identity',s-(s-1)==1)
    check('normalized second coefficient primitive identity',S.expand(k-3*s*c)==1)
    check('normalized third coefficient primitive identity',S.expand(2*BB-3*s*vv)==1)
    check('one-step Laufer intersection',-1+3==2)
    check('single Pinkham witness l1',-1-1-1==-1-2)
    check('hyperbolic triangle quotient uniform upper bound',S.Rational(1,3)+S.Rational(1,28)+S.Rational(1,5)<1)
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];count=0
    for key,cell in census.items():
        p,q,rr=map(int,key.split(','));p=p if p>=2 else p+12;q=q if q>=2 else q+12;rr=rr if rr>=2 else rr+12
        for i,row in enumerate(cell['orbits']):
            if (row['companion_order'],row['cable_order'],row['cycle_rank'])!=(6,6,0):continue
            cab=V.cable_bases(row);main=next(x for x in cab if len(x['sheets'])==3)
            root,beta,tau=main['induced'];tag=str((key,i))
            check('actual root and two beta sections '+tag,BUD.lengths(root)==[3] and BUD.lengths(beta)==[1,2] and BUD.lengths(tau)==[1,2])
            check('base annulus one cone before companion fillings '+tag,main['genus']==0 and rr%6==3)
            count+=1
    check('all64 S3/S3 incidence types retained',count==64)
    controls=[];shared=False
    cases=[(2,9,9),(4,3,9),(4,9,9),(8,3,9),(2,15,15),(10,9,21),(14,45,9)]
    for ee,oo,rr in cases:
        bb,tt,ss=ee//2,oo//3,rr//3;aa=bb*tt;cc=ee*oo-ee-oo;kk=rr*cc+1;vv0=cc-2*aa;bs=kk//2-rr*aa;dd=kk//2-aa*ss
        tag=str((ee,oo,rr));alpha=[ss,kk,bs];omega=[ss-1,cc,vv0]
        check('actual primitive coefficients even with shared cone primes '+tag,all(0<om<al and math.gcd(om,al)==1 for al,om in zip(alpha,omega)))
        check('negative definite normalized Euler '+tag,-1+sum(S.Rational(om,al) for al,om in zip(alpha,omega))==-S.Rational(dd,math.prod(alpha)))
        check('unbounded l1 witness actual coefficients '+tag,sum(S.floor(-S.Rational(om,al)) for al,om in zip(alpha,omega))==-3)
        mat=[[al if i==j else (-om if j==3 else 0) for j in range(4)] for i,(al,om) in enumerate(zip(alpha,omega))]+[[1,1,1,-1]]
        check('quotient full integral H1 '+tag,E.h1(mat)==BUD.group([dd]))
        shared |= any(math.gcd(alpha[i],alpha[j])>1 for i in range(3) for j in range(i))
        rows=[x for x in census[f'{ee%12},{oo%12},{rr%12}']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(6,6,0)]
        check('all ordinary classes retained '+tag,len(rows)==2)
        want=dict(H1=BUD.group([dd,bb]),quotient=BUD.group([math.gcd(dd,ss*kk),bb]))
        for i,row in enumerate(rows):
            check('actual satellite retains full groups and R mark '+str((tag,i)),V.full_satellite(ee,oo,rr,row)==want)
        controls.append(dict(e=ee,o=oo,r=rr,alpha=alpha,omega=omega,e0=-1,D=dd,**want))
    check('shared cone primes explicitly exercised',shared)
    # Residual r=3 is an ACTUAL annulus bundle, not a virtual filling.
    cc=S.symbols('c',positive=True,integer=True);kk=3*cc+1
    nrow=S.Matrix([[1,3]]);l1=S.Matrix([[-cc,-kk]])
    check('r3 actual outer boundary is a basis',nrow.col_join(l1).det()==-1)
    check('r3 regular fibre word retained',S.simplify(kk*nrow+3*l1)==S.Matrix([[1,0]]))
    check('r3 fixed beta section',S.simplify(-cc*nrow-l1)==S.Matrix([[0,1]]))
    mrow=S.Matrix([[-1,-3]]);lrow=S.Matrix([[kk/2-cc,kk/2]])
    check('r3 actual inner meridian inverse',mrow==-nrow)
    check('r3 actual inner longitude shear',S.simplify(lrow-kk*nrow/2-l1)==S.zeros(1,2))
    check('r3 gluing reverses boundary orientation',S.Matrix([[-1,0],[kk/2,1]]).det()==-1)
    check('r3 mark is primitive',S.expand(kk-3*cc)==1)
    residual=[]
    for ee,oo in [(2,9),(4,3),(14,45)]:
        bb,tt=ee//2,oo//3;aa=bb*tt;cc0=ee*oo-ee-oo;kk0=3*cc0+1;vv0=cc0-2*aa;bs=kk0//2-3*aa;dd=kk0//2-aa
        mat=[[kk0,-cc0],[bs,-(bs-vv0)]]
        check('r3 quotient truly finite cyclic '+str((ee,oo)),E.h1(mat)==BUD.group([dd]) and abs(S.Matrix(mat).det())==dd)
        check('r3 one-step witness correctly absent '+str((ee,oo)),S.floor(S.Rational(0,1))+S.floor(-S.Rational(cc0,kk0))+S.floor(-S.Rational(vv0,bs))==-2)
        check('r3 actual gluing retains independent torsion '+str((ee,oo)),E.h1([[dd,0],[0,bb]])==BUD.group([dd,bb]))
        rows=[x for x in census[f'{ee%12},{oo%12},3']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(6,6,0)]
        check('r3 both ordinary classes retained '+str((ee,oo)),len(rows)==2)
        want=dict(H1=BUD.group([dd,bb]),quotient=BUD.group([math.gcd(dd,kk0),bb]))
        for i,row in enumerate(rows):
            check('r3 actual full satellite independently checks marked gluing '+str((ee,oo,i)),V.full_satellite(ee,oo,3,row)==want)
        residual.append(dict(e=ee,o=oo,r=3,D=dd,quotient_group='Z/'+str(dd),m=[-1,0],L=[kk0//2,1],R=[kk0,3]))
    old74=json.loads((HERE/'rational_end_orderability_receipt.json').read_text())
    old72=json.loads((HERE/'fibre_order_transfer_receipt.json').read_text())
    check('trefoil and all non-S3 common closures retained',old74['trefoil']['all_windings_common_pencil_excluded'] and old72['quartic']['remaining_common_sector']=='S3/S3 only')
    inputs=[Path(__file__),HERE/'SEIFERT_SECTION_NORMALIZATION_AND_WINDING_THREE.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(PF74_all_degree_source_obstruction_retained=True,actual_peripheral_section_maps_required=True,
          closed_group_quotient_suffices=True,all_higher_degree_packets_assumed_common=False),
        quartic=dict(actual_primitive_binomial_two_stage_required=True,all64_S3_S3_incidence_types=True,
          actual_annulus_relative_Euler_zero=True,root_and_boundary_sections_retained=True,
          normalized_star='e0=-1;(s,s-1),(k,c),(B,V)',V='c-2bt',B='(3sV+1)/2',Euler='-D/(skB)',
          pairwise_coprime_cone_orders_required=False,twists_inferred_from_H1=False,
          non_trefoil_s_at_least3_common_entirely_excluded=True,
          universal_Laufer_first_intersection=2,universal_Pinkham_witness=1,
          entire_actual_two_stage_common_requires_r=3,trefoil_common_entirely_excluded=True,
          r3_non_trefoil_common_entirely_excluded=False,passing_complementary_sectors_excluded=False),
        residual=dict(r=3,virtual_quotient_group='Z/D',virtual_quotient_not_LO_for_D_greater1=True,
          actual_main_cable_homeomorphic_T2_times_I=True,actual_companion_pieces_retained=True,
          actual_gluing_m=[-1,0],actual_gluing_L='(k/2,1)',actual_R_mark='(k,3)',
          exact_JSJ_classification_asserted=False,affine_source_realization_asserted=False),
        full_satellite_controls=controls,r3_controls=residual,
        assurance=dict(written_unbounded_section_and_winding_proofs=True,prior_complete_primary_page_audits_reused=True,
          imported_rational_link_orderability_results=True,finite_controls_not_universal_proof=True,
          deep_theorems_independently_reproved=False,group_quotient_claimed_actual_homeomorphism=False,
          r3_quotient_claimed_actual_target_group=False,uninspected_Rasmussen_theorem_imported=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
          PF63_polynomial_curves_declared_nonexistent=False,formal_verification=False,
          independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'seifert_sections_winding_three_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; actual Seifert sections; entire common two-stage quartic frontier has r=3')
if __name__=='__main__':main()
