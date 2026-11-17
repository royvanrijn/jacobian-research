"""PF2RULED73: exact identities, incidence exhaustion and independent controls.

Universal surface, linking and group arguments are in the proof note.
Finite controls do not establish them or certify affine realizability.
"""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=load('base73','verify_seifert_base_genus.py');E=V.E
B=load('budget73','verify_peripheral_cone_budget.py')
C=load('cluster73','verify_common_pencil_mark.py')

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['fibre_order_transfer_receipt.json','fold_schur_s4_receipt.json',
           'peripheral_cone_budget_receipt.json','seifert_base_genus_receipt.json',
           'common_pencil_mark_receipt.json','oriented_pattern_linking_receipt.json',
           'binomial_periodic_receipt.json','marked_end_lattice_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/lifted-surface-s3/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('retained primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];count=0
    for key,cell in census.items():
        p,q,r=map(int,key.split(','));p=p if p>=2 else p+12;q=q if q>=2 else q+12;r=r if r>=2 else r+12
        for i,row in enumerate(cell['orbits']):
            if (row['companion_order'],row['cable_order'],row['cycle_rank'])!=(6,6,0):continue
            tag=str((key,i));a,b,beta,tau=map(tuple,row['rho']);M=E.power(beta,r);H=E.power(a,p)
            comp=V.companion_bases(row,p,q,r);cab=V.cable_bases(row)
            check('S3 companion and cable sizes '+tag,sorted(len(x['sheets']) for x in comp)==[1,3] and sorted(len(x['sheets']) for x in cab)==[1,3])
            cc=next(x for x in comp if len(x['sheets'])==3);pp=next(x for x in cab if len(x['sheets'])==3)
            intersection=set(cc['sheets'])&set(pp['sheets'])
            check('distinct fixed sheets and central transposition '+tag,len(intersection)==2 and {j for j in range(4) if M[j]!=j}==intersection)
            edges=[len(set(x['sheets'])&set(y['sheets'])) for x in comp for y in cab if set(x['sheets'])&set(y['sheets'])]
            check('complete path with edge degrees1,1,2 '+tag,sorted(edges)==[1,1,2] and H==tuple(range(4)) and all(x['genus']==0 for x in comp+cab))
            ee,oo=(p,q) if p%2==0 else (q,p)
            check('actual exponents and root orders '+tag,ee%2==0 and oo%6==3 and r%6==3 and B.lengths(E.mul(E.inv(tau),beta))==[1,3])
            count+=1
    check('all64 S3/S3 classes exhausted',count==64)
    b,t,s,u,v=S.symbols('b t s u v',integer=True);A=b*t;c=6*b*t-2*b-3*t;k=3*s*c+1;D=k/2-A*s
    check('positive determinant formula',S.expand(2*D-(s*(16*b*t-6*b-9*t)+1))==0)
    check('C lower bound',S.expand((16*b*t-6*b-9*t)-(7*t-6)-(b-1)*(16*t-6))==0)
    check('lifted longitude saturation identity',S.expand(2*S.Symbol('L')+2*A*S.Symbol('m')+S.Symbol('L')+A*S.Symbol('m')-3*(S.Symbol('L')+A*S.Symbol('m')))==0)
    check('cable meridian surface and fibre balances',-s+s==0 and S.expand((k/2)*(-s)+3*s*(k/2)-s*k)==0)
    bez=-s*u+(k/2)*v-1;Q=A*v-u
    check('oriented boundary determinant is negative Bezout',S.expand(u*s-v*k/2+bez+1)==0)
    check('eta decomposition first coordinate',S.expand(D*u+Q*k/2-A-A*bez)==0)
    check('eta decomposition second coordinate',S.expand(D*v+Q*s-1-bez)==0)
    check('linking numerator inverse identity',S.expand(Q*s-(1-D*v)-bez)==0)
    beta=S.symbols('beta',integer=True)
    check('cubic negative-star rescaling congruence',S.expand(9*t*t*beta*(1+2*t*beta)+t-t*(3*t*beta+1)*(6*t*beta+1))==0)
    controls=[]
    cases=[(2,3,3),(2,3,15),(4,3,3),(14,3,3),(2,9,3),(2,3,63),(4,9,15),(8,3,9)]
    for ee,oo,rr in cases:
        bb,tt,ss=ee//2,oo//3,rr//3;aa=bb*tt;kk=rr*(ee*oo-ee-oo)+1;dd=kk//2-aa*ss
        tag=str((ee,oo,rr));gg=math.gcd(dd,ss*kk)
        check('coprime actual parameters '+tag,math.gcd(ee,oo)==1 and math.gcd(bb,3)==1 and dd>=1 and math.gcd(dd,ss)==1 and math.gcd(ss,kk)==1)
        check('reduced companion torsion '+tag,math.prod([bb,2*bb,tt])//math.lcm(bb,2*bb,tt)==bb)
        check('closed companion Euler determinant '+tag,abs(-S.Rational(1,2*bb*tt)*bb*2*bb*tt)==bb)
        check('three torsion saturation allowed '+tag,math.gcd(3,bb)==1)
        # Actual gluing with variables (g,m,z), L=-A*m and z of order b.
        mat=[[-ss,-1,0],[kk//2,aa,0],[0,0,bb]]
        want=dict(H1=B.group([dd,bb]),quotient=B.group([gg,bb]))
        check('unrestricted integral gluing '+tag,dict(H1=E.h1(mat),quotient=E.h1(mat,[[ss*kk,0,0]]))==want)
        check('new marked cyclic quotient gcd gate '+tag,math.gcd(gg,bb)==math.gcd(dd,bb))
        uu,vv,gcd=S.gcdex(-ss,kk//2);qq=aa*int(vv)-int(uu)
        check('actual oriented surgery coordinates '+tag,gcd==1 and qq*ss%dd==1%dd and S.Matrix([[int(uu),kk//2],[int(vv),ss]]).det()==-1)
        if bb>1:
            bet=(-pow(3*tt,-1,bb))%bb
            check('actual cubic oriented rescaling '+tag,(9*tt*tt*bet*(1+2*tt*bet)+tt)%bb==0)
        rows=[x for x in census[f'{ee%12},{oo%12},{rr%12}']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(6,6,0)]
        check('both actual classes retained '+tag,len(rows)==2)
        for i,row in enumerate(rows):
            actual=V.full_satellite(ee,oo,rr,row)
            check('full ordinary satellite comparison '+str((ee,oo,rr,i)),actual==want)
            if i==0:check('reverse spanning tree control '+tag,V.full_satellite(ee,oo,rr,row,True)==actual)
            controls.append(dict(e=ee,o=oo,r=rr,D=dd,b=bb,mark=ss*kk,class_index=i,**actual))
    for bb in [2,4,5,7]:
        # Complex t=1 model is the negative A_(b-1) resolution chain.
        chain=S.diag(*([-2]*(bb-1)))
        for i in range(bb-2):chain[i,i+1]=chain[i+1,i]=1
        check('A-chain outward negative-inverse sign calibration '+str(bb),S.simplify(-chain.inv()[0,0]+S.Rational(1,bb))==1)
    # Trefoil small cases and the even-D full-marked obstruction.
    for dd in range(1,9):
        residue=4%dd
        check('trefoil small determinant common excluded '+str(dd),not(0<2*residue<dd and residue%2==0))
    dd,nr,cv=S.symbols('D nr cv',integer=True)
    check('even determinant contradiction identity',S.expand(4*(nr-cv*cv)-(4*nr-(2*cv)**2))==0)
    check('even D16 self-linking countercontrol',S.Matrix([[5,2],[2,4]]).det()==16 and all((x*x-5)%16 for x in range(16) if (2*x+2)%16==0))
    J=S.symbols('J',integer=True,nonnegative=True);j=2*J+1
    nr=3*j*j+7*j+5;cc=3*j+4;delta=3*j*j+4*j+4;x=(delta-cc)/2
    check('unbounded odd family Gram determinant',S.expand(4*nr-cc*cc-delta)==0)
    check('unbounded odd family determinant mod8',S.Poly(S.expand(delta-3),J).all_coeffs()==[12,20,8])
    # J(J+1) even: delta-3 = 8*(J+1)+12*J*(J+1).
    check('unbounded parity certificate',S.expand(delta-3-8*(J+1)-12*J*(J+1))==0)
    check('both source relations in target cyclic group',S.expand(2*nr+cc*x-delta*(cc+1)/2)==0 and S.expand(2*cc+4*x-2*delta)==0)
    check('marked linking mixed difference is integral',S.expand(2*x+cc-delta)==0)
    check('marked linking second diagonal difference is integral',S.expand(x*x-nr-delta*(delta-2*cc-1)/4)==0)
    check('cyclic source chain determinant formula',S.expand((3*(j+2)-2)*(j+1)-3*j-delta)==0)
    check('target third cone determinant',S.expand(3*(2*delta-1)-(6*delta-2)/2-(3*delta-2))==0)
    check('triangle quotient hyperbolic bound',3*S.Rational(1,21)<1)
    z,w,T,W=S.symbols('z w T W');clusters=[]
    for jj in [1,3,5]:
        R=C.local_cluster(2*z**(-2*(jj+1))+z**(-jj)+z**(1-jj)*T,z**(jj+1),0,z,T)
        U=C.local_cluster(2*w**-2+W,w,0,w,W);out=C.assemble([R,U]);mat=out['full'];labels=out['labels']
        want=S.Matrix([[3*jj*jj+7*jj+5,3*jj+4],[3*jj+4,4]]);de=3*jj*jj+4*jj+4
        tag=str(jj)
        check('whole cluster Gram '+tag,out['source']==-want and out['multiplicities']*out['multiplicities'].T==want)
        check('whole leaves and canonical volume '+tag,out['degrees']==[1,1] and [out['orders'][i] for i in out['index']]==[1,0])
        rows=R['rows'];chain=[labels[rr['key']] for rr in rows[:3]]+[labels[rr['key']] for rr in reversed(rows[4:-1])]
        check('exact shared and crossing counts '+tag,len(rows)==jj+5 and rows[3]['key']==U['terminal'] and [rr['multiplicity'] for rr in rows]==[jj+1]*3+[1]*(jj+2))
        check('exact source T chain weights '+tag,[mat[i,i] for i in chain]==[-2,-2,-(jj+2)]+[-2]*jj)
        hx,hy=labels['Hx'],labels['Hy'];retained=chain+[hx,hy]
        check('all infinity vertices retained '+tag,set(retained)==set(range(mat.rows))-set(out['index']))
        expected={(min(chain[i],chain[i+1]),max(chain[i],chain[i+1])) for i in range(len(chain)-1)}|{tuple(sorted([chain[1],hx])),tuple(sorted([hx,hy]))}
        actual={(i,h) for i in retained for h in retained if i<h and mat[i,h]!=0}
        check('exact source T plumbing tree '+tag,actual==expected and mat[hx,hx]==-2 and mat[hy,hy]==0)
        # Nonabelian elimination is proved in Section8; this checks its chain relations.
        chainmat=mat.extract(chain,chain)
        powers=[1,-int(chainmat[0,0])]
        for i in range(1,len(chain)-1):powers.append(-int(chainmat[i,i])*powers[-1]-powers[-2])
        last=-int(chainmat[-1,-1])*powers[-1]-powers[-2]
        check('chain group elimination terminal relation '+tag,last==de and abs(chainmat.det())==de)
        si,ki,third=2*de-1,6*de-2,3*de-2
        check('actual infinite target quotient orders '+tag,de>=11 and S.Rational(1,si)+S.Rational(1,ki)+S.Rational(1,third)<1)
        clusters.append(dict(j=jj,Gram=[[int(a) for a in row] for row in want.tolist()],delta=de,source_pi1='Z/'+str(de),target_triangle_quotient=[si,ki,third]))
    inputs=[Path(__file__),HERE/'LIFTED_SURFACE_BALANCES_AND_S3_MARKED_LINKING.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_cover_degrees=True,oriented_meridian_disk_balance=True,
          oriented_Seifert_surface_balance=True,individual_longitude_vanishing_asserted=False,
          torsion_linking_orthogonal_replacement=True,actual_primitive_boundary_maps_required=True,
          boundary_torsion_coordinates_zero_required=True,interior_mark_torsion_retained=True,
          peripheral_compatible_piece_epimorphism_induces_glued_group_quotient=True,
          knot_exterior_longitude_filling_group_quotient=True,group_equality_for_replacement_asserted=False),
        quartic=dict(actual_primitive_binomial_two_stage_required=True,all64_S3_S3_classes=True,
          full_H1='Z/D plus Z/b',oriented_linking='negative inverse(s)/D orthogonal negative t/b',
          ramification_mark='(sk*g,0)',marked_quotient='Z/gcd(D,sk) plus Z/b',
          new_necessary_gate='gcd(D,b)=1',trefoil_common_necessary_r_mod48=15,
          trefoil_common_necessary_r_min=63,unbounded_full_marked_linking_matching_cluster_family=True,
          same_explicit_family_excluded_by_full_end_group=True,
          actual_target_group_surjects_virtual_filling_group=True,
          virtual_filling_claimed_homeomorphic_to_actual_target=False,
          all_sources_for_those_target_keys_excluded=False,entire_S3_S3_common_sector_excluded=False),
        full_satellite_controls=controls,cluster_and_group_controls=clusters,
        assurance=dict(written_unrestricted_surface_linking_and_group_proofs=True,
          prior_complete_primary_page_audits_reused=True,finite_controls_not_universal_proof=True,
          source_affineness_proved=False,finite_Keller_map=False,full_quartic_exclusion=False,
          higher_degree_exclusion=False,PF63_polynomial_curves_declared_nonexistent=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,
          novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'lifted_surface_s3_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; full S3/S3 marked linking; unbounded paired-family end-group mismatch')
if __name__=='__main__':main()
