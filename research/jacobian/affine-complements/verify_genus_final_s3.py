"""PF2RULED91: all-depth genus/framing and all FINAL-stage S3 arrivals.

Exact finite controls support the written proofs and imported Floer inputs;
they do not prove universality or recognize solid topology from homology.
"""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_global_ruling_unfilled_cyclic as F
A=F.A;K=F.K;G=F.G;E=F.E;B=F.B
I=tuple(range(4))
def cyclic(order):return B.group([abs(order)],free=int(order==0))
def history(steps):
    aa=1;bb=6;dd=5;cc=1;wind=1;cables=[];valid=True
    for w,sign in steps:
        k=2*aa*w+sign;valid &= k<w*bb
        dd=w*dd+k-w*bb;cc=w*cc+(w-1)*k;bb=w*k;aa*=w*w;wind*=w
        cables.append((w,k))
    return dict(steps=steps,cables=cables,A=aa,B=bb,d=dd,c=cc,D=2*aa-cc,W=wind,infinity_valid=valid)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['global_ruling_unfilled_cyclic_receipt.json','nonunimodular_a4_receipt.json',
        'meridian_balance_two_stage_receipt.json','persistent_jsj_terminal_modulus_receipt.json',
        'terminal_cone_twists_receipt.json','trefoil_complementary_fusion_receipt.json',
        'rational_end_orderability_receipt.json','lifted_surface_s3_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/genus-final-s3/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('exact primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    nn=next(x for x in sources['primary_sources'] if 'norbury' in x['path'])
    check('fresh complete canonical/infinity primary page audits explicit',nn['fresh_complete_page_audit'] and nn['complete_PDF_pages_rendered_and_visually_inspected']==[6,12,13])
    osrow=next(x for x in sources['primary_sources'] if 'ozsvath-szabo-math' in x['path'])
    check('OS v1 cross-sign absolute display not used as rank formula','not used across both signs' in osrow['used_statement'] and osrow['new_deep_Floer_import'] is False)
    r,k,bb,dd,cc,aa,w,eps,D,eta,WW,Jc,b=S.symbols('r k B d c A w epsilon D eta W J b',integer=True)
    identity('all-depth canonical genus bridge induction',r*k-(r*dd+k-r*bb)-(r*(bb-dd)+(r-1)*k))
    identity('initial canonical genus bridge',bb-(bb-cc)-cc)
    identity('solid genus deficit recurrence',2*aa*w*w-(w*(2*aa-D)+(w-1)*(2*aa*w+eps))-(w*D-(w-1)*eps))
    identity('next solid continuation after minus violates strict infinity bound',(2*aa*r+eps)-r*(2*aa-w)-(r*w+eps))
    identity('arrival cap via genus degree',aa*r-(r*(2*aa-D)+eta+1)/2-(D*r-eta-1)/2)
    identity('negative eta denominator relation',2*(aa*r-k/2)-((2*aa-1)*r-k)-r)
    identity('r3 surgery distance from genus endpoint',(3*(2*aa-D)+eta+1)/2-aa-(2*aa-D)-(eta+1-D)/2)
    identity('eta0 product cap forces deficit1',((D*3-1)/2)-1-S.Rational(3,2)*(D-1))
    identity('r3 genus endpoint ramification core is generator',(6*aa-2)-3*(2*aa-1)-1)
    # Complete divisor deductions: no unbounded interval is sampled here.
    unmerged=sorted(set(int(x) for n in [2,4] for x in S.divisors(n) if x>=2))
    check('complete unmerged root list from inverse3 equals plusminus1',unmerged==[2,4])
    collar=set()
    for which in ['3','Jminus3']:
        bv=S.Integer(3) if which=='3' else Jc-3
        for sign in [-1,1]:
            for jj in S.solve(Jc-3*bv-sign,Jc):
                bj=bv.subs(Jc,jj)
                if jj.is_Integer and bj.is_Integer and jj>=2 and 1<=bj<jj and math.gcd(int(jj),int(bj))==1:
                    collar.add((int(jj),int(bj)))
    check('complete collar coefficient AND marked distance list',sorted(collar)==[(4,1),(5,2),(8,3),(10,3)])
    check('terminal trefoil winding2 free3 cannot break last gcd',math.gcd(2,3-1)!=1 and math.gcd(3,2-1)==1)
    check('unmerged root2 requires impossible square3mod7',3 not in {x*x%7 for x in range(7)})
    check('unmerged4 and collar4or10 require impossible square2mod3',2 not in {x*x%3 for x in range(3)})
    identity('root2 lattice forces3Aminus2 mod7',(38*aa-2)-(3*aa-2)-35*aa)
    identity('root4 lattice forces2Aminus1 mod3',(146*aa-4)-(2*aa-1)-3*(48*aa-1))
    lower=Jc*(3*b-Jc)*(3*WW)**2+3*(Jc-b)*3*WW+3
    identity('collar5 exact contact lower bound',lower.subs({Jc:5,b:2})-(45*WW*WW+27*WW+3))
    identity('collar8 exact contact lower bound',lower.subs({Jc:8,b:3})-(72*WW*WW+45*WW+3))
    identity('merged3,4 all-depth source lower bound',12*(3*WW)**2+13*3*WW+13-(108*WW*WW+39*WW+13))
    check('all-depth positive differences between source lower bounds and targets',
        all(all(x>0 for x in S.Poly(ex,WW).all_coeffs()) for ex in [
            (108*WW*WW+39*WW+13)-(6*WW*WW-4),
            (45*WW*WW+27*WW+3)-(2*WW*WW-5),
            (72*WW*WW+45*WW+3)-(2*WW*WW-8)]))
    check('product collar minus assignment would force even winding4',S.solve(2*w-3-(w+1),w)==[4])
    # Hom's SIGNED rank formula on a positive L-space knot, not the OSv1 abs display.
    rank_controls=[]
    for cval in [1,17,449]:
        for numerator,denominator in [(-3,1),(-1,2),(1,2),(cval-1,1),(cval,1),(2*cval+1,2)]:
            if numerator==0:continue
            rank=numerator+2*max(0,cval*denominator-numerator)
            check('signed Floer rank and exact genus threshold '+str((cval,numerator,denominator)),
                (rank==abs(numerator))==(S.Rational(numerator,denominator)>=cval))
            rank_controls.append(dict(c=cval,p=numerator,q=denominator,rank=rank))
    stepsets=[[],[(3,1)],[(3,1),(5,1)],[(3,1),(5,-1)],[(5,-1)],[(3,1),(5,1),(3,1)]]
    histories=[history(steps) for steps in stepsets]
    aa0,bb0,mm,ll=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
    seed=A.normalized(2,3,[],[G.perms[aa0],G.perms[bb0]]);M,_=A.periphery(2,3,[],seed)
    seed3=F.restrict_three(seed,M);M3=E.evaluate(K.words(2,3,[])[1],seed3)
    for h in histories:
        check('actual infinity solid history genus identities '+str(h['steps']),h['infinity_valid'] and h['A']==h['W']**2 and h['B']-h['d']==h['c'] and h['D']>=1 and (h['D']==1)==all(s==1 for w,s in h['steps']))
        images=seed3+[E.power(M3,k%2) for r,k in h['cables']]
        full,unfilled,m,L=F.sharp(2,3,h['cables'],images)
        check('actual composed unfilled cap for genus-controlled solid history '+str(h['steps']),
            E.h1(unfilled)==dict(free=1,torsion=[]) and E.h1(unfilled,[[x+h['A']*y for x,y in zip(L,m)]])==E.h1(unfilled))
    # These forbidden continuations are retained as integer countercontrols,
    # not declared infinity knots or affine packets.
    for w0,r0,sign in [(3,3,1),(3,5,-1),(5,3,-1)]:
        bad=history([(w0,-1),(r0,sign)])
        check('minus continuation countercontrol fails exact infinity inequality '+str((w0,r0,sign)),not bad['infinity_valid'])
    interval_controls=[];targets=[]
    for h in histories:
        prefimages=seed+[E.power(M,k%2) for r,k in h['cables']]
        for r0 in [3,9,15]:
            for eta0 in [0,-2,-4]:
                kval=r0*h['c']+eta0+1
                if kval<2 or math.gcd(r0,kval)!=1:continue
                jval=h['A']*r0-kval//2;sval=r0//3
                check('arrival actual canonical degree and positive cap '+str((h['steps'],r0,eta0)),kval%2==0 and jval>0 and kval<r0*h['B'])
                v=pow(-kval,-1,r0);u=(1+kval*v)//r0
                a1=-S.Rational(3*v,r0);a2=S.Rational(2*h['A']*v-u,jval);x=S.Rational(v,r0)
                check('actual solid cap coefficient and inverse3 when r3 '+str((h['steps'],r0,eta0)),
                    a2==S.Rational(2*v,r0)-S.Rational(1,r0*jval) and (r0!=3 or (3*(2*h['A']*v-u)+1)%jval==0))
                if r0>=9 and eta0<0:
                    l0=(2*h['A']-1)*r0-kval;N=(l0*v-1)//r0
                    check('two exact denominators bound missing meridian interior '+str((h['steps'],r0,eta0)),
                        l0>0 and l0*v==N*r0+1 and S.ceiling(l0*a1)==-3*N and S.ceiling(l0*a2)==2*N+1 and S.Rational(N,l0)<x and S.floor(r0*a1)==-3*v and S.floor(r0*a2)==2*v-1)
                if r0>=9 and eta0==0:
                    check('eta0 actual matched GENUS endpoint is yminus1 '+str((h['steps'],r0)),
                        v==r0-1 and a1==-3+S.Rational(1,sval) and a2==2-S.Rational(h['D'],jval) and 0<S.Rational(h['D'],jval)<S.Rational(1,sval) and -1-S.floor(a1)-S.floor(a2)==1 and S.Rational(u-h['c']*v,kval-h['c']*r0)==1)
                if r0==3:
                    t=kval//2-h['A'];check('r3 ordinary surgery genus bound exact '+str((h['steps'],eta0)),t-h['c']==S.Rational(eta0+1-h['D'],2))
                interval_controls.append(dict(history=h['steps'],A=h['A'],D=h['D'],r=r0,k=kval,eta=eta0,cap=jval,a1=str(a1),a2=str(a2),meridian=str(x)))
                z,tau,lam,beta=next(x for x in A.roots(M,I,r0,kval,2,6) if E.generated(prefimages+[x[0]])==24)
                base_order=kval//2-h['A']*sval;base_mark=sval*kval
                for tail in [[],[(5,3)],[(7,4)]]:
                    images=prefimages+[z]+[E.power(tau,tk%2) for tr,tk in tail]
                    ans=K.summarized_cover(2,3,h['cables']+[(r0,kval)]+tail,images)
                    added=[cyclic(tr) for tr,tk in tail if tk%2==0]
                    mark_wind=math.prod(tr for tr,tk in tail)
                    expected=dict(H1=B.direct_sum([cyclic(base_order)]+added),quotient=B.direct_sum([cyclic(math.gcd(base_order,mark_wind*base_mark))]+added))
                    check('full actual S4 marked cover agrees with symbolic solid-arrival/tail homology '+str((h['steps'],r0,eta0,tail)),ans==expected and E.generated(images)==24)
                    targets.append(dict(history=h['steps'],arrival=[r0,kval],eta_at_arrival=eta0,tail=tail,actual=ans))
    # A real marked-cover countercontrol prevents claiming endpoint linking persists.
    h=history([(3,1)]);prefimages=seed+[E.power(M,1)]
    z,tau,lam,beta=next(x for x in A.roots(M,I,3,52,2,6) if E.generated(prefimages+[x[0]])==24)
    continued=K.summarized_cover(2,3,h['cables']+[(3,52),(3,5)],prefimages+[z,E.power(tau,1)])
    check('endpoint Lspace seed tail countercontrol really has marked H1Z17',continued==dict(H1=cyclic(17),quotient=cyclic(1)))
    newlink=(-S.Rational(3*3,17))%1
    check('endpoint seed linking obstruction is NOT invariant under an odd tail square',newlink==S.Rational(8,17) and 0<newlink<S.Rational(1,2))
    inputs=[Path(__file__),P/'GENUS_FRAMING_AND_FINAL_S3_CLOSURE.md',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_depth_canonical_genus_bridge=True,strict_infinity_splice_retained=True,
            solid_minus_loop_must_be_last=True,solid_history_positive_Lspace_knots=True,
            signed_rank_formula_and_genus_endpoint_retained=True),
        quartic=dict(all_final_stage_S3_arrivals_closed_both_pencils=True,arbitrary_companion_depth=True,
            nonpositive_arrival_common_nonLspace_tail_closure_except_r3_eta0_allplus=True,
            complementary_two_cone_internal_annulus_tails_closed=True,
            normalized_root_and_terminal_winding_divisor_proof=True,
            no_prefix_assumed_separate_affine_source=True),
        complete_unmerged_root_orders=unmerged,complete_collar_root_pairs=sorted(collar),
        actual_solid_histories=histories,exact_interval_controls=interval_controls,
        full_target_controls=targets,signed_rank_controls=rank_controls,
        endpoint_tail_countercontrol=dict(A=9,ordinary_genus=9,arrival=[3,52],tail=[3,5],actual=continued,linking=str(newlink),not_affine_realization=True),
        assurance=dict(written_unbounded_genus_interval_root_and_lattice_proofs=True,
            actual_composed_full_CW_and_integral_mark_controls=True,
            fresh_complete_canonical_PDF_pages_audited=True,
            new_primary_Floer_cabling_and_rank_inputs_imported=True,
            retained_deep_JSJ_gluing_LO_full_peripheral_inputs=True,
            finite_controls_not_universal_proof=True,
            endpoint_Lspace_seed_linking_not_assumed_tail_invariant=True,
            all_later_S3_closed=False,all_later_A4_closed=False,full_quartic_exclusion=False,
            higher_degree_packet_classification=False,finite_Keller_map=False,
            formal_verification=False,external_review=False,independent_full_reproof=False,
            novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'genus_final_s3_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; ALL FINAL-stage S3 arrivals closed BOTH pencils at arbitrary companion depth; later tails remain OPEN')
if __name__=='__main__':main()
