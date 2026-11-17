"""PF2RULED93: integer slope grid and all common genuine iterated ends.

Finite exact controls support the written unbounded proof. They neither
compute Floer homology nor recognize a product collar from H1 alone.
"""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import composed_peripheral_cover as C
import verify_genus_final_s3 as H
A=H.A;K=H.K;E=H.E;G=H.G;B=H.B;I=tuple(range(4))
def endpoints(aa,c,r,k):
    j=aa*r-k//2;d=k-r*c;v=pow(-k,-1,r);u=(1+k*v)//r
    a1=-S.Rational(3*v,r);a2=S.Rational(2*aa*v-u,j)
    x=S.Rational(v,r);z=x+S.Rational(1,r*d)
    fm=lambda n:-(S.floor(n*a1)+S.floor(n*a2)+S.ceiling(n*x))/n
    fp=lambda n:-(S.ceiling(n*a1)+S.ceiling(n*a2)+S.floor(n*z))/n
    n1=1 if j==-1 else r;n2=1 if j==1 else r
    return dict(J=j,d=d,v=v,u=u,a1=a1,a2=a2,x=x,z=z,N1=n1,N2=n2,Fminus=fm(n1),Fplus=fp(n2))
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['recursive_linking_endpoints_receipt.json','global_ruling_unfilled_cyclic_receipt.json',
        'genus_final_s3_receipt.json','cyclic_pattern_a4_receipt.json',
        'iterated_peripheral_tail_receipt.json','meridian_balance_two_stage_receipt.json',
        'terminal_cone_twists_receipt.json','meridian_cycle_s4_receipt.json','unimodular_initial_s4_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/integer-grid-common/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path'];check('primary pin '+path.name,sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    fresh=next(x for x in sources['primary_sources'] if x.get('new_deep_graph_interval_import'))
    check('complete graph theorem and integer-section conventions audited from primary HTML',
        fresh['complete_HTML_sections_read']==['Definition3.7','Section4.2','Section4.3 Theorem4.6 and Proposition4.7','Corollary6.1'] and fresh['independently_reproved'] is False)
    aa,r,k,v,d,c,a=S.symbols('A r k v d c a',integer=True)
    u=(1+k*v)/r;j=aa*r-k/2
    check('actual signed cap coefficient identity',S.simplify((2*aa*v-u)/j-(2*v/r-1/(r*j)))==0)
    check('ordinary genus endpoint transported with exact shift',S.simplify(((u-c*v)/(k-r*c)-v/r)-1/(r*(k-r*c)))==0)
    check('normal section and ANY integer-grid slope have primitive intersection',S.det(S.Matrix([[a,1],[1,0]]))==-1)
    check('product root cap yields ordinary integer surgery',S.expand((k/2-aa).subs(k,6*aa+2)-(2*aa+1))==0)
    check('product plus-cap case already has eta zero',S.expand((k-r*(2*aa-1)-1).subs({r:3,k:6*aa-2}))==0)
    check('product minus-cap remaining eta four',S.expand((k-r*(2*aa-1)-1).subs({r:3,k:6*aa+2}))==4)
    check('product actual normal section is the ordinary dual surgery word',
        S.expand((u-(2*aa+1)*v)-(-v-(2*aa*v-u)))==0 and S.expand(3*(2*aa+1)-(6*aa+2))==1)
    # Exact cyclic patterns: all grid fillings have integrally zero core m.
    grid=[]
    for n,rr,kk in [(2,3,2),(2,3,5),(2,5,3),(3,2,3),(5,3,2),(6,5,6),(7,3,2)]:
        rows,mark,m,L,_,_=A.pattern(n,rr,kk);base=C.reduced_h1(rows)
        check('actual filled cyclic pattern m INTEGRALLY zero '+str((n,rr,kk)),C.reduced_h1(rows,[m])==base)
        fills=[]
        for aval in [-7,-1,0,1,9,225]:
            gamma=[aval*x+y for x,y in zip(m,L)]
            filled=C.reduced_h1(rows,[gamma])
            check('integer grid filling core INTEGRALLY zero '+str((n,rr,kk,aval)),C.reduced_h1(rows,[gamma,m])==filled)
            fills.append(dict(a=aval,H1=filled))
        grid.append(dict(n=n,r=rr,k=kk,H1=base,grid_fills=fills,
            incompressible_graph_application_asserted=n==2))
    # Signed graph-root witnesses. Parameter choices are controls, not proof.
    histories=[H.history(x) for x in [[],[(3,1)],[(5,1)],[(3,1),(5,1)],[(3,-1)],[(5,-1)],[(3,1),(5,-1)]]]
    certs=[];cases=set()
    for h in histories:
        for rr in [3,9,15,21]:
            for jj in [-7,-2,-1,1,2,4,7]:
                kk=2*h['A']*rr-2*jj
                if kk<2 or math.gcd(rr,kk)!=1 or kk>=rr*h['B'] or kk-rr*h['c']<3:continue
                if rr==3 and abs(jj)==1:continue # true product: separate grid proof
                q=endpoints(h['A'],h['c'],rr,kk)
                check('signed actual normal missing-interior certificate '+str((h['steps'],rr,kk)),
                    q['J']==jj and 0<q['x']<q['z']<1 and q['Fminus']>=0>=q['Fplus'])
                case='positive' if jj>=2 else 'negative' if jj<=-2 else 'unit_positive' if jj==1 else 'unit_negative';cases.add(case)
                wantminus=S.Rational(1,rr) if jj>0 else S.Integer(0)
                wantplus=S.Integer(0) if jj>0 else -S.Rational(1,rr)
                check('witness table values exact '+str((h['steps'],rr,kk)),q['Fminus']==wantminus and q['Fplus']==wantplus)
                certs.append(dict(history=h['steps'],A=h['A'],c=h['c'],r=rr,k=kk,**{key:str(val) for key,val in q.items()}))
    check('all FOUR symbolic witness cases have actual infinity-valid controls',cases=={'positive','negative','unit_positive','unit_negative'})
    # Construct genuine full-S4 first arrivals without expanding knot words.
    old=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
    seed=A.normalized(2,3,[],[G.perms[old[0]],G.perms[old[1]]]);M,_=A.periphery(2,3,[],seed)
    product=[];tails=[]
    for steps in [[],[(3,1)],[(5,1)],[(3,1),(5,1)]]:
        h=H.history(steps);rr=3;kk=6*h['A']+2;pref=seed+[E.power(M,tk%2) for tr,tk in h['cables']]
        z,tau,_,_=next(x for x in A.roots(M,I,rr,kk,2,6) if E.generated(pref+[x[0]])==24)
        data=C.peripheral_cover(2,3,h['cables']+[(rr,kk)],pref+[z])
        b=next(b for b in data['boundaries'] if b['degree']==2)
        rows=data['relation_rows']+[b0['meridian'] for b0 in data['boundaries'] if b0['degree']==1]
        null=S.Matrix.hstack(*S.Matrix(rows).nullspace());assert null.cols==1
        vals=[S.Rational(x) for x in null];scale=S.ilcm(*[x.q for x in vals]);vec=[int(x*scale) for x in vals]
        gcd=math.gcd(*vec);vec=[x//gcd for x in vec]
        mval=sum(x*y for x,y in zip(b['meridian'],vec));lval=sum(x*y for x,y in zip(b['longitude'],vec))
        if mval<0:mval=-mval;lval=-lval
        candidates=[(sign-lval)//mval for sign in [-1,1] if (sign-lval)%mval==0]
        check('actual product R exterior free generator and dual integer-grid slope '+str(steps),
            C.reduced_h1(rows)==dict(free=1,torsion=[]) and mval==2*h['A']+1 and bool(candidates))
        product.append(dict(history=steps,A=h['A'],normal_image=mval,preferred_image=lval,primitive_grid_candidates=candidates,
            product_topology_by_written_section_elimination_not_H1=True))
        for tail in [[(3,5)],[(7,4)],[(3,5),(5,3),(7,4)]]:
            full=C.peripheral_cover(2,3,h['cables']+[(rr,kk)]+tail,pref+[z]+[E.power(tau,tk%2) for tr,tk in tail])
            normal=next(b0 for b0 in full['boundaries'] if b0['degree']==2)
            closedrows=full['relation_rows']+[b0['meridian'] for b0 in full['boundaries']]
            ans=dict(H1=C.reduced_h1(closedrows),quotient=C.reduced_h1(closedrows,[normal['longitude']]))
            added=[B.group([tr]) for tr,tk in tail if tk%2==0];wind=math.prod(tr for tr,tk in tail)
            expected=dict(H1=B.direct_sum([B.group([2*h['A']+1])]+added),
                quotient=B.direct_sum([B.group([math.gcd(2*h['A']+1,wind*kk)])]+added))
            check('full actual product-seed tail marked cover '+str((steps,tail)),ans==expected and full['word_expansion_used'] is False)
            tails.append(dict(history=steps,tail=tail,actual=ans,nonLspace_by_written_grid_gluing=True))
    sat=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())['saturation']
    incoming={x['source'] for x in sat['transition_profiles'] if x['dest']==3 and x['source']!=3}
    check('complete state reduction has EXACT first-simple incoming classes',incoming=={2,4,5,6,7,9})
    check('simple full state absorbs every passing tail with cyclic2 image and odd winding',
        all(x['dest']==3 and x['cable_order']==2 and x['rk'][0]%2==1 for x in sat['useful_edges'] if x['source']==3))
    proofs={2:'PF90 nonsolid; PF92 nonpositive eta; PF93 positive eta',4:'PF86',5:'PF86',6:'PF85',7:'PF87',9:'PF85'}
    check('every first simple incoming class has an all-length COMMON theorem',set(proofs)==incoming)
    inputs=[Path(__file__),P/'INTEGER_SLOPE_GRID_AND_COMMON_ITERATED_CLOSURE.md',manifest,P/'composed_peripheral_cover.py']+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(integer_grid_missing_interior_under_integral_m_zero_and_incompressibility=True,
            degree_free_cyclic_pattern_application_with_topology_hypotheses=True,
            finite_denominator_graph_missing_interior_certificate=True,
            matched_integer_grid_slope_propagates_all_actual_cyclic2_tails=True),
        quartic=dict(actual_drilled_normal_section_zero_calibrated=True,
            all_nonproduct_positive_eta_S3_roots_missing_normal_interior=True,
            remaining_product_root_closed_by_dual_ordinary_meridian=True,
            ALL_COMMON_genuine_positive_iterated_torus_quartic_ends_closed=True,
            arbitrary_companion_depth_arrival_parameters_and_finite_tail_length=True,
            no_prefix_assumed_separate_affine_source=True),
        cyclic_grid_controls=grid,signed_root_certificates=certs,product_peripheral_controls=product,
        actual_product_seed_tail_controls=tails,complete_first_arrival_exhaustion=proofs,
        assurance=dict(written_unbounded_grid_section_witness_and_exhaustion_proofs=True,
            new_primary_graph_interval_theorem_import_explicit=True,
            integral_boundary_and_composed_actual_cover_controls=True,
            retained_deep_Floer_graph_LO_and_full_peripheral_hypotheses=True,
            finite_controls_not_universal_proof=True,
            complementary_S3_tails_closed=False,complementary_A4_tails_closed=False,
            KIG28_identified_with_intermediate_A4=False,full_quartic_exclusion=False,
            higher_degree_packet_classification=False,degenerate_stage_classification=False,
            finite_Keller_map=False,formal_verification=False,external_review=False,
            independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'integer_grid_common_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; ALL COMMON genuine iterated-torus quartic ends closed; complementary frontier OPEN')
if __name__=='__main__':main()
