"""Exact finite controls for PF2RULED54; no knot or curve passing them is declared Keller."""
from pathlib import Path
import hashlib,itertools,json,math,os,platform
import sympy as S

HERE=Path(__file__).resolve().parent;checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):return tuple(a.index(i) for i in range(len(a)))
def group(gens,n):
    seen={tuple(range(n))};todo=list(seen)
    while todo:
        a=todo.pop()
        for b in gens:
            c=compose(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return seen
I=(1,2)
GL=[(a,b) for a in [1,2,3] for b in [1,2,3] if a!=b]
def act(A,v):return (A[0] if v&1 else 0)^(A[1] if v&2 else 0)
def mul(A,B):return (act(A,B[0]),act(A,B[1]))
def inv(A):return next(B for B in GL if mul(A,B)==I)
def transform(A,z):
    return ((z[0] if A[0]&1 else 0)^(z[1] if A[1]&1 else 0),(z[0] if A[0]&2 else 0)^(z[1] if A[1]&2 else 0))
def word_linear(word,mats):
    A=I;z=(0,0)
    for letter in word:
        j=abs(letter)-1;B=mats[j];v=(1<<(2*j),1<<(2*j+1))
        if letter<0:B=inv(B);v=transform(B,v)
        vv=transform(A,v);z=(z[0]^vv[0],z[1]^vv[1]);A=mul(A,B)
    return A,z
def rank2(rows):
    pivots={}
    for row in rows:
        while row:
            bit=row.bit_length()-1
            if bit in pivots:row^=pivots[bit]
            else:pivots[bit]=row;break
    return len(pivots)
def affine(A,v):return tuple(act(A,i)^v for i in range(4))
def project(p):
    v=p[0]
    return (p[1]^v,p[2]^v)
def perm_word(word,gens):
    out=tuple(range(len(gens[0])))
    for j in word:out=compose(out,gens[j-1] if j>0 else inverse(gens[-j-1]))
    return out
def cocycle_control(mats,relations,meridians):
    rows=[]
    for rel in relations:
        A,z=word_linear(rel,mats);assert A==I;rows.extend(z)
    for mu in meridians:
        A,z=word_linear(mu,mats);assert A!=I and mul(A,A)==I
        Az=transform(A,z);rows.extend([z[0]^Az[0],z[1]^Az[1]])
    count=2*len(mats);solutions=[z for z in range(1<<count) if all((r&z).bit_count()%2==0 for r in rows)]
    cob=[]
    for v in [1,2]:
        z=0
        for j,A in enumerate(mats):z|=(v^act(A,v))<<(2*j)
        cob.append(z)
    cobspace={0,cob[0],cob[1],cob[0]^cob[1]}
    assert cobspace<=set(solutions)
    lifts=[]
    for z in solutions:
        gens=[affine(A,(z>>(2*j))&3) for j,A in enumerate(mats)]
        for rel in relations:assert perm_word(rel,gens)==tuple(range(4))
        for mu in meridians:
            p=perm_word(mu,gens);assert compose(p,p)==tuple(range(4)) and sum(p[i]!=i for i in range(4))==2
        order=len(group(gens,4));assert order==(6 if z in cobspace else 24)
        lifts.append(dict(cocycle=z,image_order=order,coboundary=z in cobspace))
    return dict(cocycle_dimension=count-rank2(rows),coboundary_dimension=rank2(cob),parabolic_H1_dimension=count-rank2(rows)-rank2(cob),solutions=lifts)
def main():
    preds=['unramified_orevkov_receipt.json','torus_infinity_receipt.json','trace_dual_pencil_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for path,digest in r['source_sha256'].items():check('predecessor seal '+name+' '+path,sha(HERE/path)==digest)
    manifest=HERE/'sources/quartic-colouring/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('primary source '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    perms=list(itertools.permutations(range(4)));id4=tuple(range(4))
    check('affine F2 plane realizes all S4 permutations',{affine(A,v) for A in GL for v in range(4)}==set(perms))
    check('linear quotient has all six S3 elements',len(GL)==6 and len({tuple(act(A,v) for v in [1,2,3]) for A in GL})==6)
    kernel=[p for p in perms if project(p)==I]
    check('quotient kernel is exactly the Klein four translations',len(kernel)==4 and all(p==affine(I,p[0]) for p in kernel))
    check('S4 linear projection respects every product',all(project(compose(a,b))==mul(project(a),project(b)) for a in perms for b in perms))
    trans=[p for p in perms if sum(p[i]!=i for i in range(4))==2]
    check('every single transposition remains a nontrivial linear transposition',all(project(p)!=I and mul(project(p),project(p))==I for p in trans))
    for A in GL:
        if A!=I and mul(A,A)==I:
            lifts=[affine(A,v) for v in range(4)]
            check('order2 parabolic lift condition for linear transposition '+str(A),all((v^act(A,v)==0)==(compose(lifts[v],lifts[v])==id4) for v in range(4)))
            check('exactly two transposition lifts and two four-cycle lifts '+str(A),sum(p in trans for p in lifts)==2 and sum(compose(p,p)!=id4 for p in lifts)==2)
    # Identify the three reflection colours in S3 directly.
    colours=[(1,0,2),(2,1,0),(0,2,1)]
    check('Fox relation is exactly transposition conjugation',all(compose(compose(colours[b],colours[a]),colours[b])==colours[(2*b-a)%3] for a in range(3) for b in range(3)))
    fox=S.Matrix([[2,-1,-1],[-1,2,-1],[-1,-1,2]])
    check('trefoil has determinant3 and a nonconstant Fox colouring',fox[:2,:2].det()==3 and all(v%3==0 for v in fox*S.Matrix([0,1,2])))
    # The local double-value V4 becomes one C2, killed on the sign surface.
    for a in trans:
        for b in trans:
            if a!=b and compose(a,b)==compose(b,a):
                H=group([a,b],4);even=[p for p in H if p in kernel]
                check('disjoint fold local sign kernel is killed in S3 '+str((a,b)),project(a)==project(b) and len(H)==4 and len(even)==2 and all(project(p)==I for p in even))
    # Complete linear lift systems on two classical knot-group presentations.
    trefoil=cocycle_control([(2,1),(1,3)],[[1,2,1,-2,-1,-2]],[[1],[2]])
    check('trefoil passes Fox3 but has zero parabolic S4 lift space',trefoil['parabolic_H1_dimension']==0 and all(row['image_order']==6 for row in trefoil['solutions']))
    aa=(1,2,0,3);bb=(1,2,3,0)
    torus34=cocycle_control([project(aa),project(bb)],[[1,1,1,-2,-2,-2,-2]],[[1,-2]])
    check('T34 has a nonzero parabolic space and actual transposition S4 lifts',torus34['parabolic_H1_dimension']>0 and any(row['image_order']==24 for row in torus34['solutions']))
    check('explicit T34 S4 representation has the exact meridian',len(group([aa,bb],4))==24 and perm_word([1,-2],[aa,bb]) in trans)
    # Closed cable determinant recurrence against exact Alexander polynomials.
    t=S.symbols('t')
    def torus(w,p):
        p=abs(p);return S.cancel((t**(w*p)-1)*(t-1)/((t**w-1)*(t**p-1)))
    def cable_det(D,w,p):return abs(p) if w%2==0 else w*D if p%2==0 else D
    cable_rows=[]
    for stages in [[(3,4)],[(3,5)],[(3,4),(2,5)],[(3,4),(2,3)],[(3,5),(3,2)],[(3,4),(3,5)],[(3,4),(2,-5)],[(2,3),(3,2),(2,7)]]:
        poly=S.Integer(1);D=1
        for w,p in stages:
            assert math.gcd(w,p)==1
            poly=S.cancel(torus(w,p)*poly.subs(t,t**w));D=cable_det(D,w,p)
        actual=abs(int(poly.subs(t,-1)))
        check('all-stage cable determinant exact '+str(stages),poly.is_polynomial(t) and poly.subs(t,1)==1 and actual==D)
        cable_rows.append(dict(stages=stages,determinant=D,quartic_det_gate=D%3==0,polynomial_realization_asserted=False))
    check('final even winding can erase companion three-colourability',cable_rows[2]['determinant']==5 and cable_rows[2]['quartic_det_gate'] is False)
    check('odd winding can add a factor three',cable_rows[4]['determinant']==3 and cable_rows[4]['quartic_det_gate'])
    # Oka's smooth C* curve and a genuine S3 cubic, not a Keller source.
    x,u,s=S.symbols('x u s');f=(x*u-2)**2-u**3;W=s**3-3*u*s-2*(x*u-2)
    check('Oka cubic has exact sign discriminant',S.expand(S.discriminant(W,s)+108*f)==0)
    check('Oka cubic is primitive linear in x',S.degree(W,x)==1 and S.gcd(-2*u,s**3-3*u*s+4)==1)
    check('Oka branch curve is smooth',list(S.groebner([f,S.diff(f,x),S.diff(f,u)],x,u))==[1])
    xx=t+2/t**2;uu=t**2
    check('Oka curve has the exact Laurent Cstar parametrization',S.cancel(f.subs({x:xx,u:uu},simultaneous=True))==0)
    check('Oka normalization parameter has a regular inverse away from u0',S.cancel(((x*u-2)/u).subs({x:xx,u:uu},simultaneous=True)-t)==0 and f.subs(u,0)==4)
    check('Oka smooth sign cover control is not an A1-normalized branch curve',S.degree(S.fraction(S.cancel(xx))[1],t)>0 and uu==t*t)
    inputs=[Path(__file__),HERE/'QUARTIC_COLOURING_AND_ETALE_TRIPLE_GATE.md',manifest]+[HERE/name for name in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      scope=dict(actual_surviving_quartic_packet=True,all_source_pole_orders=True,all_target_image_degrees=True,all_common_or_complementary_pencils=True,other_unramified_prime_retained=True,actual_global_S4_required=True,determinant_gate_asserted_for_degrees_at_least_five=False),
      infinity=dict(primary_epimorphism_without_general_position_at_infinity=True,actual_S4_survives=True,quotient='S4/V4=S3',nonconstant_Fox3_required=True,necessary='3 divides abs(Delta_K(-1))',prime_to_three_knots_excluded=True,cabling_formula='D_new=abs(p) if winding even; winding*D if winding odd and p even; D if both odd',all_finite_cabling_lengths=True,cable_controls=cable_rows),
      lifting=dict(module='V=F2^2 with GL2(F2)=S3',parabolic_condition='(I+phi(mu))*z(mu)=0',quotient_by_coboundaries=True,surjective_S4_lift_iff_nonzero_parabolic_H1=True,H1_S3_V_zero=True,linear_once_phi_fixed=True,knot_infinity_solution_implies_global_solution=False,trefoil=trefoil,torus34=torus34),
      sign_surface=dict(model='zeta^2=G_R',connected_etale_cyclic_three_cover_required=True,etale_at_finite_singularities=True,finite_double_branch_local_image_C2_after_S3_quotient=True,nonzero_antiinvariant_H1_et_mu3_required=True,unit_term_in_Kummer_sequence_retained=True,Pic3_nonzero_requires_units_all_cubes=True,cubic_subfield_of_quartic_asserted=False,simply_connected_from_smoothness=False),
      negative_controls=dict(three_colourability_not_S4_lift=True,topological_lift_not_Keller=True,even_winding_erases_inner_colourability=True,smooth_affine_sign_surface_can_have_etale_triple=True,Oka_branch_is_Cstar_not_A1=True),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,imported_primary_theorems_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,finite_Keller_pair_supplied=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'quartic_colouring_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; quartic determinant exclusion, linear S4 lifts and etale sign triple gate')
if __name__=='__main__':main()
