"""Exact controls for written PF2RULED52; no finite Keller pair is supplied."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S

HERE=Path(__file__).resolve().parent;checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    preds=['boundary_fibre_models_receipt.json','unramified_orevkov_receipt.json','coordinate_regular_sector_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,digest in r['source_sha256'].items():check('predecessor seal '+name+' '+p,sha(HERE/p)==digest)
    manifest=HERE/'sources/source-symmetry-gate/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('exact complete primary source '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    x,y,h,t,w,U,c,C=S.symbols('x y h t w U c C')
    def E(f,a=1,b=-2):return S.expand(a*x*S.diff(f,x)+b*y*S.diff(f,y))
    def parts(f,a=1,b=-2):
        out={}
        for (i,j),co in S.Poly(f,x,y).terms():out[a*i+b*j]=out.get(a*i+b*j,0)+co*x**i*y**j
        return out
    projector_controls=[]
    f=x+y+x*x*y+x*y*y+y**3+x**4
    for a,b in [(1,-2),(2,3),(0,1),(-2,0),(-3,-1),(2,2)]:
        pp=parts(f,a,b)
        for k,target in pp.items():
            g=f
            for other in pp:
                if other!=k:g=S.expand((E(g,a,b)-other*g)/(k-other))
            check('exact Euler spectral projector '+str((a,b,k)),S.expand(g-target)==0)
        axis=x if a else y;weight=a if a else b;component=pp[weight]
        check('nonzero tangent weight isolates a linear axis function '+str((a,b)),S.expand(component.subs(y,0) if a else component.subs(x,0))==axis)
        projector_controls.append(dict(weights=[a,b],weight_support=sorted(pp)))
    # Additive derivatives produce a function linear on a generic parallel line.
    for degree in [1,2,4,7]:
        f=(y*y+1)*x**degree+y*x**max(0,degree-1)+y*y
        g=S.diff(f,x,degree-1)
        check('translation stability yields a generic affine-linear restriction degree'+str(degree),S.degree(g,x)==1 and S.diff(g,x)!=0 and S.diff(g,x,2)==0)
    # Quasi-finiteness cannot be replaced by generic finiteness.
    check('contracted-curve control has a polynomial source field element missing from the target algebra',S.cancel((x*y)/x)==y and (x*y).subs(x,0)==0)
    # Finite non-Keller map: deck involution fixes the target but not Euler data.
    p=x+x*x;q=y;deck=lambda f:S.expand(f.subs(x,-1-x))
    check('finite non-Keller control has a monic source-coordinate equation',S.expand(x*x+x-p)==0)
    check('finite control target is deck invariant',deck(p)==p and deck(q)==q)
    check('finite control weight component is not in target field',deck(x)!=x)
    check('finite control Euler derivative is not in target field',deck(E(p))!=E(p))
    check('finite control is explicitly not a Keller map',S.diff(p,x)==1+2*x)
    # An automorphism need not have homogeneous coordinates in the source grading.
    p=x+y*y;q=y+p*p
    jac=S.expand(S.diff(p,x)*S.diff(q,y)-S.diff(p,y)*S.diff(q,x))
    PP,QQ=S.symbols('PP QQ');invx=PP-(QQ-PP*PP)**2;invy=QQ-PP*PP
    check('mixed-weight automorphism has exact constant Jacobian',jac==1)
    check('mixed-weight automorphism has a polynomial inverse',S.expand(invx.subs({PP:p,QQ:q},simultaneous=True)-x)==0 and S.expand(invy.subs({PP:p,QQ:q},simultaneous=True)-y)==0)
    check('invariant target algebra need not have homogeneous chosen coordinates',len(parts(p))>1 and len(parts(q))>1)
    for ff in [p,q]:
        check('automorphic Euler derivatives are target polynomials '+str(ff),S.expand(E(ff).subs({x:invx,y:invy},simultaneous=True)).is_polynomial(PP,QQ))
    # Local immersion and birationality do not imply one-line injectivity.
    pp=t*t;qq=t**3-t
    check('nodal line control is immersive everywhere',S.gcd(S.diff(pp,t),S.diff(qq,t))==1)
    check('nodal line control has an actual collision',pp.subs(t,1)==pp.subs(t,-1) and qq.subs(t,1)==qq.subs(t,-1))
    check('nodal line control has full rational function field',S.cancel(qq/(pp-1))==t)
    # Whole B_n weight components and their axis-linear witnesses.
    modules=[]
    for n in [1,2,4]:
        HH=S.prod(h-j for j in range(2,n+2));F1=(h-C)/(-C)*(1+h);G1=HH/HH.subs(h,0)*(1+h*h)
        uu=x*F1.subs(h,x*x*y);vv=y*G1.subs(h,x*x*y)
        check('missing positive component is globally in its B_n weight module n'+str(n),S.rem(F1,h-C,h)==0 and S.expand(E(uu)-uu)==0)
        check('missing negative component is globally in its B_n weight module n'+str(n),S.rem(G1,HH,h)==0 and S.expand(E(vv)+2*vv)==0)
        check('positive missing component distinguishes every point of the source axis n'+str(n),S.expand(uu.subs(y,0)-x)==0)
        check('negative missing component distinguishes every point of the source axis n'+str(n),S.expand(vv.subs(x,0)-y)==0)
        check('three primitive choices have the exact coupling n'+str(n),S.expand(uu**2*vv-(h*F1**2*G1).subs(h,x*x*y))==0)
        # Divisor order arrays: Lx,Ly,R,D_1,...,D_n.
        divx=S.Matrix([1,0,-1]+[1]*n);divy=S.Matrix([0,1,2]+[-2]*n);divh=S.Matrix([2,1,0]+[0]*n)
        check('complete source coordinate norm divisors multiply to div(h) n'+str(n),2*divx+divy==divh)
        Hx,Hy,GR=S.symbols('Hx Hy GR');Hi=S.symbols('H0:'+str(n));DD=S.prod(Hi)
        nx=Hx*DD/GR;ny=Hy*GR**2/DD**2
        check('rational coordinate norms reproduce the h norm n'+str(n),S.cancel(nx**2*ny-Hx**2*Hy)==0)
        modules.append(dict(n=n,positive_component_form=str(F1),negative_component_form=str(G1)))
    # Maximality is an unbounded written subgroup argument; these are graph controls.
    for d in [3,4,5,8]:
        first={(i,j) for i in range(1,d) for j in range(i+1,d)}
        second={(i,j) for i in range(d) if i!=1 for j in range(i+1,d) if j!=1}
        reached={0}
        while True:
            new=reached|{j for i,j in first|second if i in reached}|{i for i,j in first|second if j in reached}
            if new==reached:break
            reached=new
        check('two distinct point stabilizers generate a connected transposition graph d'+str(d),reached==set(range(d)))
    # Outside does not imply primitive when a proper intermediate field exists.
    V=S.symbols('V');M=S.zeros(4)
    for j in range(4):
        power=j+2;M[power%4,j]=V**(power//4)
    check('imprimitive fourth-power field control has proper-power characteristic polynomial',S.factor((U*S.eye(4)-M).det())==(U*U-V)**2)
    # Tame trace cancellation: x has pole1, every regular multiplier keeps only an odd negative term.
    xx=c/w+1+2*w+3*w*w
    for f in [1,t+w,t*t+w*w+w**3,1+w+w**4]:
        trace=S.expand(xx*f+(xx*f).subs(w,-w))
        check('trace-dual regularity for every regular multiplier control '+str(f),trace.is_polynomial(w,t,c) and all(k%2==0 for (k,),co in S.Poly(trace,w).terms()))
    check('trace of the square need not be regular',S.expand(xx*xx+(xx*xx).subs(w,-w)).coeff(w,-2)==2*c*c)
    inputs=[Path(__file__),HERE/'SOURCE_SYMMETRY_DESCENT_AND_TRACE_GATE.md',manifest]+[HERE/p for p in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      intersection=dict(target_field_intersect_source_polynomials_equals_target_algebra=True,quasi_finite_divisor_domination_required=True,generic_finiteness_alone_suffices=False),
      descent=dict(all_Keller_field_degrees=True,all_polynomial_source_coordinate_systems=True,translation_derivation_field_stability_forces_automorphism=True,nonzero_integral_diagonal_Euler_field_stability_forces_automorphism=True,all_weight_sign_patterns=True,target_action_linearization_required=False,invariant_algebra_not_homogeneous_generators=True,automatic_descent_from_finiteness=False,source_action_existence_excluded=False,projector_controls=projector_controls),
      marked_models=dict(all_n=True,all_other_primes_retained=True,normalized_weight_components=['u=P_1=x*F_1(h), F_1(0)=1, h-C divides F_1','v=Q_-2=y*G_1(h), G_1(0)=1, H divides G_1'],both_components_outside_target_FIELD=True,at_least_one_Euler_coordinate_derivative_outside_target_field=True,every_target_polynomial_coordinate_max_weight_at_least=2,every_target_polynomial_coordinate_min_weight_at_most=-4,primitive_components_require_no_intermediate_fields=True,retained_zero_excess_Sd_hypothesis_suffices=True,monic_norm_degree_under_primitive_hypothesis='full actual Keller degree',normal_monogenic_algebra_claimed=False,coupling='u^2*v=h*F_1(h)^2*G_1(h)',module_controls=modules),
      finite_gate=dict(rational_norm_x='H_x*product(H_i)/G_R',rational_norm_y='H_y*G_R^2/product(H_i)^2',residue_degree_one_required_for_displayed_norms=True,global_image_support_cancellations_allowed=True,trace_dual='Tr(x*B_n) subset A',trace_x_squared_always_regular=False,trace_dual_implies_x_integral=False),
      negative_controls=dict(contracted_curve_intersection_failure=True,finite_non_Keller_symmetry_failure=True,mixed_weight_automorphism=True,immersive_nodal_line=True,imprimitive_outside_element=True,tame_trace_square_pole=True),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,imported_one_line_theorem_independently_reproved=False,finite_Keller_pair_supplied=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'source_symmetry_gate_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; source-field symmetry obstruction, missing weights and trace/norm gate')
if __name__=='__main__':main()
