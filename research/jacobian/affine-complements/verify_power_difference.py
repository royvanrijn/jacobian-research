"""Exact controls for PF2RULED56; no coefficient or exponent census proves the theorem."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent;checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    preds=['quartic_colouring_receipt.json','sign_square_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor sealed input '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/power-difference/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']+m['background_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('primary/background bytes '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    p,q,a,b=S.symbols('p q a b',positive=True);N=p*q;k=p+q
    # These identities are in unrestricted symbolic exponents and ray weights.
    check('all origin rays in pa<=qb sector',S.expand(N*(a+b)-k*p*a-p*(q*b-p*a))==0)
    check('all origin rays in qb<=pa sector',S.expand(N*(a+b)-k*q*b-q*(p*a-q*b))==0)
    check('all infinity rays in (q-p)b<=qa sector',S.expand(N*a+(q*k-2*N)*b-k*(q-p)*b-p*(q*a-(q-p)*b))==0)
    check('all infinity rays in qa<=(q-p)b sector',S.expand(N*a+(q*k-2*N)*b-k*q*a-q*((q-p)*b-q*a))==0)
    check('generic infinity logarithmic pole bound witness',S.expand(q*k-3*N+N-q*(q-p))==0)
    check('retained exceptional order equals polynomial pullback exponent',N-k==p*q-p-q)
    for pp,qq in [(2,3),(2,5),(3,4),(3,5),(5,7)]:
        check('strict exponent hypothesis control '+str((pp,qq)),S.gcd(pp,qq)==1 and pp*qq>pp+qq and S.Rational(1,pp)+S.Rational(1,qq)<1)
    x,y,t,z,u,v,c,d,g,R,Q=S.symbols('x y t z u v c d g R Q')
    check('quadratic-cubic Taylor unit-ideal identity',S.expand(((1+g*R)**2-(1+g*Q)**3)-g*(2*R-3*Q+g*(R**2-3*Q**2-g*Q**3)))==0)
    gx,gy,Rx,Ry,Qx,Qy=S.symbols('gx gy Rx Ry Qx Qy')
    jac=S.det(S.Matrix([[R*gx+g*Rx,R*gy+g*Ry],[Q*gx+g*Qx,Q*gy+g*Qy]]))
    witness=R*(gx*Qy-gy*Qx)+Q*(Rx*gy-Ry*gx)+g*(Rx*Qy-Ry*Qx)
    check('contracted reduced curve divides coefficient-map Jacobian',S.expand(jac-g*witness)==0)
    # Universal cubic index and trace discriminant under basis change.
    M=S.Matrix([[0,0,-d],[1,0,-c],[0,1,0]]);Mt=u*M+v*M*M;e=S.Matrix([1,0,0]);I=u**3+c*u*v**2+d*v**3
    check('universal binary cubic determinant index',S.expand(S.Matrix.hstack(e,Mt*e,Mt*Mt*e).det()-I)==0)
    Ft=Mt.charpoly(t).as_expr();disc=S.discriminant(Ft,t)
    check('universal cubic index-square discriminant identity',S.expand(disc-(-4*c**3-27*d**2)*I**2)==0)
    scalar=S.symbols('scalar');Ms=Mt+scalar*S.eye(3)
    check('index independent of scalar translation',S.expand(S.Matrix.hstack(e,Ms*e,Ms*Ms*e).det()-I)==0)
    gram=S.diag(1,2,2*g)
    change=S.Matrix([[1,1,1],[1,0,g],[0,1,0]])
    check('normal local cubic product has reduced trace discriminant',gram.det()==4*g)
    check('local generator is valid around branch but has a global index divisor',S.expand((change.T*gram*change).det()-4*g*(1-g)**2)==0 and change.det().subs(g,0)!=0 and change.det().subs(g,1)==0)
    # Exact normalization, immersion and inner-singularity control.
    G=x*x-(x+y*y)**3;xx=(1-t*t)**3;yy=t*(1-t*t);B=x+y*y
    check('A1 A5 torus curve parametrization',S.expand(G.subs({x:xx,y:yy},simultaneous=True))==0)
    check('normalization parameter integral',S.expand((t**3-t+y).subs(y,yy))==0)
    check('normalization parameter has rational inverse',S.cancel((y*(x+y*y)/x).subs({x:xx,y:yy},simultaneous=True)-t)==0)
    check('degree4 auxiliary projection for absolute irreducibility',S.degree((xx+yy*yy).expand(),t)==4 and S.degree(((u-y*y)**2-u**3),y)==4)
    check('normalization is immersive everywhere',S.gcd(S.diff(xx,t),S.diff(yy,t))==1)
    check('one identified pair maps to the sole common coefficient zero',xx.subs(t,1)==xx.subs(t,-1)==yy.subs(t,1)==yy.subs(t,-1)==0 and list(S.groebner([x,B],x,y))==[x,y*y])
    check('both branches smooth with opposite cubic contact coefficients',S.diff(yy,t).subs(t,1)!=0 and S.diff(yy,t).subs(t,-1)!=0 and S.limit(xx/yy**3,t,1)==1 and S.limit(xx/yy**3,t,-1)==-1)
    F=z**3-3*B*z-2*x
    check('A1 torus cubic has exact discriminant',S.expand(S.discriminant(F,z)+108*G)==0)
    check('common zero gives a triple-root fibre',F.subs({x:0,y:0})==z**3)
    check('cubic quotient smooth over the A5 branch singularity',S.diff(F,x).subs({x:0,y:0,z:0})==-2)
    check('all finite plane branch singularities supported at origin',list(S.groebner([G,S.diff(G,x),S.diff(G,y)],x,y))==[2*x-3*y**4,y**5])
    # Dropped-hypothesis controls prevent unwarranted full exclusions.
    Oka=(x*u-2)**2-u**3;OkaF=z**3-3*u*z-2*(x*u-2)
    check('Cstar branch supports the forbidden-in-A1 monogenic model',S.expand(S.discriminant(OkaF,z)+108*Oka)==0 and list(S.groebner([x*u-2,u],x,u))==[1])
    check('noncoprime exponent control really has A1 branch and no coefficient basepoint',S.expand((x+1)**2-(x-1)**2)==4*x and list(S.groebner([x+1,x-1],x))==[1])
    check('reducible branch control cannot invoke the theorem',S.factor(x*x-1)==(x-1)*(x+1))
    check('ordinary commutative cube expansion uses coefficient3',S.expand((R+Q)**3-R**3-Q**3-3*R*Q*(R+Q))==0 and S.expand((R+Q)**3-R**3-Q**3-2*R*Q*(R+Q))!=0)
    inputs=[Path(__file__),HERE/'POWER_DIFFERENCE_AND_CUBIC_INDEX_GATE.md',manifest]+[HERE/n for n in preds]+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']+m['background_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      power_difference=dict(all_polynomial_degrees=True,all_coprime_exponents_at_least_two=True,irreducible_reduced_G_required=True,A1_normalization_required=True,conclusion='V(A,B) nonempty for G=A^p-B^q',dependent_image_case_proved=True,dominant_image_case_proved=True,all_positive_valuation_rays_proved=True,retained_interior_exceptional_checked=True,log_subadditivity_conjecture_used=False),
      actual_quartic=dict(entire_packet=True,all_source_orders=True,all_target_degrees=True,all_rulings=True,both_omitted_primes_retained=True,cubic_in_S4_closure=True,cubic_subfield_of_quartic_asserted=False,sign_triple_etale_at_singularities_required=True,local_cubic_algebra='R x R[z]/(z^2-G)',all_branch_fibres_partition=[2,1],trace_discriminant_reduced_G=True,locally_monogenic_everywhere=True,globally_monogenic_excluded=True),
      index=dict(determinant_line_trivial=True,primitive_integral_element_has_nonconstant_index=True,discriminant_identity='disc(F_t)=constant*G*I(t)^2',support_allowed_to_coincide_with_G=True,free_module_basis_implies_algebra_generator=False,monogenicity_forced_by_marked_source=False),
      controls=dict(immersive_A1_torus_curve_with_inner_A5=True,inner_common_zero_forces_triple_root=True,inner_sign_triple_only_quasietale=True,Cstar_Oka_monogenic_allowed=True,noncoprime_22_countercontrol=True,reducible_G_countercontrol=True,source_printed_cube_coefficient2_not_used=True),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,Pic3_vanishing_proved=False,finite_Keller_pair_supplied=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'power_difference_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(checks)} controls; all-degree power difference and essential cubic index')
if __name__=='__main__':main()
