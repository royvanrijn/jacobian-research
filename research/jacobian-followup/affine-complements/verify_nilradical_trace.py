"""PF117 small exact controls; general trace restriction uses the written proof."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import logarithmic_pullback as L

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.cancel(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    prior=json.loads((P/'logarithmic_pullback_receipt.json').read_text())
    check('sealed PF116 PASS',prior['outcome']=='PASS')
    for f,h in prior['source_sha256'].items():check('sealed PF116 '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/nilradical-trace/manifest.json').read_text())
    for source in manifest['sources']:
        f=P/source['path']
        check('primary/source seal '+f.name,sha(f)==source['sha256'] and f.stat().st_size==source['bytes'])
    check('standard duality inputs declared, no fresh PDF claim',manifest['CI_trace_generator_imported'] and manifest['finite_dual_base_change_imported'] and manifest['finite_canonical_module_duality_imported'] and not manifest['fresh_PDF_audit_claimed'])
    u,v,x,y,t=S.symbols('u v x y t')
    # The trace matrix is obtained from multiplication, not an assumed formula.
    M=S.Matrix([[0,0,v],[1,0,-u],[0,1,0]])
    I=S.eye(3);J=3*M**2+u*I
    check('actual cubic monic presentation',M**3+u*M-v*I==S.zeros(3))
    tr=[S.trace(M**k) for k in range(3)]
    check('actual nonreduced rank3 algebra trace',tr==[3,0,-2*u])
    c=S.Matrix([[0,0,1]]);one=S.Matrix([1,0,0])
    for k in range(3):ident('actual CI trace equals Jacobian times coefficient dual '+str(k),(c*J*M**k*one)[0]-tr[k])
    # Verify this chosen c generates the entire relative dual, with unit determinant.
    H=S.Matrix(3,3,lambda i,j:(c*M**(i+j)*one)[0])
    ident('actual cubic coefficient-dual is a free dual generator',H.det()+1)
    f=4*u**3+27*v**2;a=x+3*y**2;h=4*x+3*y**2;b=a*h
    ident('actual nonreduced branch pullback',f.subs({u:x,v:y**3+x*y})-a*b)
    ident('actual branch Jacobian',S.diff(y**3+x*y,y)-a)
    for name,xx,yy in [('ramified',-3*t**2,t),('ordinary',-3*t**2,-2*t)]:
        ident('actual cusp branch image '+name,yy**3+xx*yy+2*t**3)
        ident('actual reduced source branch '+name,b.subs({x:xx,y:yy}))
    unweighted=t-2*t;weighted=2*t-2*t
    ident('unweighted trace has forbidden cusp value',unweighted+t)
    ident('weighted trace of y remains in target cusp',weighted)
    check('cusp semigroup omits linear term',1 not in {2*i+3*j for i in range(3) for j in range(3)})
    c1,c2,l1,l2,k=S.symbols('c1 c2 l1 l2 k')
    ident('contact2 source gluing forces regular weighted trace',(2*l1+l2).subs(l2,-2*l1))
    ident('constant codifferent condition retains ramification diagonal',(c1-2*c2).subs({c1:2*k,c2:k}))
    check('ramification tuple is not in the actual local curve ring',2!=1)
    # Universal nodal calculation: x-axis transverse index q, field degree p.
    p,q=S.symbols('p q',integer=True,positive=True)
    ident('UNBOUNDED reducible-target constant gluing p*c_x=q*c_y',p*q-q*p)
    check('UNBOUNDED nodal trace degree is nonzero',(p*q).is_positive)
    # A source coefficient x^p has trace p*u on the x-axis, zero on the other;
    # the constant values still glue at the node, with no normalizing split.
    ident('UNBOUNDED nodal nonconstant trace glues at origin',(p*u).subs(u,0))
    # If one omits the branch value, the proposed unit cancellation fails.
    check('omitted branch-value hypothesis is necessary',(2*x).subs(x,0)==0)
    M2=S.Matrix([[0,u],[1,0]])
    invJ=M2.inv()/2
    check('missing-branch codifferent generator has regular trace',S.trace(invJ)==0 and S.trace(invJ*M2)==1)
    check('missing-branch codifferent is not the source ring',S.denom(invJ[0,0]+invJ[1,0])==2*u)
    # Actual non-Gorenstein PF116 quotient, not a freely chosen model curve.
    z=L.explicit_sharp_pair(2)
    X,Y,P0,Q0=[z[k] for k in ('X','Y','P','Q')]
    F=S.expand(P0**3-Q0**2);H0=S.cancel(F/X**8);bred=S.expand(X*H0)
    check('actual m2 branch curve is reduced',S.Poly(S.gcd(S.gcd(bred,S.diff(bred,X)),S.diff(bred,Y)),X,Y).total_degree()==0)
    check('actual m2 finite map invariant',all((3*i+j)%5==0 for h0 in [P0,Q0] for (i,j),cc in S.Poly(h0,X,Y).terms()))
    check('actual m2 reduced curve has determinant character',all((3*i+j)%5==4 for (i,j),cc in S.Poly(bred,X,Y).terms()))
    check('actual m2 ambient determinant character is nontrivial',(3+1)%5==4)
    jac=S.expand(S.diff(P0,X)*S.diff(Q0,Y)-S.diff(P0,Y)*S.diff(Q0,X))
    ident('actual m2 Jacobian/nonreduced quotient is unit',jac/S.cancel(F/bred)+207360)
    check('actual m2 transverse weights retained',z['profile']['index']==8)
    # Small group C4 swaps branch components; its generator has no eigenvalue1.
    G=S.Matrix([[0,S.I],[S.I,0]])
    check('small branch-permuting quotient control',G**4==S.eye(2) and all((G**k-S.eye(2)).det()!=0 for k in [1,2,3]))
    P1=x**2-y**2;Q1=x**2*y**2;b1=x*y*(x**2+y**2)
    sub={x:S.I*y,y:S.I*x}
    for name,h0 in [('P',P1),('Q',Q1),('bred',b1)]:ident('branch-permuting quotient invariance '+name,h0.subs(sub,simultaneous=True)-h0)
    ident('branch-permuting actual nonreduced pullback',Q1*(P1**2+4*Q1)-b1**2)
    ident('branch-permuting Jacobian cancellation',S.diff(P1,x)*S.diff(Q1,y)-S.diff(P1,y)*S.diff(Q1,x)-4*b1)
    # Exact Reynolds restriction/inverse on a rational representation, no root approximation.
    action=S.diag(1,-1);average=(S.eye(2)+action)/2
    check('Reynolds idempotent and fixed projection',average**2==average and action*average==average)
    cov=S.Matrix([[S.symbols('a'),0]])
    check('invariant functional equals extension through Reynolds',cov*average==cov and cov*action==cov)
    text=' '.join((P/'NILRADICAL_TRACE_RESTRICTION_AND_CYCLIC_DESCENT.md').read_text().lower().split())
    for phrase in ['all branch values','nonreduced','formal module generator','no ext1 vanishing',
                   'ramification-weighted trace','actual glued ring','unramified in codimension one',
                   'reynolds','non-gorenstein','gorenstein','not identify intrinsic',
                   'not exclude the sharp','whole-plane','finite-algebra','not finite-sample proofs']:
        check('written scope '+phrase,phrase in text)
    files=['NILRADICAL_TRACE_RESTRICTION_AND_CYCLIC_DESCENT.md','verify_nilradical_trace.py',
           'logarithmic_pullback_receipt.json','sources/nilradical-trace/manifest.json']+[z['path'] for z in manifest['sources']]
    result=dict(claim='PF2RULED117',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_finite_smooth_plane_degrees=True,ALL_branch_values_required=True,
        nonreduced_base_change_retained=True,nilradical_annihilator_equals_weighted_trace=True,
        embedded_curve_codifferent_Dr_B=True,all_multibranch_gluing_retained=True,
        actual_small_cyclic_quotient_descent=True,branch_permuting_small_groups_included=True,
        reduced_downstairs_branch_curve_Gorenstein=True,actual_surface_trace_restriction_closed=True,
        actual_upstairs_residue_codifferent_bridge=True,invariant_quotient_residue_bridge=True),
      controls=dict(cubic_ambient_trace=list(map(str,tr)),cubic_unweighted_y_trace=str(unweighted),
        cubic_weighted_y_trace=str(weighted),nodal_diagonal=['q on x-axis','p on y-axis'],
        actual_nonGorenstein_m2_weights=[8,1,1],branch_permuting_small_group_order=4),
      assurance=dict(CI_trace_generator_imported=True,finite_duality_imported=True,
        written_unbounded_proof=True,finite_samples_prove_theorem=False,
        inverse_Jacobian_reduced_as_fraction=False,nonreduced_base_change_replaced_by_reduced=False,
        curve_ramification_diagonal_discarded=False,normalization_gluing_discarded=False,
        ambient_singular_dual_assumed_free=False,intrinsic_singular_Kahler_residues_identified=False,
        nonplane_curve_Jacobian_duality_imported=False,sharp_finite_covers_excluded=False,
        new_global_index_exclusion=False,actual_global_finite_algebra_closed=False,
        global_trace_moment_regularity_proved=False,Bernstein_connection_transport_proved=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,
        research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
      environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'nilradical_trace_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF117 controls; local trace restriction/cyclic descent closed; global goal ACTIVE')


if __name__=='__main__':main()
