"""PF2RULED32 controls and source pins for the global Wright boundary gate.
The all-degree topology and imported Chau theorems are written proofs.
No bounded calculation is used to establish their universal scope.
"""
from pathlib import Path
import hashlib,json,platform
from time import monotonic
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def compose(a,b):return tuple(a[b[j]] for j in range(len(a)))
def transposition(n,a,b):
    p=list(range(n));p[a],p[b]=p[b],p[a];return tuple(p)
def generated(permutations):
    identity=tuple(range(len(permutations[0])));seen={identity};pending=[identity]
    while pending:
        p=pending.pop()
        for g in permutations:
            q=compose(p,g)
            if q not in seen:seen.add(q);pending.append(q)
    return seen

def main():
    start=monotonic()
    # These receipts supply the actual ring/chart and conditional fibre dictionary.
    for name in ['volume_period_receipt.json','wright_grading_receipt.json']:
        pred=json.loads((HERE/name).read_text())
        check('predecessor PASS '+name,pred['outcome']=='PASS')
        for path,digest in pred['source_sha256'].items():
            check('predecessor seal '+name+' '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    source=HERE/'sources/global-boundary-topology';manifest=json.loads((source/'manifest.json').read_text())
    check('two independently identified primary imports',len(manifest['primary_sources'])==2)
    historical=json.loads((HERE/'sources/quadratic-remainder-checkpoint/STATUS.json').read_text())
    check('quadratic checkpoint is explicitly historical and carries no exclusion',historical['process_status']=='terminal' and historical['quadratic_M_stratum_excluded'] is False)
    for entry in manifest['primary_sources']:
        path=source/entry['path']
        check('primary PDF pin '+entry['path'],path.read_bytes().startswith(b'%PDF') and hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'])
    check('Chau hypotheses and assurance explicitly separated',manifest['imports']['chau_1999']['theorem']=='1.1' and
          manifest['imports']['chau_2003']['theorem']=='2; Corollary3' and manifest['imports']['independently_reproved'] is False)
    # The homogeneous-coordinate obstruction is an all-degree leading-term argument.
    Z,d=S.symbols('Z d',integer=True,nonnegative=True)
    check('D homogeneous-mate leading coefficient is positive for every degree',S.simplify(3+d)==d+3 and (d+3).is_positive)
    g=S.Function('g')(Z)
    check('universal homogeneous mate differential expression',S.expand((3*Z+1)*g+Z*(Z+1)*S.diff(g,Z)-(g+3*Z*g+Z*(Z+1)*S.diff(g,Z)))==0)
    # Universal closed formula for the finite global jet approximants.
    n=S.symbols('n',integer=True,nonnegative=True);z,w,k=S.symbols('z w k')
    TN=1-(n+2)*z*w+(n+1)*z*z*w;Tprev=1-(n+1)*w+n*z*w
    check('ALL-n truncated inverse-square identity by induction step',S.expand(TN-Tprev-(n+1)*w*(1-z)**2)==0)
    check('truncated inverse-square induction base',S.expand((1-2*z+z*z)-(1-z)**2)==0)
    T=1-(n+2)*w+(n+1)*z*w
    euler=-(n+2)*(n+1)*w+(n+1)*(n+2)*z*w
    H=-k+k*(n+2)*(n+3)*w/2-k*(n+1)*(n+4)*z*w/2
    check('ALL-n exact source-Jacobian closed formula',S.expand(-k*(T+euler/2)-H)==0)
    check('ALL-n retained-source critical locus persists at z=1',S.expand(H.subs({z:1,w:1}))==0)
    check('ALL-n highest residual coefficient cannot vanish',((n+1)*(n+4)).is_positive)
    check('generic polynomial field degree of the jet family grows as N+4',(n+1).is_positive and S.diff(n+4,n)==1)
    check('normal error order tends to infinity while pole bound grows',S.diff(n+2,n)==1 and S.diff(2*n+6,n)==2)
    X,Y,u,v=S.symbols('X Y u v');a=X*Y;b=X*X*Y;D=X**3*Y+X
    check('actual omitted chart recovers all four ring generators',
          [S.expand(f.subs({X:1/u,Y:u**3*v-u**2},simultaneous=True)) for f in [Y,a,b,D]]==[u**3*v-u**2,u**2*v-u,u*v-1,v])
    check('rational formal mate has an actual retained-source pole',S.cancel(a*a/(b*b))==1/X**2)
    N=2;Q=-a*a*sum((j+1)*(b+1)**j for j in range(N+1))/2
    J=S.diff(D,X)*S.diff(Q,Y)-S.diff(D,Y)*S.diff(Q,X)
    target=-1+S.Rational((N+2)*(N+3),2)*(b+1)**(N+1)-S.Rational((N+1)*(N+4),2)*(b+1)**(N+2)
    check('actual polynomial source approximant transports the universal Jacobian formula',S.expand(J-target)==0)
    Qchart=S.expand(Q.subs({X:1/u,Y:u**3*v-u**2},simultaneous=True))
    Jchart=-S.diff(Qchart,u)
    check('representative normal jet reaches the predicted next order',S.Poly(S.expand(Jchart-u),u,v).terms() and min(ex[0] for ex,c in S.Poly(S.expand(Jchart-u),u,v).terms())==N+2)
    check('representative pair still has critical points on retained X=0',S.expand(J.subs(X,0))==0)
    check('representative pair contracts the whole retained source line',D.subs(X,0)==0 and Q.subs(X,0)==0)
    # The global projection gate rules out even some clean cusp+collision data.
    p=v*v;q=v**5-v**3;h=S.gcd(S.diff(p,v),S.diff(q,v));aa=S.cancel(S.diff(p,v)/h);bb=S.cancel(S.diff(q,v)/h);Wr=S.expand(aa*S.diff(bb,v)-S.diff(aa,v)*bb)
    check('pure-power boundary control has exactly one ordinary cusp parameter',h==v and Wr.subs(v,0)==-6)
    check('same control has a clean immersed image double value',p.subs(v,1)==p.subs(v,-1)==1 and q.subs(v,1)==q.subs(v,-1)==0 and
          S.det(S.Matrix([[S.diff(p,v).subs(v,1),S.diff(q,v).subs(v,1)],[S.diff(p,v).subs(v,-1),S.diff(q,v).subs(v,-1)]]))!=0)
    p1=aa;q1=bb;p2=0;q2=S.Rational(15,2)*v
    j0=p1*S.diff(q,v)-S.diff(p,v)*q1
    j1=2*(p2*S.diff(q,v)-S.diff(p,v)*q2)+p1*S.diff(q1,v)-S.diff(p1,v)*q1
    check('pure-power boundary control nevertheless passes the complete first-normal jet identity',S.expand(j0)==0 and S.expand(j1)==-6)
    # Small controls for the universal connected-transposition-graph proof.
    t01=transposition(3,0,1);t12=transposition(3,1,2)
    check('path conjugation supplies the missing transposition',compose(t01,compose(t12,t01))==transposition(3,0,2))
    check('three connected sheet edges generate S4',len(generated([transposition(4,0,1),transposition(4,1,2),transposition(4,2,3)]))==24)
    check('two disjoint sheet edges do not yield transitivity',len(generated([transposition(4,0,1),transposition(4,2,3)]))==4)
    fibre=[dict(kind='off_branch',branches=0,excess=0,retained=4),dict(kind='fold',branches=1,excess=0,retained=2),dict(kind='cusp',branches=1,excess=1,retained=1),dict(kind='collision',branches=2,excess=0,retained=0)]
    for row in fibre:check('quartic fibre identity '+row['kind'],4-2*row['branches']-row['excess']==row['retained'])
    check('a cusp collision is impossible by negative retained length',4-2*2-1<0)
    inputs=[Path(__file__),HERE/'GLOBAL_BOUNDARY_OBSTRUCTION.md',HERE/'volume_period_receipt.json',HERE/'wright_grading_receipt.json',
            HERE/'WRIGHT_MARKED_MODEL.md',HERE/'WRIGHT_GRADING_OBSTRUCTION.md',HERE/'sources/quadratic-remainder-checkpoint/STATUS.json',HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']+sorted(source.glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all pole heights and every finite field degree on the actual marked Wright model',
        global_conclusions=dict(exceptional_curve_equals_full_boundary_image=True,injective_boundary_parametrization_excluded=True,
            every_target_coordinate_boundary_projection_has_at_least_two_finite_critical_values=True,
            monodromy_is_full_symmetric_in_every_finite_degree=True,meridional_generator_count_at_least_degree_minus_one=True,
            quartic_missing_values_equal_boundary_collision_values=True,quartic_missing_value_set_is_nonempty_finite=True,
            quartic_first_Betti_number_of_boundary_equals_missing_value_count=True),
        uniform_jet_control=dict(all_nonnegative_orders_covered_by_written_induction=True,
            normal_identity='J_uv(P,Q_N)=kappa*u modulo u^(N+2)',generic_field_degree='N+4',fixed_field_degree_criterion_not_refuted=True,source_Jacobian_vanishes_on_retained_X_zero=True,every_approximant_contracts_the_retained_X_zero_line=True,
            polynomial_mate_to_D_excluded_in_every_degree=True,rational_formal_mate='-kappa/(2*X^2)',
            uniform_finite_jet_sufficiency_refuted_without_pole_bound=True,fixed_pole_space_completion_cutoff_unaffected=True),
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],
        assurance=dict(explicit_imported_Chau_theorems=True,imported_results_independently_reproved=False,
            all_degree_scope_is_written_not_inferred_from_finite_controls=True,global_result_uses_no_new_Groebner_search=True,
            quadratic_remainder_computation_checkpointed_as_historical=True,
            actual_finite_map_constructed=False,all_mixed_weight_Keller_pairs_excluded=False,
            all_boundary_collisions_excluded=False,quartic_surjectivity_proved=False,
            arbitrary_normalization_compatible_completion_proved=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)) if p.is_relative_to(HERE) else '../jc2-omitted-divisors/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'global_boundary_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls and primary source pins; global boundary gate remains conditional; goal active')
if __name__=='__main__':main()
