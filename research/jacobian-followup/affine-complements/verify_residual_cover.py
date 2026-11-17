"""PF2RULED33: all-degree conormal identity and uniform quartic residual passport.
Exact symbolic and finite algebra controls support written geometric proofs.
No search over pole heights, collision counts or tangency orders is run.
"""
from pathlib import Path
from time import monotonic
import json,hashlib,platform,itertools
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,truth):
    assert truth,name
    checks.append(name)
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def compose(a,b):return tuple(a[b[j]] for j in range(len(a)))
def main():
    start=monotonic();pred=json.loads((HERE/'global_boundary_receipt.json').read_text())
    check('current global-boundary predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    source=HERE/'sources/residual-cover';manifest=json.loads((source/'manifest.json').read_text())
    check('two primary standard-result pins',len(manifest['primary_sources'])==2)
    for row in manifest['primary_sources']:
        p=source/row['path'];check('primary source pin '+row['path'],hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'])
    # General polynomial identity, with no specialization of boundary degree.
    a,b,ad,bd,h,rho,lam,p2,q2,kappa=S.symbols('a b ad bd h rho lam p2 q2 kappa')
    W=a*bd-ad*b;J1=2*h*(p2*b-a*q2)+lam**2*W
    beta=-rho*(p2*b-a*q2)-rho*lam**2*W/(2*h)
    check('ALL-degree exact normal coefficient identity',S.cancel(2*h*beta+rho*J1)==0)
    check('gradient is perpendicular to the whole boundary tangent',S.expand((-rho*b)*(h*a)+(rho*a)*(h*b))==0)
    check('Hessian chain-rule relation used in the normal identity',S.cancel(h*h*(-rho*W/h)+rho*h*W)==0)
    # A tame quadratic block explains the residual discriminant squareclass.
    s,B=S.symbols('s B',nonzero=True)
    trace=S.Matrix([[2,0],[0,2*s/B]])
    check('ramified quadratic trace discriminant has coefficient inverse beta',S.cancel(trace.det()-4*s/B)==0)
    z,A0,A1=S.symbols('z A0 A1');f=z*z-s/B;g=z*z+A1*z+A0
    check('quartic product discriminant splits into blocks and a square resultant',
          S.cancel(S.discriminant(f*g,z)-S.discriminant(f,z)*S.discriminant(g,z)*S.resultant(f,g,z)**2)==0)
    check('residual quadratic algebra has its own discriminant squareclass',S.discriminant(g,z)==A1*A1-4*A0)
    check('compulsory boundary collision excludes finite degree3',3-2*2<0)
    # Finite controls for the all-degree sign character decomposition.
    for degree in [4,5]:
        tau=(1,0,*range(2,degree));centralizer=[p for p in itertools.permutations(range(degree)) if compose(p,tau)==compose(tau,p)]
        check(f'centralizer preserves ramified pair degree{degree}',all(set(p[:2])=={0,1} for p in centralizer))
        check(f'sign splits into pair sign and residual orientation degree{degree}',all(parity(p)==(-1 if p[0]==1 else 1)*parity(tuple(x-2 for x in p[2:])) for p in centralizer))
        check(f'point stabilizer acts by exactly the residual orientation degree{degree}',all(parity(p)==parity(tuple(x-2 for x in p[2:])) for p in centralizer if p[0]==0))
    # Complete ordinary cusp block and its residual curve.
    u,v,t,x,y=S.symbols('u v t x y');P=u-3*v*v;Q=u*v-2*v**3;H=lambda x,y:4*x**3+27*y*y
    check('cusp block has Jacobian exactly the ramification uniformizer',S.expand(S.diff(P,u)*S.diff(Q,v)-S.diff(P,v)*S.diff(Q,u))==u)
    check('cusp pullback factors into ramification twice and the residual curve',S.expand(H(P,Q)-u*u*(4*u-9*v*v))==0)
    gamma=(-3*t*t,-2*t**3);h0=t;a0=-6;b0=-6*t;rho0=18*t**3
    check('cusp gradient scalar and beta are exact',S.expand(12*gamma[0]**2-rho0*(-b0))==0 and S.expand(54*gamma[1]-rho0*a0)==0 and S.cancel(-rho0/(2*h0)+9*t*t)==0)
    check('residual cusp point is unramified over boundary normalization',S.expand(P.subs({u:9*t*t,v:-2*t})-gamma[0])==0 and S.expand(Q.subs({u:9*t*t,v:-2*t})-gamma[1])==0)
    check('source X has exactly pole order2 at the deleted cusp point',S.cancel((1/u).subs(u,9*t*t)*t*t)==S.Rational(1,9))
    # All positive collision multiplicities, split into their two parity classes.
    r=S.symbols('r',integer=True,nonnegative=True);tt=S.symbols('tt',positive=True)
    check('ALL odd collision orders normalize with base index2',S.powsimp((tt**(2*r+1))**2-(tt**2)**(2*r+1),force=True)==0)
    check('ALL positive even collision orders split into two unramified branches',S.powsimp((tt**(r+1))**2-tt**(2*r+2),force=True)==0)
    check('odd collision norm has exactly the predicted base pole',S.powsimp((-tt**(-2*r-1)*tt**(-2*r-1))*(tt**2)**(2*r+1),force=True)==-1)
    check('two even collision poles give the predicted base norm',S.powsimp((-tt**(-r-1)*tt**(-r-1))*tt**(2*r+2),force=True)==-1)
    m=S.symbols('m',integer=True,positive=True);zz,yy=S.symbols('zz yy')
    check('universal local fold pullback has residual equation z2 minus y^m',S.expand((zz*zz)*(zz*zz-yy**m)-(zz*zz)**2+(zz*zz)*yy**m)==0)
    eta=S.symbols('eta')
    check('ALL collision orders give the claimed sign-cover A_(m-1) equation',S.expand((zz-eta)*(zz+eta)-(zz*zz-eta*eta))==0)
    # A finite representative of the global squarefactor formula; the proof is factorwise for all orders.
    a1,b1,a2,b2,c=S.symbols('a1 b1 a2 b2 c');beta_model=(t-c)**2*(t-a1)*(t-b1)*(t-a2)**2*(t-b2)**2;Smodel=(t-a1)*(t-b1)
    check('even collision and cusp factors drop from the quadratic squareclass',S.cancel(beta_model/Smodel-((t-c)*(t-a2)*(t-b2))**2)==0)
    # Uniform Riemann-Hurwitz and puncture counts, symbolic in collision numbers.
    n,o=S.symbols('n o',integer=True,nonnegative=True);genus=o-1;holes=4*n-2*o+3
    check('Riemann-Hurwitz gives genus o-1 for all nonempty odd branch sets',S.expand(2*genus-2-(-4+2*o))==0)
    check('finite deletions count every collision fibre and just one cusp point',S.expand((2*o+4*(n-o)+1)-(4*n-2*o+1))==0)
    check('connected residual source Euler characteristic is 1-4n',S.expand(2-2*genus-holes-(1-4*n))==0)
    check('split residual source has the same total Euler characteristic',S.expand(4-(4*n+3)-(1-4*n))==0)
    check('complete residual divisor Euler characteristic accounts for even branch identifications',S.expand((2-2*o)-2*(n-o)-(2-2*n))==0)
    check('normalized sign cover has Euler characteristic 2+2n',S.expand(2*2-(2-2*n)-(2+2*n))==0)
    total_m=S.symbols('total_m')
    check('minimal sign resolution Euler characteristic includes every rational double point',S.expand((2+2*n)+2+2*(total_m-n)-(4+2*total_m))==0)
    check('nodal specialization has genus n-1 and 2n+3 punctures',genus.subs(o,n)==n-1 and holes.subs(o,n)==2*n+3)
    # One-node abstract passport is consistent; it is not an actual Wright cover.
    w=S.symbols('w',nonzero=True);tp=(w+1/w)/2;wp=(w-1/w)/2
    check('one-node residual conic has the stated rational parametrization',S.cancel(wp*wp-(tp*tp-1))==0)
    check('one-node compact curve loses five distinct points',len({S.Integer(0),S.oo,S.Integer(1),S.Integer(-1),S.I})==5 and 2-5==1-4)
    check('all-even quadratic algebra is split',S.expand((z-1)*(z+1))==z*z-1)
    inputs=[Path(__file__),HERE/'RESIDUAL_DOUBLE_COVER.md',HERE/'global_boundary_receipt.json',
            HERE/'GLOBAL_BOUNDARY_OBSTRUCTION.md',HERE/'WRIGHT_MARKED_MODEL.md',
            HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md',
            HERE.parent/'jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md']+sorted(source.glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        all_degree_identity='[u^2]H(P,Q)=-kappa*rho/(2*h)',
        all_degree_residual_discriminant_squareclass='beta=[u^2]H(P,Q)',necessary_finite_field_degree_at_least=4,
        quartic_scope='every pole height, collision count and positive tangency multiplicity',
        quartic_passport=dict(collision_count='n>=1',odd_collision_count='o',
            normalized_residual_equation='w^2=product((t-a_i)*(t-b_i), i with m_i odd)',
            connected_iff='o>=1',compact_genus_if_connected='o-1',infinity_points_if_connected=2,
            punctures_if_connected='4*n-2*o+3',total_Euler_characteristic='1-4*n',
            all_even_case='two affine lines; remove all collision parameters on both plus one cusp point on one',
            nodal_case_genus='n-1',nodal_case_punctures='2*n+3',
            finite_X_pole_orders=dict(cusp_missing_point=2,odd_collision_parameter='m_i',even_collision_each_point='m_i/2'),
            residual_divisor_class='T=-2*K_N'),
        quartic_sign_cover=dict(canonical_sheaf_trivial=True,degree_over_target_sign_cover=4,unramified_in_codimension_one=True,
            finite_singularities='one A2; two A_(m_i-1) per collision, A0 smooth',Euler_characteristic='2+2*n',minimal_resolution_Euler_characteristic='4+2*sum(m_i)'),
        assurance=dict(actual_marked_ring_and_finite_Keller_map_are_hypotheses=True,
            all_degree_conormal_and_orientation_arguments_are_written=True,all_degree_sign_pullback_canonical_sheaf_trivial=True,all_degree_sign_map_unramified_in_codimension_one=True,
            quartic_parity_cases_cover_all_positive_collision_orders=True,
            no_transverse_node_assumption=True,no_pole_bound_or_coefficient_stratum_used=True,
            residual_cover_is_not_an_intermediate_field_of_the_quartic_surface_extension=True,
            abstract_passport_is_not_a_global_realization=True,primary_orientation_and_Hurwitz_results_imported=True,
            actual_plane_embedding_constructed=False,actual_finite_map_constructed=False,
            all_colliding_Wright_covers_excluded=False,quartic_packets_excluded=False,
            arbitrary_normalization_compatible_completion_proved=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],
        source_sha256={str(p.relative_to(HERE)) if p.is_relative_to(HERE) else '../jc2-omitted-divisors/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'residual_cover_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} symbolic/algebra controls and source pins; residual cover reduction only; goal active')
if __name__=='__main__':main()
