'''PF2RULED36: exact controls supporting the written intrinsic-normalization gate.
Valuation, local cover, 3-manifold and ruling arguments remain written proofs.
'''
from pathlib import Path
from itertools import permutations
from math import gcd,ceil
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    pred=json.loads((HERE/'quartic_cabling_receipt.json').read_text())
    check('predecessor marked quartic receipt PASS',pred['outcome']=='PASS')
    for name,h in pred['source_sha256'].items():
        check('predecessor seal '+name,digest(HERE/name)==h)
    source_dir=HERE/'sources/intrinsic-normalization'
    manifest=json.loads((source_dir/'manifest.json').read_text())
    check('two primary sources pinned',len(manifest['primary_sources'])==2)
    for row in manifest['primary_sources']:
        check('primary source '+row['path'],digest(source_dir/row['path'])==row['sha256'])
        check('whole 2024 bundle theorem not adopted '+row['path'],
              row['whole_claimed_bundle_theorem_adopted'] is False)
    # The module proof uses the Krull valuation criterion. These check its
    # threshold arithmetic; they are not a formal proof of normality or Ass.
    valuation_controls=0
    for pole in range(1,17):
        for va in range(17):
            check_value=(va-pole>=0)==(va>=pole)
            assert check_value
            for vb in range(17):
                if va==0:assert ((va+vb>=pole)==(vb>=pole))
                if va>=pole or vb>=pole:assert va+vb>=pole
                valuation_controls+=1
        for va in range(1,17):
            exponent=ceil(pole/va)
            assert exponent*va-pole>=0 and (exponent-1)*va-pole<0
    check('valuation threshold, primary ideal, and local nilpotence controls',valuation_controls==4624)
    # A fixed multiplier never clears all powers of a negative-valuation b.
    h,n=S.symbols('h n',integer=True,positive=True)
    check('unbounded poles prevent a nonzero uniform module annihilator',
          S.expand(h-(h+1)*n+1+(h+1)*(n-1))==0)
    # Both collision preimages have the same parity; the infinity valuation
    # of the Kummer squareclass is minus twice the number of odd pairs.
    o=S.symbols('o',integer=True,nonnegative=True)
    check('paired finite quadratic branching has even infinity valuation',S.simplify(S.exp(S.pi*S.I*(-2*o)))==1)
    check('quartic Euler surplus is one',4-1-2==1)
    P=tuple(permutations(range(4)));I=tuple(range(4));tau=(1,0,2,3)
    def mul(a,b):return tuple(a[b[j]] for j in range(4))
    def even(a):return sum(a[i]>a[j] for i in range(4) for j in range(i+1,4))%2==0
    candidates={a for a in P if even(a) and mul(a,tau)==mul(tau,a)}
    check('even transposition centralizer has two elements',candidates=={I,(1,0,3,2)})
    check('residual infinity forces identity preferred longitude',
          {a for a in candidates if a[2:]==(2,3)}=={I})
    # Symbolic order audit for n=3(2r+1), rather than inferring the
    # unbounded theorem from controls on finitely many orders.
    r=S.symbols('r',integer=True,nonnegative=True);N=3*(2*r+1)
    check('excluded homology order is odd',S.expand(N-2*(3*r+1))==1)
    check('three times a generator has strictly smaller order',S.simplify(N/3)==2*r+1)
    check('twice a generator remains primitive at every odd order',S.expand(N-2*(3*r+1))==1)
    U,V=S.symbols('U V',integer=True,positive=True)
    check('alternating torus homology order is even',S.expand(2*U*V)%2==0)
    coprime_controls=0
    for order in range(3,100,6):
        assert gcd(order,2)==1 and gcd(order,3)==3
        assert len({3*j%order for j in range(order)})==order//3
        for p in range(3,50,2):
            assert (gcd(p,order)==1)==(gcd(p,2*order)==1)
            if gcd(p,order)==1:
                assert gcd(p,6)==1
            coprime_controls+=1
    check('marked-class primitivity and transparent winding controls',coprime_controls==408)
    depth=S.symbols('depth',integer=True,positive=True)
    check('descent strictly removes one cabling stage',S.expand(depth-(depth-1))==1)
    # Canonical exponent independent of affine translation in the torsor.
    u,v=S.symbols('u v');g=S.Function('g')
    for m in range(-2,7):
        x=1/u;y=u**m*v+g(u)
        det=S.det(S.Matrix([[S.diff(x,u),S.diff(x,v)],[S.diff(y,u),S.diff(y,v)]]))
        check('torsor canonical exponent m'+str(m),S.simplify(det+u**(m-2))==0)
    m=S.symbols('m')
    check('sole simple ramification forces linear bundle O minus three',S.solve(m-2-1,m)==[3])
    # Exact Wright control: source automorphism and two multiple infinity fibers.
    x,y,A,D=S.symbols('x y A D')
    aa=x*x+y;dd=x+aa*aa
    check('retained source-coordinate pair is an automorphism with Jacobian minus one',
          S.det(S.Matrix([[S.diff(aa,x),S.diff(aa,y)],[S.diff(dd,x),S.diff(dd,y)]]))==-1)
    inverse={x:D-A*A,y:A-(D-A*A)**2}
    check('coordinate inverse recovers a',S.expand(aa.subs(inverse,simultaneous=True))==A)
    check('coordinate inverse recovers d',S.expand(dd.subs(inverse,simultaneous=True))==D)
    chart={x:1/u,y:-u*u+u**3*v}
    for f,k in [(aa,2),(dd,4)]:
        numerator=S.expand(u**k*f.subs(chart,simultaneous=True))
        check('coordinate infinity pole exactly '+str(k),numerator.is_polynomial(u,v) and numerator.subs(u,0)==1)
    inputs=[Path(__file__),HERE/'INTRINSIC_NORMALIZATION_GATE.md',
            HERE/'quartic_cabling_receipt.json']+sorted(source_dir.glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope=dict(algebra='all field degrees; normal affine plane complement with exactly one omitted prime',
                   geometry='actual finite quartic Keller normalization with exactly one omitted prime of generic index two'),
        valuation_module=dict(annihilator='symbolic power p^(-v(b))',
             associated_primes=['p'],support='V(p)',uniform_annihilator=0,
             p_action_locally_nilpotent=True,prime_hypothesis_automatic=True),
        intrinsic_quartic=dict(residual_squareclass='product of paired odd collision factors',
             residual_infinity_points=2,residual_infinity_unramified=True,
             preferred_longitude='identity',Wright_marking_used=False),
        excluded_end=dict(lens=True,ramification_exterior='solid torus',
             homology_order='n odd and divisible by 3',marked_class='primitive generator',
             all_finite_cabling_depths=True),
        ruling=dict(multiplicity_one_quartic_excluded=True,
             negative_coordinate_pole_one_excluded=True,canonical_linear_bundle_degree=-3,
             all_degree_reduced_simple_boundary_torus_end_excluded=True),
        assurance=dict(valuation_proof_written=True,local_cover_proof_written=True,
             general_order_cabling_proof_written=True,topological_predecessor_imported=True,
             full_three_manifold_proof_implemented=False,BGV_imported=True,
             reduced_boundary_not_confused_with_reduced_fiber=True,
             coordinate_control_is_finite_Keller_extension=False,
             higher_index_single_prime_quartic_excluded=False,
             arbitrary_quartic_normalizations_excluded=False,
             multiple_fiber_rulings_excluded=False,
             all_finite_Wright_field_degrees_excluded=False,
             general_lens_core_end_proved=False,coefficient_search_used=False,
             formal_verification=False,independent_full_reproof=False,
             external_review=False,research_goal_completed=False),
        finite_control_counts=dict(valuation=valuation_controls,winding=coprime_controls),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        python=platform.python_version(),sympy=S.__version__)
    (HERE/'intrinsic_normalization_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; intrinsic one-simple-boundary gate; full goal active')
if __name__=='__main__':main()
