"""Exact small controls for formal binary GVC; universal proof is written."""
from pathlib import Path
from math import comb,factorial
import hashlib,json
import sympy as s


def field_cutoff(D):
    if not isinstance(D,int) or D<0:raise ValueError('Nonnegative integer degree required')
    return 1+(D+1)**2//4


def core_order(D,C=None):
    if not isinstance(D,int) or D<1:raise ValueError('Positive phase-degree bound required')
    C=field_cutoff(D) if C is None else C
    if not isinstance(C,int) or C<1:raise ValueError('Positive integer field coefficient required')
    return D*(C*D+1)


def minimal_support(points):
    points=set(points)
    if not points or any(len(p)!=2 or any(not isinstance(v,int) or v<0 for v in p) for p in points):raise ValueError('Nonempty nonnegative integral binary support required')
    return sorted(a for a in points if not any(b!=a and b[0]<=a[0] and b[1]<=a[1] for b in points))


def formal_geometric_action(poly,x,y,D):
    """Symbol xi+eta^(D+1)/(1-eta), evaluated exactly on ONE polynomial."""
    poly=s.sympify(poly)
    if poly.has(s.Float):raise ValueError('Exact polynomial coefficients required')
    p=s.Poly(poly,x,y)
    if p.is_zero:return s.Integer(0)
    if p.total_degree()>30:raise ValueError('Control degree cap')
    return s.expand(s.diff(poly,x)+sum(s.diff(poly,y,k) for k in range(D+1,max(D+1,int(p.degree(y))+1))))


def ordinary_bounds(D,k0,nu,h):
    if any(not isinstance(v,int) or v<1 for v in [D,k0,nu,h]):raise ValueError('Positive integer certificate inputs required')
    C=field_cutoff(D);R=core_order(D)
    if D*k0<R:raise ValueError('The proposed seed jet does not contain the core')
    E=(C*D+1)*(h*nu-1)
    return dict(E=E,K=max(k0,E),N=C+E,symbol_jet=D*max(k0,E))


def run():
    x,y=s.symbols('x y');thresholds=[]
    for D in range(1,7):
        C=field_cutoff(D);R=core_order(D)
        assert R+1>D*(C*D+1)
        thresholds.append(dict(degree=D,optimal_field_coefficient=C,field_core_jet=R,tail_minimum_order=R+1))
    envelopes=[]
    for points in [[(1,0),(0,4),(2,1),(0,7),(4,4)],[(0,8),(1,6),(3,2),(5,0),(6,1),(3,4)],[(2,3),(3,3),(4,2),(4,5)]]:
        mins=minimal_support(points)
        for slope in [s.Rational(1,3),s.Rational(1),s.Rational(3,2),s.Rational(7)]:
            assert min(slope*a+b for a,b in points)==min(slope*a+b for a,b in mins)
        envelopes.append(dict(support=points,minimal_points=mins))
    formal=[]
    for D in [1,2,3]:
        P=y**D;C=field_cutoff(D)
        pure=[]
        for m in range(1,5):
            out=P**m
            for _ in range(m):out=formal_geometric_action(out,x,y,D)
            assert out==0;pure.append(m)
        mixed=[]
        for m in [C+1,C+2]:
            out=x*P**m
            for _ in range(m):out=formal_geometric_action(out,x,y,D)
            assert out==0;mixed.append(m)
        formal.append(dict(phase_degree=D,symbol='xi+eta^(D+1)/(1-eta)',pure_indices=pure,mixed_zero_indices=mixed))
    sharp=[]
    for D in range(1,7):
        e=(D+1)//2;j=D+1-e;N=e*j+1;q=1;m=N*q
        a=(1+(e-1)*j)*q;b=j*q
        assert a+b==m and e*a==q+(e-1)*m and N*b==j*m
        assert N==field_cutoff(D)
        coefficient=comb(m,a)*factorial(q+(e-1)*m)*factorial(j*m)
        assert coefficient>0
        sharp.append(dict(phase_degree=D,e=e,j=j,last_nonzero_index=m,coefficient=str(coefficient)))
    # The symbolic ring-bound calculator uses supplied certificates. These sample
    # integers are arithmetic controls, NOT claimed universal k0/nu values.
    arithmetic=[]
    for D,nu in [(1,1),(2,2),(3,1)]:
        R=core_order(D);k0=R//D;h=comb(2+D,2)+comb(2+R,2)
        arithmetic.append(dict(degree=D,illustrative_certificate_only=True,values=ordinary_bounds(D,k0,nu,h)))
    # Existing strict-tail boundary: order>degP alone does not preserve purity.
    P=x*y**3
    action=lambda f:s.diff(f,x,2)+s.diff(f,y,5)
    boundary=s.expand(action(action(P**2)));assert boundary==2880*y
    rejected=0
    for action in [lambda:field_cutoff(-1),lambda:core_order(0),lambda:minimal_support([(0,-1)]),lambda:ordinary_bounds(3,1,1,1),lambda:formal_geometric_action(s.Float('0.5')*x,x,y,1)]:
        try:action()
        except ValueError:rejected+=1
        else:raise AssertionError('Expected guard refusal')
    return dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),thresholds=thresholds,
                finite_envelope_controls=envelopes,formal_geometric_tail_controls=formal,sharp_endpoint_controls=sharp,
                illustrative_ring_bound_arithmetic=arithmetic,insufficient_tail_boundary=str(boundary),guard_refusals=rejected,supporting_checker_hashes={},
                field_cutoff_and_core_calculator_implemented=True,ordinary_bound_arithmetic_implemented=True,
                universal_theorem_proved_by_finite_controls=False,universal_seed_or_radical_computed=False,
                least_formal_mixed_cutoff_algorithm_claimed=False,sharp_field_coefficient_transported_to_nonreduced_rings=False,
                full_exceptional_MZ_decided=False,general_formal_eventuality_to_linearity_decided=False,
                independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)


if __name__=='__main__':
    receipt=run();output=Path(__file__).with_name('formal_binary_gvc_verification.json')
    if output.exists():assert json.loads(output.read_text())==receipt,'Preserve immutable receipt'
    else:output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=True,checker_sha256=receipt['checker_sha256'],formal_models=len(receipt['formal_geometric_tail_controls']),sharp_endpoints=len(receipt['sharp_endpoint_controls']),guard_refusals=receipt['guard_refusals'])))
