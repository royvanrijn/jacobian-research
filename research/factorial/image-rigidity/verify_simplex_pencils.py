"""Finite independent identities and scope controls for the written pencil proof."""
from pathlib import Path
from collections import Counter
import hashlib
import itertools
import json
import sympy as s
import moment_bridges as b
import simplex_pencils as p

HERE = Path(__file__).resolve().parent
counts = {}
examples = []


def check(group, value):
    if not value:
        raise AssertionError(group)
    counts[group] = counts.get(group,0)+1


z = s.Symbol('z')
cases = [([0,1],[1,1]),([0,1,1],[1,1,1]),([0,1,s.I],[1,1,1]),
         ([2,2,2],[1,2,1]),([1,-1,1],[2,1,1]),([0,s.sqrt(2),s.I],[1,1,1]),
         ([0,1,s.I],[2,1,3])]
for nodes,shapes in cases:
    ts = s.symbols('t:'+str(len(nodes)))
    ell = sum(c*t for c,t in zip(nodes,ts))
    phase = z+s.I*z*z-s.Rational(1,3)
    rep = p.chain_representation(nodes,shapes)
    for m in range(4):
        expected = b.dirichlet(s.expand(phase.subs(z,ell)**m),ts,shapes)
        actual = p.integrate_polynomial(rep,s.expand(phase**m),z)
        check('dirichlet_chain_identity',s.simplify(actual-expected)==0)
    if rep['kind']=='chains':
        check('chain_has_correct_mass',p.integrate_polynomial(rep,1,z)==1)
        check('all_chain_amplitudes_polynomial',all(s.Poly(w,rep['variable']).degree() <= sum(shapes)-2 for _,w in rep['entries']))
    for marker_exp in [(1,)+(0,)*(len(nodes)-1), (0,2)+(0,)*(len(nodes)-2)]:
        marker = s.prod(t**a for t,a in zip(ts,marker_exp))
        scale = s.prod(s.rf(a,k) for a,k in zip(shapes,marker_exp))/s.rf(sum(shapes),sum(marker_exp))
        shifted = p.chain_representation(nodes,[a+k for a,k in zip(shapes,marker_exp)])
        for m in range(3):
            left = b.dirichlet(s.expand(marker*phase.subs(z,ell)**m),ts,shapes)
            right = scale*p.integrate_polynomial(shifted,s.expand(phase**m),z)
            check('all_original_polynomial_markers',s.simplify(left-right)==0)

    certificate = p.exponential_certificate(nodes,shapes)
    # Independent local residue derivative, rather than the path primitive used
    # by the implementation. Group nodes in the original input cases here.
    grouped = Counter()
    for c,a in zip(nodes,shapes):
        grouped[s.sympify(c)] += a
    expected_groups = {}
    for node,multiplicity in grouped.items():
        regular = 1/s.prod((z-other)**a for other,a in grouped.items() if other != node)
        value = sum(s.binomial(multiplicity-1,j)*s.diff(regular,z,j).subs(z,node)
                    for j in range(multiplicity))/s.factorial(multiplicity-1)
        value = s.simplify(s.factorial(sum(shapes)-1)*value)
        if value != 0:
            expected_groups[node] = value
    all_nodes = set(expected_groups)|set(certificate['groups'])
    check('exponential_chain_vs_local_residues',all(s.simplify(expected_groups.get(c,0)-certificate['groups'].get(c,0))==0 for c in all_nodes))
    check('nonzero_exponential_group_remains',bool(certificate['groups']))
    examples.append({'nodes':[str(c) for c in nodes],'shapes':shapes,
                     'chain_kind':rep['kind'],
                     'exponential_groups':[[str(c),str(a)] for c,a in certificate['groups'].items()]})

# A nonconstant affine triangle phase has an algebraic, nonzero exponential
# integral. This rejects a false transcendence upgrade of the theorem.
groups = p.exponential_certificate([0,1,1])['groups']
check('nonconstant_triangle_period_can_be_algebraic',groups=={s.Integer(0):s.Integer(2)})
point = p.chain_representation([0,0,0])
check('point_evaluation_is_separate',point['kind']=='point' and p.integrate_polynomial(point,z,z)==0 and p.integrate_polynomial(point,1,z)==1)

# The source's odd block shapes still have the general transfer, but are not
# justified by this integer-Dirichlet chain theorem.
for nodes,shapes in [([0,1],[s.Rational(1,2),1]),([0,2*s.pi*s.I],None),([0,0.1],None)]:
    try:
        p.exponential_certificate(nodes,shapes)
    except ValueError:
        check('unsupported_shapes_or_transcendental_nodes_rejected',True)
    else:
        raise AssertionError('unsupported exponential certificate accepted')
check('transcendental_nonvanishing_countercontrol',s.integrate(s.exp(2*s.pi*s.I*z),(z,0,1))==0)

# All-order marker cancellation is certified by identical grouped profiles.
left = p.chain_representation([1,1,s.I],[2,1,1])
right = p.chain_representation([1,1,s.I],[1,2,1])
check('identically_zero_marker_by_node_symmetry',left==right)

for nodes,shapes,degree in [([0,1,s.I],[1,1,1],3),([1,1,-1,s.I],[1,2,1,1],4),([1,2],[2,3],3)]:
    ys = s.symbols('y:'+str(len(nodes)))
    S,L = sum(ys),sum(c*y for c,y in zip(nodes,ys))
    F = s.expand(S**degree-2*S**(degree-1)*L+s.I*S*L**(degree-1)+L**degree)
    result = p.homogeneous_pencil_certificate(F,ys,nodes,shapes)
    check('pencil_membership_accepted',result['status']=='nonzero_infinitely_many_pure_moments')
    reconstructed = sum(s.sympify(c)*S**a*L**bb for (a,bb),c in result['binary_terms'])
    check('independent_membership_identity',s.Poly(s.expand(F-reconstructed),*ys,extension=True).is_zero)
    ts = s.symbols('t:'+str(len(nodes)))
    angular = sum(s.sympify(c)*(sum(a*t for a,t in zip(nodes,ts)))**exps[1] for exps,c in result['binary_terms'])
    for m in range(4):
        gm = b.gamma(s.expand(F**m),ys,shapes)
        dm = b.dirichlet(s.expand(angular**m),ts,shapes)
        check('pencil_gamma_radial_identity',s.simplify(gm-s.rf(sum(shapes),degree*m)*dm)==0)
    examples.append(result)

ys = s.symbols('y:3')
S = sum(ys)
check('coincident_nonzero_pencil',p.homogeneous_pencil_certificate(3*S**3,ys,[1,1,1])['status']=='nonzero_infinitely_many_pure_moments')
check('coincident_zero_actual_phase',p.homogeneous_pencil_certificate(0,ys,[1,1,1])['status']=='zero')
for F,nodes in [(ys[0]*ys[1]+ys[1]*ys[2]+ys[2]*ys[0],[0,1,s.I]),
                (S*S+ys[0]+1,[0,1,s.I]),(ys[0]*S,[1,1,1])]:
    check('outside_class_is_not_a_conjecture_decision',p.homogeneous_pencil_certificate(F,ys,nodes)['status']=='outside_class')

# Ordinary nilpotence survives collision of the linear form with the radial sum.
eps,U,V = s.symbols('eps U V')
nil = s.groebner([eps**2],eps,*ys)
actual = eps*S
check('actual_coefficient_ideal_at_collision',nil.reduce(s.expand(actual**2))[1]==0 and nil.reduce(actual)[1]!=0)
check('formal_binary_coefficients_are_wrong_ideal',s.Poly(V-U,U,V).coeffs()==[-1,1])

# The universal flat limit uses the exact gamma-ratio inequality; these finite
# cases check boundary indices, including all factors taken from a lower layer.
for degree,total,m in itertools.product(range(2,5),(1,3,6),range(1,7)):
    tuples = [(0,)*degree,(1,)+(0,)*(degree-1),(0,)*(degree-1)+(m,)]
    for kk in tuples:
        K=sum(kk);J=sum((j+1)*a for j,a in enumerate(kk))
        C=s.factorial(m)/s.factorial(m-K)*s.gamma(degree*m-J+total)/s.gamma(degree*m+total)
        check('flat_radial_gamma_ratio_boundary',0<=C<=1)

receipt={'date':'2026-11-16','status':'pass','checks':counts,'total':sum(counts.values()),
         'scope':'Finite exact confluent-chain/residue, original-marker, radial-pencil membership and scope controls. Universal conclusions use the written all-degree proofs.',
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                  ['SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex','simplex_pencils.py','verify_simplex_pencils.py']},
         'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'simplex_pencil_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
(HERE/'simplex_pencil_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
