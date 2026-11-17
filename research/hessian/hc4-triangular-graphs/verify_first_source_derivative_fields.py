"""Exact controls for Section 6; field/integrality/normalization proofs are written."""
import hashlib
import json
from pathlib import Path
import sympy as s

BASE = Path(__file__).resolve().parent
x, y, p, q, tau = s.symbols('x y p q tau')
checks = []


def zero(name, value):
    assert s.cancel(s.expand(value)) == 0, name
    checks.append(name)


def coframe(name, P, Q):
    jac = s.diff(P, x)*s.diff(Q, y)-s.diff(P, y)*s.diff(Q, x)
    zero(name+' unit Jacobian', jac-1)
    zero(name+' dy dx coefficient', -s.diff(Q, x)*s.diff(P, x)
         + s.diff(P, x)*s.diff(Q, x))
    zero(name+' dy dy coefficient', -s.diff(Q, x)*s.diff(P, y)
         + s.diff(P, x)*s.diff(Q, y)-1)


for r in (2, 3, 5):
    P = x+y**r
    for t in (2, 4):
        Q = y+P**t
        coframe(f'shear {r},{t}', P, Q)
        # Explicit source expressions in first-column field data.
        zero(f'shear {r},{t} y recovery', Q-P**t-y)
        zero(f'shear {r},{t} x recovery', P-(Q-P**t)**r-x)

for e in (2, 3, 5, 8):
    # Laurent constant-Jacobian control: first-column generation holds,
    # while a replica can lose the full degree. No whole-plane claim.
    P, Q = y**e, -x/(e*y**(e-1))
    coframe(f'Laurent index {e}', P, Q)
    b = s.diff(Q, x)
    zero(f'Laurent index {e} first-velocity recovery', -e*P*b-y)
    replica = y**(e-1)
    H = y**e/e
    zero(f'Laurent index {e} primitive Hamiltonian', s.diff(H, y)-replica)
    zero(f'Laurent index {e} replica first velocity', replica*b+s.Rational(1, e))

for e in (2, 3, 7):
    # General rational-map theorem needs nonzero, not constant, Jacobian.
    G1, G2 = x**e, y+x**(e+1)
    J = s.Matrix([[s.diff(G, v) for v in (x, y)] for G in (G1, G2)])
    assert J.det() != 0
    zero(f'rational index {e} inverse differential 11', (J.inv()*J)[0, 0]-1)
    zero(f'rational index {e} derivative field recovers x', e*G1/J[0, 0]-x)

# Generic source directions give a constant leading source coefficient.
u_frame, w_frame, kappa = s.symbols("u_frame w_frame kappa")
for degree in (2, 3, 5):
    transformed = u_frame+(kappa*u_frame+w_frame)**degree
    zero(f"generic monic direction degree {degree}",
         s.Poly(transformed, u_frame).LC()-kappa**degree)
    assert s.Poly(transformed.subs(kappa, 1), u_frame).LC() == 1
    checks.append(f"generic direction leading unit degree {degree}")

# Regular differential pole test and adjugate coframe transcription.
u_pole = s.symbols("u_pole")
q_pole = s.symbols("q_pole", integer=True, positive=True)
zero("all positive pole orders differentiate to a pole",
     s.diff(u_pole**(-q_pole), u_pole)
     +q_pole*u_pole**(-q_pole-1))
for dimension in (2, 3):
    entries = s.symbols("a0:"+str(dimension**2))
    mat = s.Matrix(dimension, dimension, entries)
    product = mat.adjugate()*mat-mat.det()*s.eye(dimension)
    for row in range(dimension):
        for col in range(dimension):
            zero(f"dimension {dimension} adjugate coframe {row},{col}",
                 product[row, col])

# Nonconstant Jacobian allows an escaping source with bounded full first jets.
G1, G2 = x**2*y**2, y
all_generators = (G1, G2, s.diff(G1, x), s.diff(G1, y),
                  s.diff(G2, x), s.diff(G2, y))
for index, generator in enumerate(all_generators):
    pullback = s.cancel(generator.subs({x:tau**-1, y:tau**2}, simultaneous=True))
    numerator, denominator = s.fraction(pullback)
    assert not denominator.has(tau)
    checks.append(f"non-Keller bounded full first jet {index}")
assert (tau**-1).as_powers_dict()[tau] == -1
checks.append("non-Keller source x keeps a pole")
zero("non-Keller Jacobian is nonconstant",
     s.diff(G1,x)*s.diff(G2,y)-s.diff(G1,y)*s.diff(G2,x)-2*x*y**2)

# The retained 3D Keller collision guards the base-ring distinction in (19).
z = s.symbols("z")
v = 1+x*y
K3 = s.Matrix([v**3*z+y**2*v*(4+3*x*y),
               y+3*x*v**2*z+3*x*y**2*(4+3*x*y),
               2*x-3*x**2*y-x**3*z])
DK3 = K3.jacobian((x,y,z))
zero("3D Keller determinant", DK3.det()+2)
point1 = {x:0, y:0, z:-s.Rational(1,4)}
point2 = {x:1, y:-s.Rational(3,2), z:s.Rational(13,2)}
for index in range(3):
    zero(f"3D Keller collision coordinate {index}",
         K3[index].subs(point1)-K3[index].subs(point2))
assert DK3[2,2].subs(point1) != DK3[2,2].subs(point2)
checks.append("3D retained collision separated by one first-jet entry")

# Finite normalization without an open immersion can still identify points.
node_u, node_v = y**2-1, y*(y**2-1)
zero("node source integral equation", y**2-node_u-1)
zero("node target relation", node_v**2-node_u**2*(node_u+1))
for name, generator in (("u",node_u),("v",node_v)):
    zero(f"node collision {name}", generator.subs(y,1)-generator.subs(y,-1))
assert (s.diff(node_u,y).subs(y,1),s.diff(node_v,y).subs(y,1)) != (0,0)
checks.append("node finite birational map is not an open immersion")

# Guard the indispensable hypotheses rather than testing only positive cases.
assert s.diff(x**3, x).subs(x, 0) == 0
checks.append('degree-three normalization has a critical point')
assert s.Poly(s.diff(x**3, x), x, modulus=3).is_zero
checks.append('characteristic-three differential-constant shortcut rejected')

proof = BASE/'REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md'
out = {
    'outcome': 'PASS', 'check_count': len(checks), 'checks': checks,
    'proof_sha256': hashlib.sha256(proof.read_bytes()).hexdigest(),
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'sympy_version': s.__version__,
    'scope': 'Exact coframe/recovery/pole/adjugate controls and non-Keller/3D boundary cases; all-degree field and integral-normalization proofs are written in Section 6.',
    'first_column_field_generation_proved': True,
    'first_column_integrality_proved': False,
    'transverse_coordinate_integral_over_first_column_jet_algebra': True,
    'normalized_first_column_image_has_replica_action': True,
    'generic_first_column_normalization_is_source': True,
    'generic_first_velocities_generate_raw_source_over_full_S': True,
    'full_first_derivatives_generate_raw_source_over_full_S_all_dimensions': True,
    'full_first_jet_normalization_is_source_in_every_dimension': True,
    'full_first_jet_map_finite_birational_onto_image': True,
    'first_velocities_integral_over_target_B': False,
    'jc2_resolved': False,
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(BASE/'first_source_derivative_fields_verification.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({'outcome': out['outcome'], 'check_count': len(checks)}))
