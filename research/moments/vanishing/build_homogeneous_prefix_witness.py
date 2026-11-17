"""Refine a retained numerical candidate and emit exact rational isolating data.

The separate verifier is the acceptance gate. Midpoint Newton steps here
are only a way to propose a center and a preconditioner.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from flint import acb, acb_mat, ctx
from verify_homogeneous_prefix_witness import evaluate, gradient_matrix

P = Path(__file__).resolve().parent
ctx.prec = 320
exploration = json.loads((P/'homogeneous_prefix_witness_exploration.json').read_text())
candidate = exploration['candidates'][0]
point = [acb(str(v[0]), str(v[1])) for v in candidate['parameters']]


def midpoint(z):
    return acb(z.real.mid(), z.imag.mid())


for _ in range(5):
    values, J = evaluate(point, maximum=12)
    step = J.solve(acb_mat([[value] for value in values]))
    point = [midpoint(point[i]-step[i, 0]) for i in range(6)]


def encode(z):
    return [z.real.mid().str(75, radius=False), z.imag.mid().str(75, radius=False)]


values, J = evaluate(point, maximum=12)
identity = acb_mat([[int(i == j) for j in range(6)] for i in range(6)])
inverse = J.solve(identity)
gradient = gradient_matrix(point)
approx = np.array([[complex(gradient[i, j]) for j in range(18)] for i in range(15)])
selected = []
residual = approx.copy()
for _ in range(15):
    norm = np.linalg.norm(residual, axis=0)
    norm[selected] = -1
    column = int(np.argmax(norm))
    assert norm[column] > 1e-10
    selected.append(column)
    unit = residual[:, column]/norm[column]
    residual -= np.outer(unit, unit.conj()@residual)

certificate = {
    'family': exploration['family'], 'coefficient_order': ['b', 'c', 'd', 'e', 'f', 'h'],
    'equations': [2, 4, 6, 8, 10, 12], 'normalization': 'M_n/(3n+1)!!',
    'center': [encode(z) for z in point], 'complex_sup_radius': '1/10000000000000000000000000000000000000000',
    'preconditioner': [[encode(inverse[i, j]) for j in range(6)] for i in range(6)],
    'claimed_eta_upper': '1e-60', 'claimed_q_upper': '1e-20',
    'normalized_M14_real_interval': ['2e-11', '4e-11'],
    'normalized_M14_imag_interval': ['-9e-11', '-8e-11'],
    'smooth_gradient_minor_columns': selected,
    'provenance': {'exploration_sha256': hashlib.sha256((P/'homogeneous_prefix_witness_exploration.json').read_bytes()).hexdigest(),
                   'builder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    'status': 'Rational candidate certificate; accepted only after the separate verifier passes.'}
output = P/'homogeneous_prefix_witness_certificate.json'
assert not output.exists(), 'Preserve the existing certificate; write a new version instead'
output.write_text(json.dumps(certificate, indent=2)+'\n')
print(json.dumps({'certificate': output.name, 'center': certificate['center'], 'unverified': True}, indent=2))
