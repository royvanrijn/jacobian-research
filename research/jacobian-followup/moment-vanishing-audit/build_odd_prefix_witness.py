"""Refine one bounded-search proposal; acceptance requires the separate checker."""
from pathlib import Path
import hashlib
import json
import numpy as np
from flint import acb, acb_mat, ctx
from verify_odd_prefix_witness import evaluate, gradient_matrix

P = Path(__file__).resolve().parent
ctx.prec = 384
exploration = json.loads((P/'odd_prefix_witness_exploration.json').read_text())
point = [acb(str(v[0]), str(v[1])) for v in exploration['candidates'][0]['parameters']]


def midpoint(z):
    return acb(z.real.mid(), z.imag.mid())


for step_number in range(5):
    values, J = evaluate(point, maximum=18)
    step = J.solve(acb_mat([[value] for value in values]))
    point = [midpoint(point[i]-step[i, 0]) for i in range(9)]
    print(json.dumps({'Newton_step': step_number+1, 'step_upper': max(float(step[i, 0].abs_upper()) for i in range(9))}), flush=True)


def encode(z):
    return [z.real.mid().str(85, radius=False), z.imag.mid().str(85, radius=False)]


_, J = evaluate(point, maximum=18)
identity = acb_mat([[int(i == j) for j in range(9)] for i in range(9)])
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
    'family': exploration['family'], 'coefficient_order': exploration['coefficient_order'],
    'equations': list(range(2, 19, 2)), 'normalization': 'M_n/(3n+1)!!',
    'center': [encode(z) for z in point], 'complex_sup_radius': '1e-40',
    'preconditioner': [[encode(inverse[i, j]) for j in range(9)] for i in range(9)],
    'claimed_eta_upper': '1e-65', 'claimed_q_upper': '1e-15',
    'normalized_M20_real_interval': ['2e-21', '3e-21'],
    'normalized_M20_imag_interval': ['1e-21', '2e-21'],
    'smooth_gradient_minor_columns': selected,
    'provenance': {'exploration_sha256': hashlib.sha256((P/'odd_prefix_witness_exploration.json').read_bytes()).hexdigest(),
                   'builder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    'status': 'Rational proposal only; acceptance requires the separate interval verifier.'}
output = P/'odd_prefix_witness_certificate.json'
assert not output.exists(), 'Preserve the existing proposal; write a new version instead'
output.write_text(json.dumps(certificate, indent=2)+'\n')
print(json.dumps({'certificate': output.name, 'center': certificate['center'], 'unverified': True}, indent=2))
