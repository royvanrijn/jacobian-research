"""Check the literal linear Frobenius certificate with stdlib Fraction only."""
from fractions import Fraction
import json
from pathlib import Path

data = json.loads((Path(__file__).resolve().parent/'jet-certificate.json').read_text())
rows = [[Fraction(*q) for q in row] for row in data['rows']]
weights = [Fraction(*q) for q in data['weights']]
target = [Fraction(*q) for q in data['target']]
assert len(rows) == len(weights) == 18
assert all(len(row) == len(target) == 28 for row in rows)
assert [sum(c*row[j] for c, row in zip(weights, rows)) for j in range(len(target))] == target
assert {data['variables'][j]: value for j, value in enumerate(target) if value} == {
    'K322': Fraction(1), 'K313': Fraction(-1)}
print('FROBENIUS_LITERAL_CERTIFICATE_PASS')
