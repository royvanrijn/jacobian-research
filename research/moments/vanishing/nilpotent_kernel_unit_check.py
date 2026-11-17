"""Direct sparse rational verification of the exported kernel-chart unit.

The exported coefficients have thousands of decimal digits. Parse them
as polynomial data, without evaluating a giant Python expression, and
multiply exactly with FLINT. No Groebner or modular inference is used.
"""
from pathlib import Path
import hashlib
import json
import re
import time

import flint

HERE = Path(__file__).resolve().parent


def parse_polynomial(text, variables, context):
    text = text.replace(' ', '').replace('**', '^')
    data = {}
    names = {name: i for i, name in enumerate(variables)}
    for chunk in re.findall(r'[+-]?[^+-]+', text):
        coefficient = flint.fmpq(-1 if chunk.startswith('-') else 1)
        chunk = chunk.lstrip('+-')
        powers = [0]*len(variables)
        for factor in chunk.split('*'):
            if re.fullmatch(r'\d+(?:/\d+)?', factor):
                coefficient *= flint.fmpq(factor)
            else:
                match = re.fullmatch(r'([A-Za-z][A-Za-z0-9_]*)(?:\^(\d+))?', factor)
                assert match and match[1] in names, factor
                powers[names[match[1]]] += int(match[2] or 1)
        key = tuple(powers)
        data[key] = data.get(key, flint.fmpq(0))+coefficient
    return context.from_dict(data)


def verify():
    start = time.monotonic()
    chart = json.loads((HERE/'nilpotent_kernel_linear_charts.json').read_text())['charts']['lowest_minus2']
    variables = chart['variables']
    context = flint.fmpq_mpoly_ctx.get(tuple(variables))
    log = (HERE/'nilpotent_kernel_lowest_minus2_unit_witness.log').read_text()
    lines = log.split('WITNESS_BEGIN\n')[1].split('WITNESS_END')[0].strip().splitlines()
    assert len(lines) == len(chart['moments']) == 7
    total = context.constant(0)
    counts = []
    for line, row in zip(lines, chart['moments']):
        multiplier = parse_polynomial(line, variables, context)
        source = parse_polynomial(row['primitive'], variables, context)
        total += multiplier*source
        counts.append(len(multiplier))
    assert total == context.constant(1)
    receipt = {
        'status': 'PASS', 'chart': 'lowest_minus2', 'prefix': 10,
        'method': 'Direct exact rational polynomial multiplication with FLINT; no standard-basis oracle.',
        'multiplier_nonzero_term_counts': counts,
        'witness_sha256': hashlib.sha256((HERE/'nilpotent_kernel_lowest_minus2_unit_witness.log').read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds': round(time.monotonic()-start, 3),
    }
    (HERE/'nilpotent_kernel_lowest_minus2_unit_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
