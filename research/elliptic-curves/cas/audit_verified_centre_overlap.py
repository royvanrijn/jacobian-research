#!/usr/bin/env python3
"""Compare initial selected centres of two replayed runs, up to elliptic negation.

This audits selection overlap only, not point-search coverage or quartic boxes.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def audit(old, new):
    bindings = {}

    def read(path):
        raw = path.read_bytes()
        bindings[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    def load(folder):
        terminal = read(folder / 'terminal.json')
        verified = read(folder / 'verified.json')
        if (verified['status'] != 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY'
                or verified['terminal_sha256'] != bindings[str((folder / 'terminal.json').resolve())]):
            raise ArithmeticError('sealed independent replay required')
        stage = terminal['stages'][0]
        selection_path = folder / f"epoch-{stage['epoch']:02d}/landscape/selection.json"
        selection = read(selection_path)
        if stage['selection_sha256'] != bindings[str(selection_path.resolve())]:
            raise ArithmeticError('initial selection seal differs')
        return tuple(map(F, terminal['curve'])), selection

    model, left = load(old)
    other_model, right = load(new)
    if model != other_model or left['basis'] != right['basis']:
        raise ArithmeticError('same curve and initial basis required')
    a1, a2, a3, a4, a6 = model

    def keys(selection):
        result = set()
        for row in selection['centres']:
            x, y = map(F, row['point'])
            if y*y + a1*x*y + a3*y != x*x*x + a2*x*x + a4*x + a6:
                raise ArithmeticError('selected centre is off curve')
            result.add((x, min(y, -y-a1*x-a3)))
        return result

    left_keys, right_keys = keys(left), keys(right)
    return {
        'status': 'PASS_EXACT_SELECTED_CENTRE_OVERLAP',
        'old_selected_centres_up_to_sign': len(left_keys),
        'new_selected_centres_up_to_sign': len(right_keys),
        'overlap_up_to_sign': len(left_keys & right_keys),
        'new_centres_outside_old_selection': len(right_keys - left_keys),
        'bindings': bindings,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'claim_boundary': 'Exact comparison of two selected centre sets only. '
            'Not a comparison of every historical centre, completed quartic boxes, '
            'or returned point directions; no rank or absence claim.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('old', 'new', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.old.resolve(), args.new.resolve())
    if args.output.exists():
        if json.loads(args.output.read_text()) != result:
            raise ArithmeticError('retained overlap audit differs')
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['overlap_up_to_sign'])
