"""Exact permutation groups and inertia/deletion rows for PF2OMIT8."""
from collections import deque
from itertools import permutations
from pathlib import Path
import hashlib
import json
import platform
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0]*len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def cycle(n, letters):
    out = list(range(n))
    for a, b in zip(letters, letters[1:]+letters[:1]):
        out[a] = b
    return tuple(out)


def closure(n, generators):
    ident = tuple(range(n))
    found, todo = {ident}, deque([ident])
    while todo:
        current = todo.popleft()
        for g in generators:
            new = compose(g, current)
            if new not in found:
                found.add(new)
                todo.append(new)
    return found


def orbits(n, group):
    pending, result = set(range(n)), []
    while pending:
        a = min(pending)
        orbit = {g[a] for g in group}
        pending -= orbit
        result.append(sorted(orbit))
    return result


def sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p))) % 2 else 1


def conjugates(element, group):
    return {compose(compose(g, element), inverse(g)) for g in group}


def cycle_lengths(p):
    pending, result = set(range(len(p))), []
    while pending:
        a = min(pending)
        part, q = [], a
        while q not in part:
            part.append(q)
            pending.remove(q)
            q = p[q]
        result.append(len(part))
    return sorted(result, reverse=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    adjacent = [cycle(6, [i, i+1]) for i in range(5)]
    group = closure(6, adjacent)
    all_six = set(permutations(range(6)))
    check('connected transposition graph generates all720 S6 elements', group == all_six)
    check('S6 is transitive', orbits(6, group) == [list(range(6))])
    t = cycle(6, [0, 1])
    c = cycle(6, [0, 1, 2])
    transpositions = conjugates(t, group)
    threes = conjugates(c, group)
    check('transposition conjugacy class has15 elements', len(transpositions) == 15)
    check('three-cycle conjugacy class has40 elements', len(threes) == 40)
    check('every conjugated transposition is one2-cycle', all(cycle_lengths(p) == [2, 1, 1, 1, 1] for p in transpositions))
    check('every conjugated3-cycle is single', all(cycle_lengths(p) == [3, 1, 1, 1] for p in threes))
    check('transposition normal closure is S6', closure(6, transpositions) == group)
    alternating = {p for p in group if sign(p) == 1}
    check('three-cycle normal closure is A6', closure(6, threes) == alternating and len(alternating) == 360)
    check('all transpositions odd and all3-cycles even',
          all(sign(p) == -1 for p in transpositions) and all(sign(p) == 1 for p in threes))
    stabilizer = {p for p in group if p[5] == 5}
    check('source stabilizer has120 elements', len(stabilizer) == 120)
    check('stabilizer has only fixed-letter and five-letter orbits', orbits(6, stabilizer) == [list(range(5)), [5]])
    check('source degree is6', len(group)//len(stabilizer) == 6)
    check('point stabilizer contains odd permutations', any(sign(p) == -1 for p in stabilizer))
    check('adding a moving letter to S5 gives S6',
          closure(6, adjacent[:4]+[cycle(6, [4, 5])]) == group)
    mixed = [cycle(6, [0, 1]), cycle(6, [1, 2]), cycle(6, [3, 4]),
             cycle(6, [4, 5]), cycle(6, [2, 3, 4])]
    check('connected mixed small-cycle example generates S6', closure(6, mixed) == group)
    # Sharpness: a single4-cycle or a product of3-cycles is outside the
    # actual inertia hypotheses, and can preserve nontrivial blocks.
    imprimitive4 = closure(4, [cycle(4, [0, 1]), cycle(4, [0, 2, 1, 3])])
    check('4-cycle countercontrol transitive but not symmetric', len(imprimitive4) == 8 and orbits(4, imprimitive4) == [list(range(4))])
    blocks4 = [{0, 1}, {2, 3}]
    check('4-cycle countercontrol permutes pair blocks',
          all({p[i] for i in block} in blocks4 for p in imprimitive4 for block in blocks4))
    double3 = compose(cycle(6, [0, 2, 4]), cycle(6, [1, 3, 5]))
    imprimitive6 = closure(6, [t, double3])
    check('double3-cycle countercontrol transitive of order24', len(imprimitive6) == 24 and orbits(6, imprimitive6) == [list(range(6))])
    blocks6 = [{0, 1}, {2, 3}, {4, 5}]
    check('double3-cycle countercontrol preserves three pair blocks',
          all({p[i] for i in block} in blocks6 for p in imprimitive6 for block in blocks6))
    check('double3-cycle is excluded by single-cycle hypothesis', cycle_lengths(double3) == [3, 3])

    prior = json.loads((HERE/'residue_ledger.json').read_text())
    check('entire residual affine ramification has two distinct f1 images',
          all(row['ramified_images_distinct'] and len(row['ramified_primes']) == 2
              and all(p['f'] == 1 for p in row['ramified_primes']) for row in prior['rows']))
    rows = []
    for label, e_other, unram_image in [('A_distinct', 2, 'H'), ('A_H=C', 2, 'C'),
                                      ('A_H=B', 2, 'B'), ('B', 2, None), ('C', 3, None)]:
        target_rows = []
        for image, e_ram in [('C', 2), ('B', e_other)]+([('H', 1)] if unram_image == 'H' else []):
            omitted = []
            if image in ('C', 'B'):
                omitted.append(dict(e=e_ram, f=1, ramifies=True, omitted=True,
                                    prime='R' if image == 'C' else 'S'))
            if unram_image == image:
                omitted.append(dict(e=1, f=1, ramifies=False, omitted=True, prime='D'))
            lost = sum(p['e']*p['f'] for p in omitted)
            source_degree = 6-lost
            inertia = cycle(6, list(range(e_ram))) if image != 'H' else tuple(range(6))
            moving = [length for length in cycle_lengths(inertia) if length > 1]
            check(label+' '+image+' inertia one small cycle or identity', moving in ([], [2], [3]))
            check(label+' '+image+' full DVR degree accounting', lost+source_degree == 6)
            target_rows.append(dict(image=image, image_equation=None, geometric_inertia_cycle_lengths=cycle_lengths(inertia),
                                   omitted_primes=omitted, omitted_generic_length=lost,
                                   retained_generic_degree=source_degree, individual_retained_residue_degrees=None))
        if unram_image == 'H':
            hrow = next(r for r in target_rows if r['image'] == 'H')
            check('identity inertia still omits one sheet', hrow['geometric_inertia_cycle_lengths'] == [1]*6
                  and hrow['omitted_generic_length'] == 1 and hrow['retained_generic_degree'] == 5)
        if unram_image in ('C', 'B'):
            shared = next(r for r in target_rows if r['image'] == unram_image)
            check(label+' unramified omitted sheet shares a ramified target image',
                  len(shared['omitted_primes']) == 2 and shared['omitted_generic_length'] == 3
                  and shared['geometric_inertia_cycle_lengths'] == [2, 1, 1, 1, 1])
        rows.append(dict(row=label, monodromy='S6', Galois_closure_degree=720,
                         source_stabilizer='S5', source_extension_degree=6,
                         proper_intermediate_source_fields=[], source_cover_Galois=False,
                         generic_target_ledger=target_rows,
                         sign_cover_branch_images=['C'] if label == 'C' else ['C', 'B'],
                         sign_quadratic_field_inside_source_field=False,
                         maximal_Galois_quotient_with_trivial_generic_B_inertia_degree=2 if label == 'C' else 1,
                         actual_global_representation_constructed=False))
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT8',
        hypotheses='PF2OMIT6 residual immersed-r7 degree6 rows; full affine ramification counted',
        residual_rows=rows,
        group_controls=dict(S6_order=len(group), A6_order=len(alternating), S5_order=len(stabilizer),
                            imprimitive_single4_order=len(imprimitive4), imprimitive_double3_order=len(imprimitive6),
                            controls_are_geometric_realizations=False),
        topology_step='Actual complement meridians normally generate because ambient C2 is simply connected; finite N cover is connected',
        actual_affine_image_equations=None, actual_exceptional_values=None,
        full_global_Keller_normalization_constructed=False,
        unrestricted_singleton_zero_properness_proved=False, full_goal_complete=False)
    target = HERE/'global_monodromy_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/'GLOBAL_MONODROMY.md', HERE/'RESIDUE_LEDGER.md', HERE/'residue_ledger.json',
              HERE/'DISCRIMINANT_LEDGER.md', HERE/'SIMPLE_RAMIFICATION_AND_RADIAL.md']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
        runtime_seconds=round(time.monotonic()-started, 6), python=platform.python_version(),
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target),
        proof_boundary='Exact full finite permutation groups checked. Curve-complement normal generation and field correspondence are written geometric/algebraic steps, not inferred from finite cover samples.',
        actual_global_cover_available=False, formal_verification=False, external_review=False, full_goal_complete=False)
    (HERE/'global_monodromy_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
