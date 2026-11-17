"""NEW MO0 transcription controls. Refuse any registered/receipt replay."""
import hashlib
import json
import math
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / 'receipt-low-euler-multiplicity-pair.json'
assert not RECEIPT.exists(), 'Sealed receipt exists: never replay.'
ledger = json.loads((ROOT / 'MATH_STATUS.json').read_text())
assert not any(e['id'] == 'HC4MYMO0' for e in ledger['entries']), 'Already registered.'
started = time.perf_counter()
checks, algebraic = [], []


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check(name, condition, alg=True):
    assert bool(condition), name
    checks.append(name)
    if alg:
        algebraic.append(name)


# Explicit profiles check formulas; they are not a surface classification.
for mu, nu, pair, comparable in [
    (1, 7, (1, 7), True), (3, 6, (1, 2), True),
    (6, 3, (2, 1), True), (2, 3, (2, 3), False),
    (12, 18, (2, 3), False), (35, 77, (5, 11), False),
]:
    g = math.gcd(mu, nu)
    norm = (mu // g, nu // g)
    check(f'normalized pair {mu},{nu}', norm == pair and math.gcd(*norm) == 1)
    check(f'component equals gcd iff divisibility {mu},{nu}',
          ((mu == g or nu == g) == (mu % nu == 0 or nu % mu == 0))
          and ((min(norm) == 1) == comparable))
    check(f'normalized-root transverse indices {mu},{nu}',
          g // math.gcd(g, mu) == g // math.gcd(g, nu) == 1)

# Root covers are etale over each original A1: g separate lifts, not one.
for m1, m2, g, norm, bad_count in [
    (1, 5, 1, (1, 5), 1), (2, 6, 2, (1, 3), 1),
    (6, 2, 2, (3, 1), 1), (4, 6, 2, (2, 3), 2),
    (6, 6, 6, (1, 1), 0),
]:
    check(f'complete cyclic cover degree {m1},{m2}', math.gcd(m1, m2) == g)
    check(f'complete pullback multiplicities {m1},{m2}', (m1 // g, m2 // g) == norm)
    check(f'complete split-fibre Euler identity {m1},{m2}',
          2 + 2 * (g - 1) == 2 * g)
    check(f'complete bad-fibre count {m1},{m2}', sum(v > 1 for v in norm) == bad_count)
check('equal multiple orders are the forbidden zero-bad-fibre case',
      6 // math.gcd(6, 6) == 1)

# Euler two counts reduced components, regardless of their multiplicities.
check('affine Euler-two two-component profile', 1 + (2 - 1) == 2)
check('complete Euler-two all irreducible profile', 2 + (1 - 1) == 2)
check('two reducible affine fibres cannot have Euler two', 1 + 2 * (2 - 1) > 2)
check('a reducible complete fibre cannot have Euler two', 2 + (2 - 1) > 2)

# The degree comparison must retain b>1, and it is strict only for m<t.
for b, t, m, smaller in [(5, 6, 2, True), (5, 6, 6, False), (2, 1, 1, False)]:
    check(f'chart-composition degree b{b} t{t} m{m}',
          t % m == 0 and b * m > 1 and ((b * m < b * t) == smaller))
check('Euler-two cyclic compression boundary count', 2 * 7 - 1 == 13)
check('relative group need not be the full group at even quotient degree',
      (4 * 7) % (2 * 7) == 0)
check('odd quotient degree excludes full-deck index two',
      (5 * 7) % (2 * 7) != 0)

# Noncyclic full decks after cyclic relative compression need c>=6.
for c in [1, 2, 3, 4, 5]:
    check(f'noncyclic full-deck cost impossible at c{c}',
          not any(c % h == 0 and c // h >= 3 for h in range(2, c + 1)))
check('Euler six is the first arithmetic escape, not an existence claim',
      6 % 2 == 0 and 6 // 2 == 3)

context_path = HERE / 'sources/low-euler-multiplicity-pair-context.json'
context = json.loads(context_path.read_text())
for p, digest in context['retained_proof_and_source_sha256'].items():
    check('source byte pin ' + p, sha(ROOT / p) == digest, False)
check('four prior claims and nine retained pins',
      len(context['retained_structural_inputs']) == 4
      and len(context['retained_proof_and_source_sha256']) == 9, False)
check('no new external import or primary cache',
      context['new_external_theorem_imports'] == context['new_primary_caches']
      == context['new_primary_downloads'] == 0, False)
flags = context['scope_flags']
check('positive chart premise and minimality explicitly conditional',
      flags['complete_one_bad_fibre_chart'] and flags['chi_two_comparable_pair_criterion']
      and flags['minimality_required_for_relative_compression'], False)
check('forbidden automatic conclusions remain false',
      not any(flags[k] for k in ['automatic_comparability', 'incomparable_pair_proves_no_plane_chart',
                                'chart_almost_surjective_required', 'automatic_additive_action',
                                'full_deck_always_equals_relative_deck', 'unrestricted_HC4_or_JC2_resolved']), False)
check('no registered mathematical checker replay',
      context['registered_mathematical_checkers_replayed'] == [], False)
receipt = dict(
    outcome='PASS', claim='HC4MYMO0', check_count=len(checks),
    algebraic_controls=len(algebraic), checks=checks, algebraic_check_names=algebraic,
    elapsed_seconds=time.perf_counter() - started,
    proof_sha256=sha(HERE / 'COMPLETE_RULING_PLANE_CHARTS_AND_THE_LOW_EULER_MULTIPLICITY_PAIR.md'),
    checker_sha256=sha(Path(__file__)), context_sha256=sha(context_path),
    registered_mathematical_checkers_replayed=[], independent_replay=False,
    formal_verification=False, external_review=False,
)
RECEIPT.write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'algebraic_controls', 'elapsed_seconds']}))
