"""Outcome-blind contracts for the eight-fibre, same-parent V3 pilot.

Only fixed discovery scores, height strata, equations and generic sections may
select inputs. No old point-search ledger, exceptional point or catalogue rank
is read. Numerical and resource failures are NOT negative rank observations.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import re

from v3_warm_support import atomic, bindings, curve_tuple, point_tuple, read, require, sha

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
ART = ROOT/'artifacts/generated-results/elliptic-curves'
D = LOCAL/'det1092-v3-eight-pilot'
AUTO = LOCAL/'det1092-v3-eight-controller'
CALIBRATION = LOCAL/'adaptive-visibility-cascade-v3'
SELECTION = LOCAL/'det1092-record-scale-selection-v1/selection-result.json'
PACKED_SELECTION = ART/'det1092_record_scale_intake_v1/selection-result.json.gz'
REDUCED_PARENT = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
ORBITS = ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
LATTICE = ART/'curve302_parent_degree2_multisection_lattice_v1.json'
# From the retained intake manifest and Git objects at d440079b. A different
# local experiment must not be silently substituted for the original 48.
SELECTION_SHA256 = '1dbce6d7bc287176afeae198a60fe9cae0cc8c1db9c42ffe97f71cdeaeb6b869'
PARENT_BLOB = 'ba867908f343c74c4ad21768e46c89d101c123f5'
ORBITS_BLOB = '9d1e6d7de520b96c866c11e05a672750320f9b8f'
LATTICE_BLOB = 'c35936478de14bcc5a9a431db65c138f0ce38944'
DOMAIN = 'det1092-v3-eight-pilot-v1'
FIELDS = ('id','family','parameter','model','height_bin','stratum','score_units',
          'j_numerator_bits','j_denominator_bits')
RESOURCE = {'prepare': 3600, 'preflight': 180, 'search': 7200, 'replay': 7200,
            'rss_bytes': 3221225472}
CLAIM = ('Prospective V3 execution on eight previously selected determinant-1092 '
         'fibres, starting only from 17 generic sections. A finite no-gain result '
         'is not an upper rank bound, saturation proof or absence-of-jumps proof. '
         'Same-parent transfer does not establish arbitrary-parent sensitivity. '
         'No automatic catalogue submission, retuning or 40-fibre expansion.')


def own_sources():
    names = ('det1092_v3_contract.py', 'det1092_v3_worker.py', 'det1092_v3_replay.py',
             'run_det1092_v3_trial.py', 'v3_warm_support.py', 'v3_warm_engine.py',
             'v3_warm_replay.py', 'v3_transfer_contract.py',
             'certify_compact_r17_candidates.py', 'research_runtime/supervisor.py')
    return {str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}


def git_blob(path):
    path = Path(path)
    h = hashlib.sha1(b'blob '+str(path.stat().st_size).encode()+b'\0')
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def selection_input():
    path = SELECTION if SELECTION.is_file() else PACKED_SELECTION
    raw = path.read_bytes()
    if path.suffix == '.gz':
        raw = gzip.decompress(raw)
    require(hashlib.sha256(raw).hexdigest() == SELECTION_SHA256,
            'the frozen 48-fibre discovery selection differs; do not substitute/refill')
    return json.loads(raw), path


def poly(values, t):
    value = F(0)
    for c in reversed(values):
        value = value*t + F(c)
    return value


def evaluate(record, t):
    denominator = poly(record['denominator'], t)
    require(denominator != 0, 'generic section has a pole; no replacement candidate is authorized')
    return poly(record['numerator'], t)/denominator


def j_invariant(model):
    _, _, _, a, b = curve_tuple(model)
    disc = 4*a**3+27*b*b
    require(disc != 0, 'singular specialization')
    return 6912*a**3/disc


def on_curve(model, point):
    _, _, _, a, b = curve_tuple(model)
    x, y = point
    return y*y == x*x*x+a*x+b


def specialize(parent, row):
    """Exact reduced-parent evaluation; preserve the ordered generic marking."""
    t = F(row['parameter']); q = t.denominator
    aa = [evaluate(r, t) for r in parent['a_invariants']]
    require(len(aa) == 5 and aa[:3] == [0, 0, 0], 'reduced generic model is not short')
    model = (F(0), F(0), F(0), aa[3]*q**8, aa[4]*q**12)
    require(model == curve_tuple(row['model']), 'selected equation and generic specialization disagree')
    points = tuple((evaluate(x, t)*q**4, evaluate(y, t)*q**6)
                   for x, y in parent['basis_weierstrass_coordinates'])
    require(len(points) == 17 and all(on_curve(model, p) for p in points),
            'generic seed point membership failed')
    j_invariant(model)
    return model, points


def score_key(row):
    t = F(row['parameter'])
    return (-row['score_units'], t.denominator, t.numerator, row['id'])


def hash_key(row):
    return (hashlib.sha256((DOMAIN+':'+row['id']).encode()).hexdigest(), row['id'])


def choose(selection):
    """Four strongest; one central moderate and one SHA remainder per band.

    Missing/duplicate/malformed inputs fail, not refill. Projection to FIELDS
    makes extra annotations incapable of steering the selection.
    """
    require(selection.get('status') == 'PASS', 'unfinished original selection')
    original = selection['selected']
    require(len(original) == 48, 'expected exactly 48 frozen fibres')
    rows = []
    for raw in original:
        row = {k: raw[k] for k in FIELDS}
        require(re.fullmatch(r'scale-\d{7}', row['id']) is not None, 'invalid case id')
        require(row['family'] == 'det1092-reduced' and type(row['score_units']) is int,
                'invalid family/discovery score')
        require(row['height_bin'] in (10, 11) and row['stratum'] in ('strong','moderate','lower_fixed'),
                'invalid original stratum')
        j = j_invariant(row['model'])
        require(abs(j.numerator).bit_length() == row['j_numerator_bits'] and
                j.denominator.bit_length() == row['j_denominator_bits'] and
                row['j_numerator_bits']//64 == row['height_bin'], 'arithmetic j-height changed')
        rows.append(row)
    require(len({r['id'] for r in rows}) == len({F(r['parameter']) for r in rows}) ==
            len({j_invariant(r['model']) for r in rows}) == 48, 'duplicate fibre/id/j-invariant')
    for band in (10, 11):
        counts = {s: sum(r['height_bin'] == band and r['stratum'] == s for r in rows)
                  for s in ('strong','moderate','lower_fixed')}
        require(counts == {'strong':16, 'moderate':4, 'lower_fixed':4}, 'original 16/4/4 strata changed')
    pilot = [dict(r, pilot_role='strong') for r in sorted(rows, key=score_key)[:4]]
    used = {r['id'] for r in pilot}
    for band in (10, 11):
        moderate = sorted((r for r in rows if r['stratum'] == 'moderate' and
                           r['height_bin'] == band), key=score_key)
        row = moderate[len(moderate)//2]
        require(row['id'] not in used, 'strong/moderate overlap; frozen strata inconsistent')
        pilot.append(dict(row, pilot_role='middle')); used.add(row['id'])
    for band in (10, 11):
        row = min((r for r in rows if r['height_bin'] == band and r['id'] not in used), key=hash_key)
        pilot.append(dict(row, pilot_role='sha-control')); used.add(row['id'])
    reserve = sorted((r for r in rows if r['id'] not in used), key=hash_key)
    require(len(pilot) == 8 and len(reserve) == 40, 'pilot size changed')
    return pilot, reserve


def trial_policy(original):
    """No shortlist, height, CVP or chart-order tuning. Only target32/metadata differ."""
    for key in ('anchors_per_shell','canonical_per_shell','exact_cvp_node_limit',
                'height','seconds_per_chart','gp_sha256','max_epochs','max_charts'):
        require(key in original, 'missing frozen V3 key: '+key)
    p = dict(original)
    p.update(schema='det1092-v3-eight-job.v1', initial_rank=17, generic_rank=17,
             family='det1092', generic_determinant=1092, target_rank=32,
             calibration_only=False, oracle_boundary=CLAIM, scope=CLAIM)
    return p


def validate_roster(root=ROOT, directory=D):
    r = read(directory/'roster.json')
    require(r.get('status') == 'READY_EIGHT_FIBRES' and len(r['pilot']) == 8 and len(r['reserve']) == 40,
            'trial preparation is incomplete')
    selection = read(directory/'frozen-selection.json')
    pilot, reserve = choose(selection)
    require(r['pilot'] == pilot and r['reserve'] == reserve, 'frozen pilot/reserve order differs')
    bindings(root, r['inputs']); bindings(root, r['implementation_sources'])
    for row in pilot:
        folder = directory/row['id']; entry = r['jobs'][row['id']]
        require(sha(folder/'protocol.json') == entry, 'job protocol changed')
        p = read(folder/'protocol.json')
        bindings(root, p['inputs']); bindings(root, p['sources'])
        require(p['implementation_sources'] == r['implementation_sources'], 'worker sources differ')
        require(p['initial_rank'] == 17 and p['target_rank'] == 32, 'trial rank targets changed')
    return r



def preflight_word(path, rank):
    """The complete 1092 table includes orbit zero; never use it as a centre."""
    require(17 <= rank <= 32, 'preflight rank outside declared pilot')
    with Path(path).open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if int(row['minimum_norm']) not in (8, 10):
                continue
            word = list(map(int, row['parent_MW17_w'].split()))
            require(len(word) == 17 and any(word), 'invalid nonzero generic shell witness')
            return word + [0]*(rank-17)
    raise ValueError('no nonzero norm-8/10 preflight centre in the bound orbit table')



def result_summary(results):
    """No censoring or upper-rank inference from missing/unchanged lower bounds."""
    clean = len(results) == 8 and all(r['stop_reason'] in
             ('COMPLETE_FINITE_NO_GAIN','TARGET_LOWER_BOUND_REACHED') for r in results)
    return {'completed': len(results), 'all_eight_uncensored_terminal': clean,
            'fibres_with_certified_gain': sum(r['gain'] > 0 for r in results),
            'expansion_eligible': clean and any(r['gain'] > 0 for r in results),
            'automatic_expansion': False, 'claim_boundary': CLAIM}
