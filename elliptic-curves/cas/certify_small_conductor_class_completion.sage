#!/usr/bin/env sage-python
"""Certify the sixteen-generator span by a quadratic-character explicit formula."""
import argparse
import json
from pathlib import Path
import runpy
from sage.all import RealIntervalField, RealBallField
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT, ART, cert = batch.ROOT, batch.ART, batch.cert
SOURCE = Path(__file__).resolve()
RUNNER = ROOT/'elliptic-curves/cas/pursue_small_conductor_class_target_residual.sage'
TRIANGLE = ROOT/'elliptic-curves/cas/certify_small_conductor_smaller_base_v2.sage'
PARITY = ROOT/'elliptic-curves/cas/certify_small_conductor_selmer_rank_target.py'
NOTE = ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md'
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-class-completion-v1'
OUT = ART/'small_conductor_class_completion_v1.json'
CHARACTERS = ART/'small_conductor_class_characters_v1.json'
LOWER = ART/'small_conductor_class_lower16_v1.json'
BASE = ART/'small_conductor_norm_batch_relations_v1.json'
LAST = ART/'small_conductor_class_target_residual_wave_002_v1.json'


def record(x):
    return {'lower':str(x.lower()), 'upper':str(x.upper())}


def interval(R, x):
    return R(x['lower'], x['upper'])


def reduce_row(row, basis):
    while row:
        pivot = row.bit_length()-1
        if pivot not in basis:
            return row
        row ^= basis[pivot]
    return 0


def full_basis(rows):
    basis = {}
    for row in rows:
        row = reduce_row(row, basis)
        if row:
            basis[row.bit_length()-1] = row
    return basis


def anchor_normal_form(m, column, mask):
    row = 1 << column
    if column in m.canonical:
        row ^= m.canonical[column]
    inside, outside = row & m.inside_mask, row & ~m.inside_mask
    while outside:
        pivot = outside.bit_length()-1
        if pivot not in m.outside:
            return None
        o, i = m.outside[pivot]
        outside ^= o
        inside ^= i
    while inside & ~mask:
        pivot = (inside & ~mask).bit_length()-1
        if pivot not in m.inside:
            return None
        inside ^= m.inside[pivot]
    return inside


def conservative_character_margin(columns, known, discriminant, degree, r1, T, bits):
    """Generic sufficient test; caller certifies prime completeness and membership.

    For a quadratic character trivial on H, known prime classes have sign +1;
    unknown prime classes have worst sign -1 only at odd powers. The coarse
    archimedean upper bound is BDyDF (13), omitting its negative tail.
    """
    R = RealIntervalField(bits)
    L = R(T).log()
    signed_sum = R(0)
    terms = 0
    for c, col in enumerate(columns):
        q = col['p']**col['f']
        power, exponent = q, 1
        if q >= T:
            continue
        logq = R(q).log()
        while power < T:
            sign = 1 if c in known or exponent % 2 == 0 else -1
            signed_sum += sign*logq*(1-R(power).log()/L)/R(power).sqrt()
            terms += 1
            power *= q
            exponent += 1
    rhs = (R(discriminant).log()-degree*(R.euler_constant()+(8*R.pi()).log())
           -r1*R.pi()/2+(degree*R.pi()**2/2
           +4*r1*R(RealBallField(bits).catalan_constant()))/L)
    margin = 2*signed_sum-rhs
    return {'precision_bits':bits, 'prime_power_terms':terms,
            'worst_signed_prime_sum_interval':record(signed_sum),
            'archimedean_rhs_upper_bound_interval':record(rhs),
            'conservative_margin_interval':record(margin),
            'positive_margin_certified':bool(margin.lower()>0)}


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve protocol')
    sources = [SOURCE, RUNNER, TRIANGLE, PARITY, NOTE, CHARACTERS, LOWER, BASE, LAST,
               ART/'small_conductor_selmer_rank_target_v1.json']
    checkpoint(D/'protocol.json', {
        'schema':'elliptic-curves.small-conductor-class-completion-protocol.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sources},
        'residual_wave':2, 'cutoff':50000, 'precision_bits':256,
        'series_terms':4096, 'independent_precision_bits':160,
        'workers':1, 'wall_seconds_per_stage':600, 'rss_bytes':1610612736,
        'gate':'User goal16. The audited formal quotient is18, with16 independent anchors. Test whether its proved prime memberships already force the anchor span to generate Cl(K)/2 under GRH, allowing worst-case signs for every unresolved prime.',
        'selection_disclosure':'Exploratory interval evaluations at37638,40000,50000,75000 were positive. This certificate fixes50000 with a large margin and independently checks a more conservative archimedean bound.',
        'assumption':'GRH for nontrivial quadratic characters of the ordinary ideal class group of this field, in fact only those trivial on the anchor span.',
        'verification':'Replay all inherited arithmetic and principal relations; check prime membership with two distinct elimination orders; enumerate every prime-ideal power below the cutoff; certify positivity with MPFI intervals and two archimedean formulas.',
        'claim_boundary':'An analytic generating-set certificate may close the upper bound without reducing the old formal relation quotient. No missing principal relation or unconditional rank upper bound is invented.'})
    print('PREPARED CHARACTER COMPLETION AT50000', flush=True)


def expected():
    p = cert.read(D/'protocol.json')
    for name, h in p['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen source differs: '+name)
    runner = runpy.run_path(str(RUNNER))
    runner['audit_all'](p['residual_wave'], write=False)
    m = runner['prior_state'](p['residual_wave']+1)
    matrix_report = m.report()
    if matrix_report['quotient_dimension'] != 18:
        raise ArithmeticError('fixed formal quotient differs')
    chars, lower = cert.read(CHARACTERS), cert.read(LOWER)
    anchors = sorted(a['column'] for a in chars['anchors'])
    if len(anchors) != 16 or chars['anchor_character_rank'] != 16:
        raise ArithmeticError('sixteen independent anchors required')
    if lower['unconditional_class_two_rank_lower_bound'] != 16:
        raise ArithmeticError('lower bound differs')
    mask = sum(1 << c for c in anchors)
    # Independent full-coordinate row-space checks do not use the inside/outside
    # elimination order. Canonical rows are already normalized out of m.rows.
    basis = full_basis(m.rows)
    if len(basis) != matrix_report['all_relation_rank']:
        raise ArithmeticError('independent full matrix rank differs')
    augmented = full_basis([*m.rows, *(1 << a for a in anchors)])
    if len(augmented)-len(basis) != 16:
        raise ArithmeticError('formal anchor independence differs')
    known, unknown = [], []
    T = p['cutoff']
    if not 1 < T <= 400000:
        raise ArithmeticError('cutoff exceeds certified prime-ideal coverage')
    # The inherited audit reconstructs the entire base above all rational
    # primes <=400000 by exact maximal-order prime decomposition.
    for c, col in enumerate(m.base['columns']):
        if col['p']**col['f'] >= T:
            continue
        normal = anchor_normal_form(m, c, mask)
        row = (1 << c) ^ m.canonical.get(c, 0)
        independently_known = reduce_row(row, augmented) == 0
        if independently_known != (normal is not None):
            raise ArithmeticError('prime classification differs between elimination orders')
        if normal is None:
            unknown.append({'column':c, **col})
        else:
            if reduce_row(row ^ normal, basis):
                raise ArithmeticError('claimed anchor representation is not a proved relation')
            known.append({'column':c, 'anchor_bits':sum(((normal >> a)&1) << i for i,a in enumerate(anchors))})
    known_set = {r['column'] for r in known}
    R = RealIntervalField(p['precision_bits'])
    L = R(T).log()
    triangle = runpy.run_path(str(TRIANGLE))['triangle'](m.base, T, p['precision_bits'], p['series_terms'])
    penalty = R(0)
    odd_terms = 0
    for col in unknown:
        q = col['p']**col['f']
        power, exponent = q, 1
        logq = R(q).log()
        while power < T:
            if exponent % 2:
                penalty += 4*logq*(1-exponent*logq/L)/R(power).sqrt()
                odd_terms += 1
            power *= q
            exponent += 1
    margin = interval(R, triangle['margin_interval'])-penalty
    if not margin.lower()>0:
        raise ArithmeticError('character exclusion is not certified')
    independent = conservative_character_margin(m.base['columns'], known_set,
        int(m.base['field_discriminant']), 3, 3, T, p['independent_precision_bits'])
    if not independent['positive_margin_certified']:
        raise ArithmeticError('independent conservative bound did not certify generation')
    if triangle['prime_ideals_contributing'] != len(known)+len(unknown):
        raise ArithmeticError('incomplete prime classification')
    if triangle['prime_power_terms'] != independent['prime_power_terms']:
        raise ArithmeticError('prime-power counts disagree')
    signed = interval(R, triangle['prime_sum_interval'])-penalty/2
    if not signed.overlaps(interval(R, independent['worst_signed_prime_sum_interval'])):
        raise ArithmeticError('independent signed sums disagree')
    parity = json.loads(json.dumps(runpy.run_path(str(PARITY))['expected']()))
    if parity != cert.read(ART/'small_conductor_selmer_rank_target_v1.json'):
        raise ArithmeticError('rank implication replay differs')
    if parity['brumer_kramer_offset'] != 7 or parity['selmer_dimension_parity'] != 0:
        raise ArithmeticError('rank implication hypotheses differ')
    result = {
        'schema':'elliptic-curves.small-conductor-class-completion.v1', 'status':'PASS',
        'sources':{**p['sources'],str((D/'protocol.json').relative_to(ROOT)):cert.hashed(D/'protocol.json')},
        'field_discriminant':m.base['field_discriminant'], 'degree':3, 'real_embeddings':3,
        'residual_wave':p['residual_wave'], 'formal_relation_matrix':matrix_report,
        'anchor_columns':anchors, 'anchor_character_rank':16,
        'cutoff':T, 'known_prime_ideal_count':len(known), 'unknown_prime_ideal_count':len(unknown),
        'known_prime_anchor_coordinates':known, 'unknown_prime_ideals':unknown,
        'membership_verification':'All known representations lie in the exactly audited normalized principal-relation row space. A separate full-coordinate augmented-row-space test agrees for every prime ideal below the cutoff. Uncancelled outside coordinates remain unknown.',
        'triangle':triangle, 'unknown_odd_prime_power_terms':odd_terms,
        'worst_case_penalty_interval':record(penalty), 'corrected_margin_interval':record(margin),
        'independent_conservative_test':independent,
        'assumption':p['assumption'],
        'generation_argument':'If H is a proper subspace of Cl(K)/2, a nontrivial quadratic ideal-class character is trivial on H. At known primes its values are1; at unknown primes even powers are1 and odd powers are>=-1. Thus2*S_chi-C >= triangle_margin-4*S_unknown_odd >0. BDyDF2008 equation(3), no pole for nontrivial chi, gives2*S_chi-C=-sum_rho Phi(rho)<=0 under GRH, because the triangular test has nonnegative Fourier transform. Contradiction proves H=Cl(K)/2.',
        'references':['https://doi.org/10.1090/S0025-5718-07-02003-0',
                      'https://repositorio.uchile.cl/handle/2250/154648',
                      'https://arxiv.org/pdf/1607.02430'],
        'unconditional_class_two_rank_lower_bound':16,
        'conditional_on_grh_class_two_rank_upper_bound':16,
        'conditional_on_grh_class_two_rank':16,
        'certified_generating_quotient_dimension_under_grh':16,
        'conditional_on_grh_selmer_dimension':22,
        'unconditional_rank_lower_bound':22, 'unconditional_rank_upper_bound':None,
        'unconditional_exact_rank':None, 'conditional_on_grh_exact_rank':22,
        'rank_argument':'22 independent rational points give rank>=22. The proved Brumer-Kramer offset7 gives Sel_2<=16+7=23. Proved even2-Selmer parity gives Sel_2<=22, hence rank=Sel_2 dimension=22 under the class-character GRH hypothesis.',
        'claim_boundary':'The old formal relation quotient remains18; an independent analytic generation certificate supplies a16-dimensional generating quotient. No extra principal relations, new rational points, full class-group structure, or unconditional curve-rank upper bound are claimed.'}
    print('CHARACTER MARGIN', result['corrected_margin_interval'], flush=True)
    print('CLASS2=16; SELMER2=22; CURVE RANK=22 UNDER GRH', flush=True)
    return result


def launch(stage):
    p = cert.read(D/'protocol.json')
    out = D/(stage+'.supervisor.json')
    if out.exists():
        raise FileExistsError('preserve stage')
    r = run(['/home/royvanrijn/.local/bin/sage','-python',str(SOURCE),stage],
            limits=Limits(p['wall_seconds_per_stage'],p['rss_bytes']),cwd=ROOT,
            log_path=D/(stage+'.log'),checkpoint_path=out)
    print(stage, r['outcome'], r['returncode'], flush=True)
    if r['outcome']!='completed' or r['returncode']!=0:
        raise SystemExit(1)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','build','check','launch-build','launch-check'])
    args = parser.parse_args()
    if args.stage=='prepare':
        prepare()
    elif args.stage.startswith('launch-'):
        launch(args.stage[7:])
    else:
        result = expected()
        if args.stage=='check':
            if cert.read(OUT)!=result:
                raise ArithmeticError('completion certificate differs')
        else:
            if OUT.exists():
                raise FileExistsError('preserve certificate')
            checkpoint(OUT,result)
