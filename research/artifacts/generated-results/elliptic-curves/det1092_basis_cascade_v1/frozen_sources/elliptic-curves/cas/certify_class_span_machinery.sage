#!/usr/bin/env sage-python
"""Fixed small-field controls and adaptation of the frozen MW16 proof."""
import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys
import unittest
from sage.all import RealIntervalField, pari
from sage.version import version
import class_span_grh as engine
from class_span_fixtures import CASES, fixture
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
D = ROOT/'artifacts/local/elliptic-curves/class-span-machinery-v1'
FIXTURES = ART/'class-span-v1'
OUT = ART/'class_span_machinery_v1.json'
MWOUT = ART/'class_span_mw16_regression_v1.json'
SOURCE = Path(__file__).resolve()


def hashed(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def prepare():
    sources = [SOURCE, Path(engine.__file__), SOURCE.with_name('class_span_fixtures.py'),
               SOURCE.with_name('verify_class_span_grh.py'),
               ROOT/'elliptic-curves/tests/test_class_span_grh.py',
               ROOT/'elliptic-curves/notes/CLASS_SPAN_GRH_MACHINERY.md']
    if (D/'protocol.json').exists():raise FileExistsError('preserve protocol')
    checkpoint(D/'protocol.json', {'schema':'number-fields.class-span-validation-protocol.v1',
        'sources':{str(p.relative_to(ROOT)):hashed(p) for p in sources},
        'cases':CASES, 'cutoff':100, 'precision_bits':192, 'workers':1,
        'wall_seconds_per_stage':600, 'rss_bytes':1610612736,
        'scope':'Seven fixed small fields in five signatures and degrees2,3,4; two proper-span negative controls; one missing-formal-relation positive control; exact malformed-input tests. Optional adaptation replays only the existing MW16 proof at3/17. No new large relation collection or curve sweep.',
        'oracle_boundary':'Small bnf computations are unconditionally certified and used only to build calibration witnesses. Production verification uses nf and exact principal ideal products, never a class-group oracle.',
        'mw16_boundary':'The frozen MW16 checker remains unchanged. Its exact rows and independent-anchor proof feed the general low-level algebra and analytic test; this adapter retains that prior proof as a dependency.'})


def protocol():
    p = read(D/'protocol.json')
    for name,h in p['sources'].items():
        if hashed(ROOT/name) != h:raise ArithmeticError('frozen source differs: '+name)
    return p


def small(check=False):
    p = protocol()
    records = []
    for name, coefficients in p['cases']:
        if check:
            document = read(FIXTURES/(name+'.input.json'))
            _, oracle = fixture(coefficients, p['cutoff'])
        else:
            document, oracle = fixture(coefficients, p['cutoff'])
        result = engine.verify_document(document)
        if result['status'] != 'CERTIFIED_UNDER_GRH' or result['class_two_rank_upper_bound_under_grh'] != oracle['unconditional_class_two_rank']:
            raise ArithmeticError('small-field bound differs from certified oracle')
        if check:
            if result != read(FIXTURES/(name+'.certificate.json')):raise ArithmeticError('fixture replay differs')
        else:
            for suffix,data in [('input',document),('certificate',result)]:
                path = FIXTURES/(name+'.'+suffix+'.json')
                if path.exists():raise FileExistsError('preserve fixture')
                checkpoint(path,data)
        records.append({'case':name,'signature':document['field']['signature'],
                        'oracle':oracle,'upper_bound_under_grh':result['class_two_rank_upper_bound_under_grh'],
                        'input_sha256':engine.digest(document),'certificate_sha256':engine.digest(result)})
    negatives = []
    for name in ['imaginary_c2','imaginary_c2_squared']:
        doc = read(FIXTURES/(name+'.input.json'));doc['anchors'] = []
        result = engine.verify_document(doc)
        if result['status']!='UNKNOWN' or result['class_two_rank_upper_bound_under_grh'] is not None:
            raise ArithmeticError('proper span incorrectly certified')
        negatives.append({'case':name,'result':result})
    doc = read(FIXTURES/'imaginary_trivial.input.json')
    doc['relations'] = [r for r in doc['relations'] if all(doc['columns'][c]['p']!=97 for c,e in r['factorization'])]
    incomplete = engine.verify_document(doc)
    if incomplete['formal_quotient_dimension']<=0 or incomplete['class_two_rank_upper_bound_under_grh']!=0:
        raise ArithmeticError('missing-formal-relation control differs')
    suite = unittest.defaultTestLoader.discover(str(ROOT/'elliptic-curves/tests'), pattern='test_class_span_grh.py')
    tests = unittest.TextTestRunner(verbosity=1).run(suite)
    if not tests.wasSuccessful():raise ArithmeticError('regression tests failed')
    result = {'schema':'number-fields.class-span-machinery-validation.v1','status':'PASS',
        'sources':{**p['sources'],str((D/'protocol.json').relative_to(ROOT)):hashed(D/'protocol.json')},
        'software':{'sage':version,'pari':str(pari('version()'))},'positive_cases':records,
        'negative_controls':negatives,'incomplete_formal_presentation_control':incomplete,
        'tests_run':tests.testsRun,
        'claim_boundary':'Validated reusable sufficient generation test for Cl(K)/2 under quadratic class-character GRH. Matching toy oracle bounds validate the implementation; production output contains no lower bound or exact elliptic rank.'}
    if check:
        if read(OUT)!=result:raise ArithmeticError('validation replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve validation')
        checkpoint(OUT,result)
    print('GENERAL CLASS-SPAN VALIDATION PASS:',tests.testsRun,'TESTS, SEVEN FIELDS',flush=True)


def mw16(check=False):
    protocol()
    source = SOURCE.with_name('certify_small_conductor_class_completion.sage')
    old = runpy.run_path(str(source))
    proof = old['expected']()
    if proof!=read(old['OUT']):raise ArithmeticError('frozen MW16 proof differs')
    runner = runpy.run_path(str(old['RUNNER']))
    m = runner['prior_state'](3)
    span = engine.RelationSpan(len(m.base['columns']))
    for row in [*m.canonical.values(), *m.rows]:span.add_relation(row)
    algebra = span.analyze([1 << c for c in proof['anchor_columns']])
    margin = engine.quadratic_margin(int(proof['field_discriminant']), [3,0],
        [c['p']**c['f'] for c in m.base['columns']], algebra['known_coordinates'],50000,160)
    old_known = {x['column']:x['anchor_bits'] for x in proof['known_prime_anchor_coordinates']}
    known = {c:b for c,b in algebra['known_coordinates'].items() if m.base['columns'][c]['p']**m.base['columns'][c]['f'] < 50000}
    if known != old_known or algebra['anchor_image_dimension_upper_bound'] != 16 or not margin['positive_margin_certified']:
        raise ArithmeticError('general machinery does not reproduce MW16 span')
    R=RealIntervalField(160)
    a,b=margin['margin_interval'],proof['independent_conservative_test']['conservative_margin_interval']
    if not R(a['lower'],a['upper']).overlaps(R(b['lower'],b['upper'])):
        raise ArithmeticError('general and frozen formula intervals disagree')
    result={'schema':'number-fields.class-span-mw16-regression.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashed(p) for p in [SOURCE,Path(engine.__file__),source,old['OUT'],D/'protocol.json']},
        'known_prime_count':len(known),'anchor_dimension':16,'analytic_test':margin,
        'class_two_rank_upper_bound_under_grh':16,'prior_unconditional_lower_bound':16,
        'prior_curve_rank_under_grh':22,
        'claim_boundary':'General algebra and formula reproduce the frozen MW16 proof after full prior replay. The old proof supplies arithmetic and independence; this is a regression, not a new curve result.'}
    if check:
        if read(MWOUT)!=result:raise ArithmeticError('MW16 adapter replay differs')
    else:
        if MWOUT.exists():raise FileExistsError('preserve MW16 regression')
        checkpoint(MWOUT,result)
    print('GENERAL MACHINERY MW16 PASS: UPPER16, MATCHING4740 MEMBERSHIPS',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','build','check','mw16','check-mw16','launch-build','launch-check','launch-mw16','launch-check-mw16'])
    args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage.startswith('launch-'):
        p=protocol();stage=args.stage[7:]
        if (D/(stage+'.supervisor.json')).exists():raise FileExistsError('preserve stage')
        result=run([sys.executable,str(SOURCE),stage],limits=Limits(p['wall_seconds_per_stage'],p['rss_bytes']),
            cwd=ROOT,log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'))
        print(stage,result['outcome'],result['returncode'],flush=True)
        if result['outcome']!='completed' or result['returncode']!=0:raise SystemExit(1)
    elif args.stage in ['build','check']:small(args.stage=='check')
    else:mw16(args.stage=='check-mw16')


if __name__=='__main__':main()
