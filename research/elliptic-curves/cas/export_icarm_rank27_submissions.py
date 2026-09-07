#!/usr/bin/env python3
"""Export seven manual submissions and replay exact models and point proofs."""
import argparse
import json
from pathlib import Path
import zipfile

import certify_compact_r17_candidates as cert
import memory_rank_certificate as memory

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'icarm_rank27_submissions_v1'
INTAKE = ROOT / 'artifacts/local/elliptic-curves/icarm-rank27-submission-v1'
INDEX = ART / 'new_high_rank_curve_index_v22.json'
AUDIT = ART / 'icarm626_publication_audit_v1.json'
MODELS = [
    'paired_high_rank_minimal_proofs_v2.json',
    'retention_high_rank_minimal_proofs_v1.json',
    'next24_high_rank_minimal_proofs_v4.json',
    'full11952_high_rank_models_v1.json',
    'new_mw16_rank27_minimal_proof_v1.json',
]
NOTES = {
    40: 'EXTENDED_PRIME_RANK27_DISCOVERIES_2026-09-06.md',
    41: 'EXTENDED_PRIME_RANK27_DISCOVERIES_2026-09-06.md',
    71: 'WIDER_RETENTION_DISCOVERIES_2026-09-06.md',
    72: 'WIDER_RETENTION_DISCOVERIES_2026-09-06.md',
    48: 'NEXT24_RANK27_DISCOVERY_2026-09-06.md',
    186: 'FULL11952_NEW_RANK27_2026-09-06.md',
    90: 'NEW_MW16_RANK27_2026-09-06.md',
}


def dump(value):
    return json.dumps(value, indent=2, sort_keys=True) + '\n'


def expected_files():
    inventory = cert.read(INDEX)
    selected = cert.read(AUDIT)['optional_rank27_submission_shortlist']
    database, intake = map(cert.read, (INTAKE / 'database.json', INTAKE / 'metadata.json'))
    if cert.hashed(INTAKE / 'database.json') != intake['sha256']:
        raise ArithmeticError('catalogue intake hash differs')
    if database['count'] != len(database['curves']):
        raise ArithmeticError('incomplete catalogue')
    byj = {}
    for q in database['curves']:
        inv = cert.weierstrass_invariants(tuple(map(cert.F, q['ainvs'])))
        byj.setdefault(inv['c4']**3 / inv['discriminant'], []).append(q)
    models = [(name, r) for name in MODELS for r in cert.read(ART / name)['curves']]
    files, evidence, sections = {}, [], []
    paths = [Path(__file__).resolve(), Path(cert.__file__), Path(memory.__file__),
             ROOT / 'elliptic-curves/cas/elliptic_candidate_record.py',
             INDEX, AUDIT, INTAKE / 'database.json', INTAKE / 'metadata.json',
             *(ART / n for n in MODELS)]
    for identifier in selected:
        row = next(r for r in inventory['curves'] if r['id'] == identifier)
        number = int(identifier.rsplit('-', 1)[1])
        source_name, minimal = next((n, r) for n, r in models
            if (r['family'], r['parameter']) == (row['family'], row['parameter']))
        model = tuple(map(cert.F, minimal.get('minimal_curve', minimal.get('integral_curve'))))
        points = [tuple(map(cert.F, p)) for p in minimal['points']]
        short = tuple(map(cert.F, row['curve']))
        if len(points) != 27 or row['rank_lower_bound'] != 27:
            raise ArithmeticError('expected 27-point local result')
        if not all(a.denominator == 1 for a in model):
            raise ArithmeticError('nonintegral submission model')
        if not all(cert.is_on_weierstrass_curve(model, p) for p in points):
            raise ArithmeticError('submitted point off curve')
        inv = cert.weierstrass_invariants(model)
        short_inv = cert.weierstrass_invariants(short)
        if any(inv[k] != short_inv[k] for k in ('c4', 'c6', 'discriminant')):
            raise ArithmeticError('model invariants differ')
        a1, a2, a3, _, _ = model
        b2 = a1*a1 + 4*a2
        transported = [(x+b2/12, y+(a1*x+a3)/2) for x, y in points]
        if transported != [tuple(map(cert.F, p)) for p in row['points']]:
            raise ArithmeticError('transport differs from certified ordered basis')
        proof = row['rank_certificate']
        actual = memory.checked_rank(short, transported,
            [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
        if dump(actual) != dump(proof):
            raise ArithmeticError('finite independence replay differs')
        matches = [q['id'] for q in byj.get(inv['c4']**3/inv['discriminant'], [])
                   if cert.isomorphic(model, q['ainvs'])]
        if matches:
            raise ArithmeticError(f'{identifier} is already published: {matches}')
        note = 'elliptic-curves/notes/' + NOTES[number]
        commentary = (
            f'Found by Roy van Rijn with AI assistance from OpenAI Codex, as a rational '
            f'specialization at t = {row["parameter"]} of the compact elliptic-curve '
            f'family labelled {row["family"]} in our research repository. '
            f'This is a specialization of an existing family; the family label identifies '
            f'the construction and is not a claim to have originated the family.\n\n'
            f'The 27 submitted rational points are exactly certified independent, proving '
            f'rank E(Q) >= 27. The equation is a proved global minimal Weierstrass model. '
            f'The exact rank and conductor remain unknown; no record is claimed.\n\n'
            f'Repository reference: {note}; stable curve ID {identifier}. '
            f'Exact point and model certificates: artifacts/generated-results/elliptic-curves/{source_name} '
            f'and {row["source_certificate"]}. The coefficients and point coordinates '
            f'were exported from these certificates and checked by exact arithmetic. '
            f'Catalogue comparisons acknowledge ICARM, supported by NSF Grant DMS 2425401.')
        payload = {'ainvs': list(map(str, model)),
                   'points': [list(map(str, p)) for p in points], 'commentary': commentary}
        prefix = identifier + '/'
        files[prefix + 'ainvs.txt'] = ', '.join(payload['ainvs']) + '\n'
        files[prefix + 'points.txt'] = '\n'.join(', '.join(p) for p in payload['points']) + '\n'
        files[prefix + 'commentary.txt'] = commentary + '\n'
        files[prefix + 'submission.json'] = dump(payload)
        md = (f'# {identifier}: rank at least 27\n\n'
              f'Family `{row["family"]}`, parameter `{row["parameter"]}`.\n\n'
              f'## a-invariants\n\n```text\n{files[prefix+"ainvs.txt"]}```\n\n'
              f'## Points: all 27 lines\n\n```text\n{files[prefix+"points.txt"]}```\n\n'
              f'## Primes of bad reduction\n\nLeave this field blank. The complete list is unknown.\n\n'
              f'## Commentary\n\n```text\n{commentary}\n```\n')
        files[prefix + 'SUBMIT.md'] = md
        sections.append(md)
        evidence.append({'id': identifier, 'family': row['family'], 'parameter': row['parameter'],
            'model_source': source_name, 'source_curve_id': minimal['id'],
            'ainvs': payload['ainvs'], 'points': payload['points'],
            'short_model': row['curve'], 'short_points': row['points'],
            'replayed_rank_certificate_on_short_model': actual,
            'minimality_source_certificate': minimal['minimality'],
            'transport': 'X=x+(a1^2+4*a2)/12; Y=y+(a1*x+a3)/2; scale u=1',
            'catalogue_matches': matches, 'exact_conductor': 'UNKNOWN'})
        paths.append(ART / row['source_certificate'])
        print(identifier, 'PASS: 27 independent points and exact integral transport', flush=True)
    summary = ('# Seven ICARM submissions\n\n'
        'Use https://elliptic-rank.icarm.cloud/ and submit one curve at a time. '
        'Each linked sheet contains all three fields to paste: a-invariants, '
        '27 point lines, and commentary. Leave primes of bad reduction blank for every curve. '
        'All models are globally minimal; all ranks are lower bounds, not exact ranks.\n\n'
        f'No Q-isomorphism matches in the {database["count"]}-curve catalogue downloaded '
        f'{intake["retrieved_at"]}. No submissions were sent by this exporter.\n\n'
        '| Submission sheet | Family | Parameter | Rank lower bound |\n|---|---|---|---:|\n' +
        ''.join(f'| [{r["id"]}]({r["id"]}/SUBMIT.md) | {r["family"]} | '
                f'{r["parameter"]} | 27 |\n' for r in evidence) +
        '\nEach directory also contains plain-text fields and a submission.json API payload '
        '(with the optional primes key omitted). All sheets are combined in ALL_SUBMISSIONS.md. '
        'verification.json retains the exact transports, rank proofs and source hashes.\n\n'
        'Replay from the repository root:\n\n```sh\n'
        'python3 elliptic-curves/cas/export_icarm_rank27_submissions.py --check\n```\n')
    files['README.md'] = summary
    files['ALL_SUBMISSIONS.md'] = summary + '\n---\n\n' + '\n---\n\n'.join(sections)
    files['verification.json'] = dump({'status': 'PASS', 'catalogue_intake': intake,
        'catalogue_count': database['count'], 'curves': evidence,
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in sorted(set(paths))},
        'claim_boundary': 'Replays the seven existing rank proofs and exact integral transports. '
            'Minimality certificates are retained from existing proofs. No new search, '
            'factorization, rank gain, record, or external submission.'})
    return files


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = expected_files()
    archive = OUT.with_suffix('.zip')
    if args.check:
        if {str(p.relative_to(OUT)) for p in OUT.rglob('*') if p.is_file()} != set(files):
            raise ArithmeticError('unexpected package file set')
        for name, content in files.items():
            if (OUT / name).read_text() != content:
                raise ArithmeticError('export differs: ' + name)
        with zipfile.ZipFile(archive) as z:
            if set(z.namelist()) != set(files):
                raise ArithmeticError('archive file set differs')
            for name, content in files.items():
                if z.read(name) != content.encode():
                    raise ArithmeticError('archive differs: ' + name)
    else:
        if OUT.exists() or archive.exists():
            raise FileExistsError('preserve existing submission package')
        OUT.mkdir()
        with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
            for name, content in sorted(files.items()):
                path = OUT / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
                info = zipfile.ZipInfo(name, date_time=(2026, 9, 7, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(info, content.encode())
    print('SEVEN MANUAL SUBMISSION PACKETS PASS')
