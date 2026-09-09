"""Load explicit certified seed curves without launching arithmetic campaigns."""
import hashlib
import json
from pathlib import Path

import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT/'elliptic-curves/data/research_curve_supplement.json'


def load_supplement(existing):
    claims = {r['id']: r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
    rows = []
    # Bucket by exact j, then test Q-isomorphism: equal-j twists stay distinct.
    buckets = {}
    def bucket(model):
        inv = cert.weierstrass_invariants(model)
        return inv['c4']**3/inv['discriminant']
    for row in existing:
        buckets.setdefault(bucket(row['ainvs']), []).append(row)
    for entry in json.loads(MANIFEST.read_text())['curves']:
        claim = claims[entry['claim']]
        if claim['state'] != 'proved':
            raise ArithmeticError('supplement requires a proved ledger entry')
        documents = {}
        for key in ('packet', 'proof'):
            raw = (ROOT/entry[key]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != entry['sha256'][key]:
                raise ArithmeticError('supplement source hash differs: '+entry[key])
            documents[key] = json.loads(raw)
        packet, proof = documents['packet'], documents['proof']
        if entry['kind'] == 'progression':
            if proof['status'] != 'PASS_INDEPENDENT_ORACLE_FREE_RANK18_PROGRESSION':
                raise ArithmeticError('progression replay did not pass')
            if packet['n'] != '0' or packet['parameter'] != proof['base_parameter_at_u_zero']:
                raise ArithmeticError('progression seed parameter differs')
            if packet['inputs'][entry['proof']] != entry['sha256']['proof']:
                raise ArithmeticError('progression proof binding differs')
            model = packet['elliptic_a_invariants']
            points = packet['generic_points'] + [packet['seed']]
            if any(len(p) != 3 or p[2] != '1' for p in points):
                raise ArithmeticError('progression points must use affine projective chart z=1')
            points = [p[:2] for p in points]
            rank = proof['rank_at_zero']
        else:
            if proof['status'] not in ('PASS_STANDALONE_FUNNEL_M18',
                    'PASS_STANDALONE_SMALL_CONIC_M18', 'PASS_STANDALONE_FUNNEL_EPOCH'):
                raise ArithmeticError('standalone replay did not pass')
            if packet['parameter'] != proof['parameter']:
                raise ArithmeticError('seed parameter differs')
            if entry['kind'] == 'seed' and entry['sha256']['packet'] not in proof['inputs'].values():
                raise ArithmeticError('standalone proof does not bind seed packet')
            model, points, rank = packet['curve'], packet['points'], proof['rank_lower_bound']
            if proof['matrix_rank'] != rank:
                raise ArithmeticError('rank differs from certified matrix rank')
        if len(points) != rank or any(not cert.is_on_weierstrass_curve(model, p) for p in points):
            raise ArithmeticError('supplement point identities differ')
        peers = buckets.setdefault(bucket(model), [])
        duplicate = next((r for r in peers if cert.isomorphic(model, r['ainvs'])), None)
        if duplicate:
            if rank > duplicate['rank_lower_bound']:
                raise ArithmeticError('duplicate has stronger rank; explicitly select its proof')
            if duplicate not in rows:
                raise ArithmeticError('supplement duplicates base inventory; reconcile provenance')
            duplicate.setdefault('aliases', []).append(entry['id'])
            continue
        row = {k: None for k in ('conductor', 'log_conductor', 'conductor_divisor',
            'conductor_upper_bound', 'bad_primes', 'discriminant', 'log_abs_discriminant',
            'naive_height', 'faltings_height', 'original_to_minimal_isomorphism', 'conductor_certificate')}
        row.update(id=entry['id'], ainvs=model, points=points,
            rank_lower_bound=rank, local_search_rank_lower_bound=rank,
            rank_provenance='standalone_exact_finite_group_certificate',
            family='det1092', parameter=packet['parameter'],
            parameter_coordinate='original t' if entry['kind']=='progression' else 'reduced s',
            icarm_ids=[], conductor_status='UNKNOWN', known_bad_primes=[],
            rank_source_certificate=entry['packet'], rank_certificate=entry['proof'],
            rank_certificate_model='Exactly the displayed source model; no minimal-model transport asserted.',
            model_status='SOURCE_MODEL_MINIMALITY_NOT_CERTIFIED',
            canonical_source=claim['canonical_source'], status_claim=entry['claim'],
            source_sha256=entry['sha256'])
        rows.append(row)
        peers.append(row)
    return rows


def curve_page(row):
    """Keep uncomputed minimal-model metrics separate from raw source equations."""
    r = row
    lines = [f'# {r["id"]}', '',
        f'Rank **≥ {r["rank_lower_bound"]}**. Family `{r["family"]}` at {r["parameter_coordinate"]} `{r["parameter"]}`.', '',
        f'[Full data and certified points]({r["id"]}.json) · [Inventory](../../INVENTORY.md)', '',
        'Source a-invariants (global minimality has not been certified for this inventory):', '',
        '```text', ', '.join(r['ainvs']), '```', '',
        f'Conductor status: **{r["conductor_status"]}**. Minimal discriminant, naive height and Faltings height are uncomputed in this inventory.', '',
        f'Certified point count: {len(r["points"])}. This is a rank lower bound; bounded no-gain searches do not prove exact rank.', '',
        f'[Canonical proof note](../../../{r["canonical_source"]}) · [Equation and points](../../../{r["rank_source_certificate"]}) · [Independent rank proof](../../../{r["rank_certificate"]})', '',
        'This finite seed supplement makes no literature-wide novelty or conductor-record claim.']
    if r['conductor']:
        lines += ['', 'Exact conductor:', '', '```text', r['conductor'], '```', '',
                  'Complete bad primes:', '', '```text', ', '.join(r['bad_primes']), '```']
    elif r['conductor_divisor']:
        lines += ['', 'Certified conductor divisor:', '', '```text', r['conductor_divisor'], '```', '',
                  'Certified conductor upper bound:', '', '```text', r['conductor_upper_bound'], '```', '',
                  'Proved bad primes (incomplete):', '', '```text', ', '.join(r['known_bad_primes']), '```']
    if r['conductor_certificate']:
        lines += ['', f'[Conductor certificate](../../../{r["conductor_certificate"]})']
    if r.get('aliases'):
        lines += ['', 'Duplicate seed packets on this same rational-isomorphism class: '+', '.join(f'`{a}`' for a in r['aliases'])+'.']
    return '\n'.join(lines)+'\n'
