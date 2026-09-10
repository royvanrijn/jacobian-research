"""Apply selected, ledger-backed curve updates without point searches."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT/'elliptic-curves/data/research_curve_refresh.json'


def apply_refresh(rows):
    manifest = json.loads(MANIFEST.read_text())
    claims = {r['id']: r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
    documents = {}
    for path, digest in manifest['sources'].items():
        raw = (ROOT/path).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ArithmeticError('curve refresh source hash differs: '+path)
        documents[path] = json.loads(raw)

    def claim_for(entry):
        claim = claims[entry['claim']]
        if claim['state'] != 'proved' or entry['source'] not in claim['software_lock']:
            raise ArithmeticError('curve refresh requires a proved, source-bound claim')
        return claim

    for entry in manifest['curves']:
        claim = claim_for(entry)
        source = documents[entry['source']]
        family, parameter = entry['family'], entry['parameter']
        if entry['kind'] == 'foundry':
            selected = source['records'][entry['record']]['result']
            receipt = selected['certificate_replay']
            packet = selected['packet']
            digest = hashlib.sha256((json.dumps(packet,indent=2,sort_keys=True)+'\n').encode()).hexdigest()
            if (selected['status'] != 'PASS_CERTIFIED_SEARCH'
                    or receipt['status'] != 'PASS_TWO_FINITE_IMPLEMENTATIONS'
                    or digest != selected['packet_sha256'] or digest != receipt['packet_sha256']
                    or selected['family'] != family or selected['parameter'] != parameter
                    or selected['rank_lower_bound'] != packet['rank_lower_bound']
                    or receipt['rank_lower_bound'] != packet['rank_lower_bound']):
                raise ArithmeticError('foundry packet identity/replay differs')
        elif entry['kind'] == 'productive':
            if source['status'] != 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY':
                raise ArithmeticError('productive replay did not pass')
            def record(name):
                return next(v for k, v in source['records'].items() if k.endswith('/'+name))
            packet, protocol, replay = (record(n) for n in ('terminal.json', 'protocol.json', 'verified.json'))
            if (protocol['family'] != family or protocol['parameter'] != parameter
                    or replay['rank_lower_bound'] != packet['rank_lower_bound']
                    or source['final_lower_bound'] != packet['rank_lower_bound']
                    or replay['status'] != source['status']):
                raise ArithmeticError('productive terminal identity differs')
        elif entry['kind'] == 'panel':
            if (source['status'] != 'COMPLETE_BOUNDED_R17_60_PANEL'
                    or source['completed'] != 60 or source['independently_verified'] != 60):
                raise ArithmeticError('complete verified panel required')
            selected = next(r for r in source['results'] if r['id'] == entry['record'])
            if (selected['family'] != family or selected['parameter'] != parameter
                    or selected['status'] != 'PASS_INDEPENDENT_R17_60_CASE'):
                raise ArithmeticError('panel identity differs')
            packet = selected['packet']
        elif entry['kind'] == 'cohort':
            selected = next(r for r in source['records'] if r['id'] == entry['record'])
            if selected['family'] != family or selected['parameter'] != parameter:
                raise ArithmeticError('cohort identity differs')
            packet = selected['packet']
        elif entry['kind'] == 'reconciled':
            if source['status'] != 'PASS_REPLAYED_CLOUD_RECONCILIATION':
                raise ArithmeticError('cloud reconciliation did not pass')
            packet = source
        elif entry['kind'] == 'bifibration':
            frame, replay = documents[entry['frame']], documents[entry['replay']]
            if (replay['status'] != 'PASS_INDEPENDENT_BIFIBRATION_CONICS_AND_FIXED_M18_SEED'
                    or source['status'] != 'NEW_INDEPENDENT_DIRECTION'
                    or source['original_parameter'] != parameter
                    or replay['seed']['original_parameter'] != parameter
                    or source['frame_sha256'] != manifest['sources'][entry['frame']]):
                raise ArithmeticError('bifibration seed identity differs')
            packet = dict(curve=frame['curve'], points=frame['basis']+[source['terminal']],
                          rank_lower_bound=replay['seed']['augmented_rank'],
                          proof=dict(signatures=frame['records'],
                                     no_rational_2_torsion_prime=frame['no_two_torsion_prime']))
        else:
            raise ArithmeticError('unknown refresh packet kind')
        model = tuple(map(F, packet['curve']))
        points = tuple(tuple(map(F, p)) for p in packet['points'])
        proof = packet['proof']
        rank = packet['rank_lower_bound']
        if rank != len(points) or rank != entry['rank_lower_bound']:
            raise ArithmeticError('selected rank differs from packet')
        # Replay only the saved finite reductions, including torsion exclusion.
        checked_rank(model, points, [r['prime'] for r in proof['signatures']],
                     proof['no_rational_2_torsion_prime'])
        existing = next((r for r in rows if r['id'] == entry['id']), None)
        if existing:
            if rank <= existing['rank_lower_bound'] or not cert.isomorphic(model, existing['ainvs']):
                raise ArithmeticError('rank update does not strengthen the same curve')
            if model == tuple(map(F,existing['ainvs'])):
                transported = points
            else:
                if not existing.get('original_to_minimal_isomorphism'):
                    raise ArithmeticError('different source model needs an explicit transport')
                u, r, s, t = map(F, existing['original_to_minimal_isomorphism'])
                transported = [((x-r)/u**2, (y-s*(x-r)-t)/u**3) for x, y in points]
            if any(not cert.is_on_weierstrass_curve(existing['ainvs'], p) for p in transported):
                raise ArithmeticError('updated point transport failed')
            row = existing
            row['points'] = [list(map(str, p)) for p in transported]
        else:
            if any(cert.isomorphic(model, r['ainvs']) for r in rows):
                raise ArithmeticError('new refresh curve duplicates inventory')
            row = {k: None for k in ('conductor', 'log_conductor', 'conductor_divisor',
                'conductor_upper_bound', 'bad_primes', 'discriminant', 'log_abs_discriminant',
                'naive_height', 'faltings_height', 'original_to_minimal_isomorphism', 'conductor_certificate')}
            row.update(id=entry['id'], ainvs=list(map(str, model)),
                points=[list(map(str, p)) for p in points], family=family, parameter=parameter,
                parameter_coordinate='original t' if family == 'det1092' else 'native parameter',
                icarm_ids=[], conductor_status='UNKNOWN', known_bad_primes=[],
                model_status='SOURCE_MODEL_MINIMALITY_NOT_CERTIFIED',
                rank_certificate_model='Exactly the displayed source model; no minimal-model transport asserted.')
            rows.append(row)
        row.update(rank_lower_bound=rank, local_search_rank_lower_bound=rank,
            rank_provenance='standalone_exact_finite_group_certificate',
            rank_source_certificate=entry['source'], rank_certificate=entry['source'],
            canonical_source=claim['canonical_source'], status_claim=entry['claim'],
            source_sha256=manifest['sources'][entry['source']])
        if existing and existing.get('original_to_minimal_isomorphism'):
            row['rank_certificate_model'] = 'Selected source packet; points transported to the certified minimal model by the saved isomorphism.'

    for entry in manifest['conductors']:
        claim_for(entry)
        row = next(r for r in rows if r['id'] == entry['id'])
        source = documents[entry['source']]
        if entry.get('kind') == 'checkpoint':
            result = source['records'][entry['id']]['result']
            if result['status'] != 'PASS_INDEPENDENT_CONDUCTOR_REPLAY':
                raise ArithmeticError('independent conductor replay required')
            audit = result['certificate']
            if not cert.isomorphic(row['ainvs'], audit['curve']):
                raise ArithmeticError('conductor equation differs')
            remaining = int(audit['remaining_cofactor'])
            divisor = math.prod(int(r['prime'])**r['conductor_exponent'] for r in audit['local_data'])
            delta = abs(cert.weierstrass_invariants(audit['curve'])['discriminant'])
            if (remaining * math.prod(int(r['prime'])**r['discriminant_valuation'] for r in audit['local_data']) != delta
                    or str(divisor) != audit['conductor_divisor']
                    or str(divisor*remaining) != audit['conductor_upper_bound']):
                raise ArithmeticError('local conductor products differ')
            previous = row.get('conductor_certificate')
            if previous and previous != entry['source']:
                row.setdefault('prior_conductor_certificates',[]).append(previous)
            lower = math.lcm(int(row.get('conductor_divisor') or 1),divisor)
            upper = min(int(row.get('conductor_upper_bound') or divisor*remaining),divisor*remaining)
            if lower > upper:raise ArithmeticError('incompatible conductor bounds')
            row.update(conductor_divisor=str(lower),conductor_upper_bound=str(upper),conductor_certificate=entry['source'])
            new_bad = [r['prime'] for r in audit['local_data'] if r['conductor_exponent']]
            row['known_bad_primes'] = sorted(set(row['known_bad_primes']) | set(new_bad),key=int)
            if audit['status'] == 'EXACT':
                if remaining != 1 or audit['exact_conductor'] != str(divisor):
                    raise ArithmeticError('incomplete exact conductor')
                if row.get('conductor') and row['conductor'] != str(divisor):
                    raise ArithmeticError('existing exact conductor differs')
                row.update(conductor=str(divisor),conductor_status='EXACT',log_conductor=math.log(divisor),bad_primes=new_bad)
            elif audit['status'] != 'UNKNOWN' or audit['exact_conductor'] is not None:
                raise ArithmeticError('invalid partial conductor status')
            continue
        selected = next(r for r in source['rows'] if r['family'] == row['family'])
        audit = selected['audit']
        if not cert.isomorphic(row['ainvs'], audit['curve']):
            raise ArithmeticError('conductor equation differs')
        row.update(conductor_divisor=audit['conductor_lower_bound'],
            conductor_upper_bound=audit['conductor_upper_bound'], conductor_certificate=entry['source'],
            known_bad_primes=[r['prime'] for r in audit['local_data'] if r['conductor_valuation']])
        if audit['status'] == 'EXACT_CONDUCTOR_CERTIFIED':
            conductor = audit['exact_conductor']
            if (source['exact_replay']['sha256'] != selected['audit_sha256']
                    or audit['unprocessed_cofactor'] != '1'
                    or audit['conductor_lower_bound'] != conductor
                    or audit['conductor_upper_bound'] != conductor
                    or math.prod(int(r['prime'])**r['conductor_valuation'] for r in audit['local_data']) != int(conductor)):
                raise ArithmeticError('exact conductor binding differs')
            row.update(conductor=conductor, conductor_status='EXACT',
                       log_conductor=math.log(int(conductor)), bad_primes=row['known_bad_primes'])

    for entry in manifest.get('metrics', []):
        claim_for(entry)
        row = next(r for r in rows if r['id'] == entry['id'])
        source = documents[entry['source']]
        if (source['status'] != 'PASS_EXACT_MINIMAL_MODEL_AND_POINT_TRANSPORT'
                or source['id'] != entry['id']
                or source['rank_lower_bound'] != row['rank_lower_bound']
                or source['source_ainvs'] != row['ainvs']
                or source['point_transports'] != len(row['points'])
                or source['conductor'] != row.get('conductor')):
            raise ArithmeticError('selected minimal-model metric binding differs')
        if any(not cert.is_on_weierstrass_curve(source['minimal_ainvs'], point)
               for point in source['points']):
            raise ArithmeticError('selected transported point does not lie on minimal model')
        row.update(
            ainvs=source['minimal_ainvs'],
            points=source['points'],
            discriminant=source['discriminant'],
            log_abs_discriminant=source['log_abs_discriminant'],
            naive_height=source['naive_height'],
            faltings_height=source['faltings_height'],
            original_to_minimal_isomorphism=source['original_to_minimal_isomorphism'],
            model_status='GLOBAL_MINIMAL_MODEL_CERTIFIED',
            minimal_model_certificate=entry['source'],
            rank_certificate_model='Selected source packet; all points transported to the certified global minimal model by the saved isomorphism.')
    return rows
