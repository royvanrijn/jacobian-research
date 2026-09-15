"""Check a complete supplied bank selection without rerunning its optimizer.

Numeric metric construction and the frozen heuristic selection are reproduced.
Every CVP minimum and tie is independently proved in that rounded metric.
Native group arithmetic checks the selected words, and the prescribed map
recipe is reproduced. Producer node counts remain execution metadata.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import time

from cancellation_basis_epoch import read, sha, need
from finite_cancellation_corpus import canonical, digest
from verify_parity_minimum import ParityMinimumVerifier


def verify_selection(packet, case, plan, folder, selection):
    import numpy as np
    from sage.all import EllipticCurve, QQ, ZZ, matrix, pari
    from v3_warm_engine import load, CAS
    from visibility_lattice_v2 import ExactParity
    from visibility_selection_v3 import cheap_shortlist, final_shortlist, chart_profile
    start = time.process_time(); components = {}; certificates = []
    selection_bytes = (folder/'selection.json').read_bytes()
    need(json.loads(selection_bytes) == selection, 'selection object differs from bound file')
    selection_hash = digest(selection_bytes)
    n = len(packet['points']); parent = case['parent_bank']; generic = parent['dimension']
    policy = {**plan['bank_policy'], 'generic_rank': generic,
              'scaled_shells': parent['shells'], 'anchor_count': len(parent['rows'])}
    need(0 < generic <= n and n-generic <= plan['maximum_extension_dimension'], 'unsupported bank dimension')
    need(parent['status'] in ('COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'COMPLETE_FROZEN_PAIRWISE_PARENT_SUBSET', 'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET'), 'unverified parent subset')
    need(parent['generic_gram_scale'] == 2 and parent['rows'] and
         len({r['mask'] for r in parent['rows']}) == len(parent['rows']), 'invalid parent subset')
    need(packet['curve'] == case['seed']['curve'] and
         packet['points'][:len(case['seed']['points'])] == case['seed']['points'], 'changed initial subgroup')
    curve = tuple(map(F, packet['curve'])); basis = [tuple(map(F, p)) for p in packet['points']]
    need(curve[:3] == (0, 0, 0), 'short-model word normalization required')
    E = EllipticCurve(QQ, list(curve)); points = [E(list(p)) for p in basis]
    old = load('certificate_bank_fingerprints', CAS/'adaptive_visibility_cascade_v3.sage')
    geo = load('certificate_bank_metric', CAS/'prospective_half_lattice_v3.sage')
    mapper = load('certificate_bank_maps', CAS/'factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    tick = time.process_time()
    hg, asymmetry = geo.canonical_height_gram(curve, basis)
    g = matrix(ZZ, geo.rounded_gram(hg, 1000000))
    u = matrix(ZZ, pari(g).qflllgram()).transpose(); inverse = u.inverse()
    need(inverse.denominator() == 1 and abs(u.det()) == 1, 'invalid lattice transport')
    reduced = u*g*u.transpose(); verifier = ParityMinimumVerifier(reduced.rows())
    # Babai is a frozen heuristic, not the independent CVP proof. Reuse its
    # exact input conversion and floating tie rule solely to bind selection.
    heuristic = ExactParity(reduced.rows())
    need(selection['rank'] == n and selection['generic_rank'] == generic and
         selection['basis'] == packet['points'] and selection['rounded_gram'] == [list(map(int, r)) for r in g.rows()] and
         selection['LLL'] == [list(map(int, r)) for r in u.rows()] and
         selection['height_asymmetry'] == str(asymmetry), 'metric or basis binding differs')
    components['metric_cpu_seconds'] = time.process_time()-tick
    tick = time.process_time(); columns = old.fingerprints(curve, basis)
    need(selection['fingerprint_columns'] == columns, 'semantic coset keys differ')
    components['fingerprint_cpu_seconds'] = time.process_time()-tick
    native_cache = {}

    def key_word(word):
        key = tuple(map(int, word))
        if key not in native_cache:
            need(len(key) == n, 'word dimension differs')
            P = sum((c*p for c, p in zip(key, points)), E(0))
            need(not P.is_zero(), 'zero anchor')
            x, y = P.xy(); sign = -1 if y < 0 else 1
            native_cache[key] = ((F(str(x)), F(str(abs(y)))), [sign*c for c in key])
        return native_cache[key]

    tick = time.process_time(); buckets = {q//2 if q % 2 == 0 else q/2: [] for q in parent['shells']}
    for row in parent['rows']:
        word = row['word']+[0]*(n-generic); shell = row['norm']//2 if row['norm'] % 2 == 0 else row['norm']/2
        need(len(word) == n and shell in buckets and row['norm'] == 2*shell, 'parent word differs')
        norm = sum(word[i]*int(g[i,j])*word[j] for i in range(n) for j in range(n))
        buckets[shell].append({'orbit': row['mask'], 'shell': shell, 'representative': word,
                              'fingerprint': old.fp(word, columns), 'metric_norm': norm})
    anchors = []; chosen = []
    for rows in buckets.values():
        rows.sort(key=lambda r: (-r['metric_norm'], r['fingerprint']))
        anchors.extend(rows[:policy['anchors_per_shell']])
        for row in rows[:policy['canonical_per_shell']]:
            point, word = key_word(row['representative'])
            chosen.append({**row, 'point': list(map(str, point)), 'representative': word, 'lane': 'canonical'})
    components['canonical_word_cpu_seconds'] = time.process_time()-tick
    count = 1 << (n-generic)
    need(selection['extensions_per_anchor'] == count and selection['full_cosets_scored'] == count*len(anchors), 'incomplete extension census')
    need(len(selection['anchors']) == len(anchors), 'anchor count differs')
    extensions = np.asarray([[(m >> i) & 1 for i in range(n-generic)] for m in range(count)], dtype=np.int64).reshape(count,n-generic)
    extension_keys = [old.fp(row, columns[generic:]) for row in extensions]
    inv = np.asarray(inverse.rows(), dtype=np.int64) % 2
    maps_by_point = {}; cvp_cpu = 0.; map_cpu = 0.; word_cpu = 0.
    for ai, (anchor, saved) in enumerate(zip(anchors, selection['anchors'])):
        need(saved['anchor'] == anchor and saved['extension_count'] == count, 'anchor order differs')
        array = folder/f'anchor-{ai:02d}-full.npz'; maps_file = folder/f'anchor-{ai:02d}-maps.json'
        need(sha(array) == saved['full_scores_sha256'] and sha(maps_file) == saved['maps_sha256'], 'discovery payload changed')
        residues = np.zeros((count,n), dtype=np.int64)
        residues[:,:generic] = np.asarray(anchor['representative'][:generic]) % 2; residues[:,generic:] = extensions
        reduced_residues = (residues @ inv) % 2
        words, norms = heuristic.babai(reduced_residues)
        with np.load(array, allow_pickle=False) as payload:
            expected = {'residues': residues, 'reduced_residues': reduced_residues, 'babai_words': words, 'babai_norms': norms}
            need(set(payload.files) == set(expected) and all(np.array_equal(payload[k], v) for k,v in expected.items()), 'full score array differs')
        keys = [anchor['fingerprint'] ^ x for x in extension_keys]
        need(len(set(keys)) == count, 'noninjective extension keys')
        indices = cheap_shortlist(norms, keys); by_index = {r['extension']:r for r in saved['refined']}
        need(len(by_index) == len(saved['refined']) and set(by_index) == set(indices), 'shortlist omitted or added a coset')
        checked = []
        for index in indices:
            row = by_index[index]; proof = row['cvp']; tick = time.process_time()
            need(proof['certificate'] == 'complete exact rational LDL ellipsoid enumeration' and
                 type(proof['nodes']) is int and 1 <= proof['nodes'] <= policy['exact_cvp_node_limit'], 'invalid discovery metadata')
            certified = verifier.verify(list(map(int,reduced_residues[index])), proof['norm'], proof['minima'], policy['exact_cvp_node_limit'])
            need(proof['nodes'] >= certified['verification_nodes'], 'discovery count below fixed-radius nodes')
            certificates.append({'anchor': ai, 'extension': index, **certified})
            cvp_cpu += time.process_time()-tick; tick = time.process_time(); possible = []
            for vector in certified['minima']:
                w = list(map(int,(matrix(ZZ,1,n,vector)*u).row(0)))
                need([x % 2 for x in w] == list(map(int,residues[index])), 'parity transport differs')
                possible.append(key_word(w))
            point, word = min(possible)
            expected = {'orbit': anchor['orbit'], 'shell': anchor['shell'], 'extension': index,
                'fingerprint': keys[index], 'metric_norm': certified['norm'], 'babai_norm': int(norms[index]),
                'representative': word, 'point': list(map(str,point)), 'lane': 'parity', 'cvp': proof}
            word_cpu += time.process_time()-tick; tick = time.process_time()
            if point not in maps_by_point:
                maps_by_point[point] = mapper.mapping(curve, basis, {'representative': word})
            expected.update(chart_profile(maps_by_point[point])); expected['multiplicity'] = len(certified['minima'])//2
            need(expected == row, 'refined word, profile or multiplicity differs')
            map_cpu += time.process_time()-tick; checked.append(expected)
        checked.sort(key=lambda r: (-r['metric_norm'], r['fingerprint']))
        need(checked == saved['refined'], 'refined order differs')
        stored_maps = read(maps_file); need(len(stored_maps) == len(checked), 'map list differs')
        for row, stored in zip(checked, stored_maps):
            point, _ = key_word(stored['centre']['representative'])
            need(list(map(str,point)) == row['point'], 'cached map has wrong centre')
            expected = maps_by_point[point]
            need({k:v for k,v in stored.items() if k != 'centre'} == {k:v for k,v in expected.items() if k != 'centre'}, 'prescribed map recipe differs')
        chosen.extend(final_shortlist(checked))
    chosen.sort(key=lambda r: (-r['metric_norm'],r['fingerprint'],r['lane'],tuple(map(F,r['point']))))
    unique = []; seen = set()
    for row in chosen:
        point = tuple(row['point'])
        if point not in seen:
            seen.add(point); unique.append(row)
    need(unique == selection['centres'], 'complete exported centre order differs')
    components.update(cvp_cpu_seconds=cvp_cpu, refined_word_cpu_seconds=word_cpu, map_recipe_cpu_seconds=map_cpu)
    return {'status':'PASS_INDEPENDENT_BANK_SELECTION_CERTIFICATE', 'rank': n,
        'seed_sha256': digest(canonical(packet)), 'parent_bank_sha256': digest(canonical(parent)),
        'selection_sha256': selection_hash, 'cvp_certificates': certificates,
        'native_word_checks': len(native_cache), 'map_recipes': len(maps_by_point),
        'complete_centres': len(unique), 'components': components, 'cpu_seconds': time.process_time()-start,
        'boundary':'Full frozen selection and prescribed maps reproduced. Independent rational fixed-radius enumeration proves all parity minima and ties in the rounded metric. Exact native group arithmetic checks all representative choices. Producer optimizer node counts remain execution metadata; canonical heights remain numerical. The caller must certify the supplied subgroup before use.'}


def rebuild(packet, case, plan, folder):
    from cancellation_scheduler_cpu import cpu
    from visibility_complement_subset import landscape
    from v3_warm_engine import certified_state
    from finite_cancellation_corpus import write
    start = cpu(); model = tuple(map(F,packet['curve'])); basis = [tuple(map(F,p)) for p in packet['points']]
    certified_state(model,basis,packet['proof']); write(folder/'seed.json',packet)
    parent = case['parent_bank']; policy = {**plan['bank_policy'], 'generic_rank':parent['dimension'],
        'scaled_shells':parent['shells'],'anchor_count':len(parent['rows'])}
    need(len(basis)-parent['dimension'] <= plan['maximum_extension_dimension'], 'extension limit exceeded')
    tick = cpu(); selection = landscape(model,basis,set(),folder/'producer',policy,parent); producer_cpu = cpu()-tick
    proof = verify_selection(packet,case,plan,folder/'producer',selection)
    write(folder/'selection-proof.json',proof)
    bank = {'seed':packet,'centres':selection['centres'][:plan['maximum_centres']],'rank':len(basis)}
    need(bank['centres'],'empty bank');write(folder/'bank.json',bank)
    receipt = {'status':'PASS_INDEPENDENT_BANK_CERTIFICATE','rank':len(basis),
        'bank_sha256':sha(folder/'bank.json'),'seed_sha256':sha(folder/'seed.json'),
        'producer_sha256':sha(folder/'producer/selection.json'), 'selection_proof_sha256':sha(folder/'selection-proof.json'),
        'parent_bank_sha256':digest(canonical(parent)), 'exported_centres':len(bank['centres']),
        'production_cpu_seconds':producer_cpu,'independent_verification_cpu_seconds':proof['cpu_seconds'],
        'cpu_seconds':cpu()-start,'boundary':proof['boundary']}
    write(folder/'verification.json',receipt)
    return bank, receipt
