"""V3 selection on an explicitly frozen exact parent subset, not a full shell bank.

Reuses its selection rules and arithmetic helpers, with an explicit parent bank
and integer exact CVP. A reference-solver mode supports independent replay.
"""
from fractions import Fraction as F
import time
import numpy as np
from visibility_selection_v3 import cheap_shortlist, final_shortlist, chart_profile
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
from v3_warm_engine import load, CAS
from research_runtime.store import checkpoint


def landscape(model, basis, tested, wd, policy, bank, *, reference=False):
    from sage.all import ZZ, matrix, pari
    old = load('future_v3_helpers', CAS/'adaptive_visibility_cascade_v3.sage')
    geo = load('future_height_geometry', CAS/'prospective_half_lattice_v3.sage')
    mapper = load('future_factor_free_mapper', CAS/'factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    rank, generic = len(basis), policy['generic_rank']
    if (not 0 < generic <= rank or bank['dimension'] != generic
            or bank['status'] not in ('COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'COMPLETE_FROZEN_PAIRWISE_PARENT_SUBSET')
            or bank['generic_gram_scale'] != 2 or bank['shells'] != policy['scaled_shells']
            or len(bank['rows']) != policy['anchor_count']
            or not bank['rows']
            or len({row['mask'] for row in bank['rows']}) != len(bank['rows'])):
        raise ArithmeticError('wrong generic dimension or uncertified parent bank')
    if (wd/'selection.json').exists():
        raise FileExistsError('preserve sealed landscape')
    wd.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    hg, asym = geo.canonical_height_gram(model, basis)
    g = matrix(ZZ, geo.rounded_gram(hg, 1000000))
    u = matrix(ZZ, pari(g).qflllgram()).transpose()
    inv = u.inverse()
    if inv.denominator() != 1:
        raise ArithmeticError('nonunimodular specialized LLL')
    reduced = u*g*u.transpose()
    solver = (ExactParity if reference else IntegerExactParity)(reduced.rows())
    columns = old.fingerprints(model, basis)
    # Half-integers are represented exactly in binary; integer labels retain
    # the historical 302 serialization for the calibration comparison.
    shell_label = lambda q: q//2 if q % 2 == 0 else q/2
    buckets = {shell_label(q): [] for q in bank['shells']}
    for row in bank['rows']:
        w = row['word'] + [0]*(rank-generic)
        shell = shell_label(row['norm'])
        if len(w) != rank or row['norm'] != 2*shell or shell not in buckets:
            raise ArithmeticError('invalid parent word')
        buckets[shell].append({'orbit': row['mask'], 'shell': shell,
                              'representative': w, 'fingerprint': old.fp(w, columns)})
    gg = np.asarray(g.rows(), dtype=np.int64)
    chosen, anchors = [], []
    for shell, rows in buckets.items():
        words = np.asarray([r['representative'] for r in rows], dtype=np.int64)
        if rank**2*int(abs(words).max())**2*int(abs(gg).max()) >= 2**62:
            raise ArithmeticError('canonical norm overflow')
        norms = np.einsum('ij,jk,ik->i', words, gg, words, optimize=True)
        for row, norm in zip(rows, norms):
            row['metric_norm'] = int(norm)
        rows.sort(key=lambda r: (-r['metric_norm'], r['fingerprint']))
        anchors += rows[:policy['anchors_per_shell']]
        retained = 0
        for row in rows:
            key, word = old.point_key(model, basis, row['representative'])
            if tuple(map(str, key)) in tested:
                continue
            chosen.append({**row, 'representative': word, 'point': list(map(str, key)), 'lane': 'canonical'})
            retained += 1
            if retained == policy['canonical_per_shell']:
                break
    count = 1 << (rank-generic)
    extensions = np.asarray([[(m >> i) & 1 for i in range(rank-generic)]
                             for m in range(count)], dtype=np.int64).reshape(count, rank-generic)
    extension_fp = [old.fp(row, columns[generic:]) for row in extensions]
    inverse = np.asarray(inv.rows(), dtype=np.int64) % 2
    records, mapping_cache = [], {}
    preparation_seconds = time.monotonic() - started
    for ai, anchor in enumerate(anchors):
        tick = time.monotonic()
        residues = np.zeros((count, rank), dtype=np.int64)
        residues[:, :generic] = np.asarray(anchor['representative'][:generic]) % 2
        residues[:, generic:] = extensions
        rp = (residues @ inverse) % 2
        words, norms = solver.babai(rp)
        keys = [anchor['fingerprint'] ^ f for f in extension_fp]
        if len(set(keys)) != count:
            raise ArithmeticError('noninjective extension fingerprints')
        array = wd/f'anchor-{ai:02d}-full.npz'
        np.savez_compressed(array, residues=residues, reduced_residues=rp,
                            babai_words=words, babai_norms=norms)
        scoring_seconds = time.monotonic() - tick
        tick = time.monotonic()
        refined = []
        for index in cheap_shortlist(norms, keys):
            proof = solver.solve(rp[index], words[index], policy['exact_cvp_node_limit'])
            possible = []
            for v in proof['minima']:
                word = list(map(int, (matrix(ZZ, 1, rank, v)*u).row(0)))
                if [x % 2 for x in word] != list(map(int, residues[index])):
                    raise ArithmeticError('CVP transport parity differs')
                key, word = old.point_key(model, basis, word)
                possible.append((key, word))
            key, word = min(possible)
            refined.append({'orbit': anchor['orbit'], 'shell': anchor['shell'],
                'extension': index, 'fingerprint': keys[index], 'metric_norm': proof['norm'],
                'babai_norm': int(norms[index]), 'representative': word,
                'point': list(map(str, key)), 'lane': 'parity', 'cvp': proof})
        cvp_and_points_seconds = time.monotonic() - tick
        refined.sort(key=lambda r: (-r['metric_norm'], r['fingerprint']))
        tick = time.monotonic()
        maps = []
        for row in refined:
            key = tuple(row['point'])
            if key not in mapping_cache:
                mapping_cache[key] = mapper.mapping(model, basis, row)
            mapping = mapping_cache[key]
            row.update(chart_profile(mapping))
            row['multiplicity'] = len(row['cvp']['minima'])//2
            maps.append(mapping)
        chosen += final_shortlist([r for r in refined if tuple(r['point']) not in tested])
        checkpoint(wd/f'anchor-{ai:02d}-maps.json', maps)
        record = {'anchor': anchor, 'extension_count': count, 'refined': refined,
                  'full_scores_sha256': old.sha(array),
                  'maps_sha256': old.sha(wd/f'anchor-{ai:02d}-maps.json')}
        checkpoint(wd/f'anchor-{ai:02d}.json', record)
        checkpoint(wd/f'anchor-{ai:02d}-timings.json', {'scoring_seconds': scoring_seconds,
            'cvp_and_points_seconds': cvp_and_points_seconds,
            'mapping_seconds': time.monotonic() - tick})
        records.append(record)
        print('FUTURE_LANDSCAPE', rank, ai+1, '/', len(anchors), flush=True)
    chosen.sort(key=lambda r: (-r['metric_norm'], r['fingerprint'], r['lane'], tuple(map(F, r['point']))))
    unique, seen = [], set(tested)
    for row in chosen:
        if tuple(row['point']) not in seen:
            unique.append(row)
            seen.add(tuple(row['point']))
    result = {'rank': rank, 'generic_rank': generic, 'basis': [list(map(str, p)) for p in basis],
        'rounded_gram': [list(map(int, row)) for row in g.rows()],
        'LLL': [list(map(int, row)) for row in u.rows()], 'fingerprint_columns': columns,
        'height_asymmetry': str(asym), 'extensions_per_anchor': count,
        'full_cosets_scored': len(anchors)*count, 'anchors': records, 'centres': unique}
    checkpoint(wd/'selection.json', result)
    checkpoint(wd/'timings.json', {'preparation_seconds': preparation_seconds,
        'landscape_wall_seconds': time.monotonic() - started, 'reference_solver': reference})
    return result
