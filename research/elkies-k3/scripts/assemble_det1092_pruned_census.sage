#!/usr/bin/env sage-python
"""Recheck complete anchor packets and classify every rootless complement.

Completeness uses the exact pruned enumeration, bound to its source and input
hashes. This independently rechecks retained lattice witnesses, not the full
enumeration itself. No partial packet can close the classification.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import runpy
from sage.all import ZZ, QQ, matrix, pari

source = Path(__file__).with_name('classify_det1092_rootless_j2_pruned.sage')
m = runpy.run_path(str(source))
old, ROOT = m['old'], m['ROOT']


def isometry(G, H):
    raw = pari(G).qfisom(pari(H))
    if raw == 0:
        return None
    U = matrix(ZZ, raw.sage())
    # PARI's returned columns express G in the basis of H.
    assert abs(U.det()) == 1 and U.transpose() * H * U == G
    U = matrix(ZZ, U.inverse())
    assert abs(U.det()) == 1 and U.transpose() * G * U == H
    return old['rows'](U)


def assemble(work):
    anchors = json.loads(old['ANCHORS'].read_text())['anchors']
    bound = m['bindings']()
    assert len(anchors) == 16
    packets = []
    for n, a in enumerate(anchors, 1):
        path = work / ('anchor-%02d.json' % n)
        d = json.loads(path.read_text())
        assert d['bindings'] == bound
        assert d['status'] == 'PASS_COMPLETE_ANCHOR'
        assert d['anchor_number'] == n and d['anchor'] == a
        assert len(d['sixths']) == len(d['seventh_chambers'])
        assert sum(x['rootless_embeddings'] for x in d['seventh_chambers']) == len(d['embeddings'])
        packets.append(d)
    catalog = json.loads(old['CATALOG'].read_text())['rooted_niemeier_lattices']
    grams = {x['label']: matrix(ZZ, x['gram']) for x in catalog}
    A = matrix(ZZ, json.loads(old['AUXILIARY'].read_text())['auxiliary']['gram'])
    known = matrix(QQ, json.loads(old['PARENT'].read_text())['generic_height_gram'])
    classes, assignments, coverage = [], [], []
    for n, d in enumerate(packets, 1):
        G = grams[d['anchor']['niemeier']]
        local_counts = Counter()
        for index, x in enumerate(d['embeddings']):
            B = matrix(ZZ, x['auxiliary_basis_in_ambient'])
            W = matrix(ZZ, x['complement_basis_in_ambient'])
            H = matrix(ZZ, x['gram'])
            assert B * G * B.transpose() == A and old['primitive'](B)
            assert W == (B * G).right_kernel_matrix()
            assert W * G * W.transpose() == H and H.det() == 1092
            assert H.nrows() == 17 and int(pari(H).qfminim(2)[0]) == 0
            match, transform = None, None
            for k, c in enumerate(classes):
                transform = isometry(matrix(ZZ, c['gram']), H)
                if transform is not None:
                    match = k
                    break
            if match is None:
                match = len(classes)
                shortest = pari(H).qfminim()
                classes.append({'class_index': match + 1, 'gram': old['rows'](H),
                    'gram_sha256': old['gram_hash'](H), 'determinant': int(H.det()),
                    'minimum': int(shortest[1]), 'minimum_vector_count': int(shortest[0]),
                    'automorphism_group_order': int(pari(H).qfauto()[0]),
                    'representative_embedding': dict(x, anchor_number=n, niemeier=d['anchor']['niemeier']),
                    'primitive_embeddings_in_weyl_cover': 0})
                transform = old['rows'](matrix.identity(ZZ, H.nrows()))
            classes[match]['primitive_embeddings_in_weyl_cover'] += 1
            local_counts[str(match + 1)] += 1
            assignments.append({'anchor_number': n, 'embedding_index': index,
                'class_index': match + 1, 'integral_isometry_from_representative': transform})
        coverage.append({'anchor_number': n, 'niemeier': d['anchor']['niemeier'],
            'd5_anchor_index': d['anchor']['anchor_index'],
            'sixth_count': len(d['sixths']), 'sixth_accounting': d['sixth_accounting'],
            'seventh_chambers': d['seventh_chambers'],
            'primitive_rootless_embeddings': len(d['embeddings']), 'class_counts': dict(local_counts)})
    matching = []
    for c in classes:
        transform = isometry(known, matrix(ZZ, c['gram']))
        c['matches_recovered_curve302_frame'] = transform is not None
        c['known_frame_integral_isometry'] = transform
        if transform is not None:
            matching.append(c['class_index'])
    attained = len(matching) == 1
    return {'schema': 'elkies-k3.det1092-pruned-rootless-j2-census.v1',
        'status': 'PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION' if attained else 'CONTRADICTION_KNOWN_FRAME_NOT_RECOVERED',
        'bindings': dict(bound, **{str(Path(__file__).relative_to(ROOT)): old['digest'](Path(__file__))}),
        'anchor_packet_hashes': {('anchor-%02d.json' % n): old['digest'](work / ('anchor-%02d.json' % n)) for n in range(1, 17)},
        'accounting': {'complete_anchor_count': len(coverage),
            'primitive_rootless_embeddings_in_weyl_cover': len(assignments),
            'rootless_integral_isometry_class_count': len(classes),
            'known_frame_matching_class_indices': matching},
        'anchor_weyl_coverage': coverage, 'rootless_classes': classes, 'embedding_class_assignments': assignments,
        'scope': {'J2_frame_classification': 'COMPLETE' if attained else 'DEBUG_REQUIRED',
            'J1_surface_automorphism_classification': 'NOT_COMPUTED',
            'rational_marked_nef_U_realizations': 'NOT_COMPUTED', 'elliptic_equations': 'NOT_COMPUTED',
            'enumeration_replay': 'Hash-bound exact pruned enumeration; retained witnesses rechecked here.',
            'embedding_count': 'Representatives in the declared Weyl cover; not the number of all embeddings or all ambient-isometry orbits.',
            'excluded_ambient_cases': 'Leech has no D5 roots; the complete supplied D5-anchor classification excludes rooted Niemeier types without D5 embeddings.'}}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, default=m['DEFAULT_WORK'])
    p.add_argument('--output', type=Path, default=ROOT/'artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json')
    args = p.parse_args()
    result = assemble(args.work)
    raw = old['canonical'](result)
    if args.output.exists():
        assert args.output.read_bytes() == raw
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('xb') as f:
            f.write(raw)
    print(result['status'], result['accounting'], flush=True)
    if result['status'] != 'PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION':
        raise SystemExit(2)
