#!/usr/bin/env sage-python
"""Frozen four-basis generic17-to-specialized31 comparison, derived from V3.
Original algorithms unchanged except exact orbit/metric transport, generic-only
seeds, arm paths, finite group table reuse and explicit resource accounting.
"""
import argparse
import os
import resource
import traceback
from visibility_selection_v3 import cheap_shortlist, final_shortlist, chart_profile
import contextlib
import csv
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from importlib.machinery import SourceFileLoader
import numpy as np
from sage.all import ZZ, matrix, pari
from research_runtime.store import checkpoint
from visibility_lattice_v2 import ExactParity

CAS = Path(__file__).resolve().parent
v1 = SourceFileLoader('unchanged_cascade_v1', str(CAS/'adaptive_visibility_cascade.sage')).load_module()
ROOT, ART, ORBITS = v1.ROOT, v1.ART, v1.ORBITS
OLD = v1.D
D = ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1'/os.environ['CASCADE_ARM']
ORBITS = D/'orbits.tsv'
v1.ORBITS = ORBITS
v1.D = D
read, sha, rel, load = v1.read, v1.sha, v1.rel, v1.load


def sources():
    return {**v1.sources(), **{rel(p):sha(p) for p in
            (Path(__file__), CAS/'visibility_lattice_v2.py', CAS/'visibility_selection_v3.py')}}



def fingerprints(model, basis):
    import audit_recorded_point_mod2_rank_v3 as mod2
    from mod2_reduction_independence import _primes_up_to
    from research_runtime.finite_reduction import ReductionCache
    from research_runtime.memory_store import MemoryFactStore
    cache = ReductionCache(MemoryFactStore())
    columns = [0]*len(basis)
    offset = 0
    pivots = {}
    for prime in _primes_up_to(1000):
        if prime == 2:
            continue
        try:
            sig = mod2.signature(cache, model, basis, prime)
        except ValueError:
            continue
        for row in sig.rows:
            mod2.insert(pivots, row)
            for j, bit in enumerate(row):
                columns[j] |= int(bit) << offset
            offset += 1
    if len(pivots) != len(basis):
        raise ArithmeticError('finite coset fingerprints are not injective')
    return columns


def fp(word, columns):
    value = 0
    for x, col in zip(word, columns):
        if int(x) % 2:
            value ^= col
    return value


def point_key(model, basis, word):
    from half_lattice_pointed_sieve import linear_combination
    point = linear_combination(model, basis, word)
    if point is None:
        raise ArithmeticError('zero centre')
    x, y = point
    sign = -1 if y < 0 else 1
    return (x, abs(y)), [int(sign*v) for v in word]


def landscape(model, basis, tested, wd, policy):
    geo = load('prospective_half_lattice_v3.sage')
    from sage.all import block_diagonal_matrix, identity_matrix
    fixture = read(D/'fixture.json')
    transport = block_diagonal_matrix(matrix(ZZ, fixture['basis_transform']), identity_matrix(ZZ,len(basis)-17))
    reference = tuple(tuple(map(F,q)) for q in fixture['reference_points']) + basis[17:]
    hg, asym = geo.canonical_height_gram(model, reference)
    reference_gram = matrix(ZZ, geo.rounded_gram(hg, 1000000))
    g = transport*reference_gram*transport.transpose()
    checkpoint(wd/'metric-transport.json', {'reference_basis':[list(map(str,q)) for q in reference], 'reference_gram':[list(map(int,r)) for r in reference_gram.rows()], 'transport':[list(map(int,r)) for r in transport.rows()]})
    u = matrix(ZZ, pari(g).qflllgram()).transpose()
    inv = u.inverse()
    if inv.denominator() != 1:
        raise ArithmeticError('LLL transport not unimodular')
    reduced = u*g*u.transpose()
    exact = ExactParity(reduced.rows())
    columns = fingerprints(model, basis)
    rank = len(basis)
    buckets = {8:[], 10:[]}
    with ORBITS.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            shell = int(row['minimum_norm'])
            if shell not in buckets:
                continue
            word = list(map(int, row['parent_MW17_w'].split())) + [0]*(rank-17)
            buckets[shell].append({'orbit':int(row['orbit_mask']), 'shell':shell,
                'representative':word, 'fingerprint':fp(word, columns)})
    gg = np.asarray(g.rows(), dtype=np.int64)
    chosen, anchors = [], []
    for shell, rows in buckets.items():
        words = np.asarray([r['representative'] for r in rows], dtype=np.int64)
        if rank**2*int(abs(words).max())**2*int(abs(gg).max()) >= 2**62:
            raise ArithmeticError('canonical norm overflow')
        norms = np.einsum('ij,jk,ik->i', words, gg, words, optimize=True)
        for row, norm in zip(rows, norms):
            row['metric_norm'] = int(norm)
        rows.sort(key=lambda r:(-r['metric_norm'], r['fingerprint']))
        anchors += rows[:policy['anchors_per_shell']]
        count = 0
        for row in rows:
            key, word = point_key(model, basis, row['representative'])
            if tuple(map(str, key)) in tested:
                continue
            chosen.append({**row, 'representative':word, 'point':list(map(str,key)), 'lane':'canonical'})
            count += 1
            if count == policy['canonical_per_shell']:
                break
    count = 1 << (rank-17)
    extensions = np.asarray([[(m>>i)&1 for i in range(rank-17)] for m in range(count)], dtype=np.int64).reshape(count, rank-17)
    extension_fp = [fp(row, columns[17:]) for row in extensions]
    inverse = np.asarray(inv.rows(), dtype=np.int64) % 2
    records = []
    mapper = load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    for ai, anchor in enumerate(anchors):
        residues = np.zeros((count, rank), dtype=np.int64)
        residues[:, :17] = np.asarray(anchor['representative'][:17]) % 2
        residues[:, 17:] = extensions
        rp = (residues @ inverse) % 2
        words, norms = exact.babai(rp)
        keys = [anchor['fingerprint'] ^ f for f in extension_fp]
        order = sorted(range(count), key=lambda i:(-int(norms[i]), keys[i]))
        # Reversing coordinate enumeration must leave the semantic shortlist intact.
        reverse = sorted(reversed(range(count)), key=lambda i:(-int(norms[i]), keys[i]))
        if order != reverse or len(set(keys)) != count:
            raise ArithmeticError('coordinate-order dependence or duplicate coset')
        np.savez_compressed(wd/f'anchor-{ai:02d}-full.npz', residues=residues,
                            reduced_residues=rp, babai_words=words, babai_norms=norms)
        refined = []
        for index in cheap_shortlist(norms, keys):
            proof = exact.solve(rp[index], words[index], policy['exact_cvp_node_limit'])
            possible = []
            for v in proof['minima']:
                word = list(map(int, (matrix(ZZ,1,rank,v)*u).row(0)))
                if [x%2 for x in word] != list(map(int,residues[index])):
                    raise ArithmeticError('wrong transported CVP parity')
                key, word = point_key(model, basis, word)
                possible.append((key, word))
            key, word = min(possible)
            refined.append({'orbit':anchor['orbit'], 'shell':anchor['shell'],
                'extension':index, 'fingerprint':keys[index], 'metric_norm':proof['norm'],
                'babai_norm':int(norms[index]), 'representative':word,
                'point':list(map(str,key)), 'lane':'parity', 'cvp':proof})
        refined.sort(key=lambda r:(-r['metric_norm'], r['fingerprint']))
        for row in refined:
            mapping = mapper.mapping(model,basis,row)
            row.update(chart_profile(mapping))
            row['multiplicity'] = len(row['cvp']['minima'])//2
        finalists = final_shortlist([r for r in refined if tuple(r['point']) not in tested])
        chosen += finalists
        record = {'anchor':anchor, 'extension_count':count, 'refined':refined,
                  'full_scores_sha256':sha(wd/f'anchor-{ai:02d}-full.npz')}
        checkpoint(wd/f'anchor-{ai:02d}.json', record)
        records.append(record)
        print('LANDSCAPE',rank,ai+1,'/32','masks',count,'exact nodes',sum(r['cvp']['nodes'] for r in refined),flush=True)
    chosen.sort(key=lambda r:(-r['metric_norm'],r['fingerprint'],r['lane'],tuple(map(F,r['point']))))
    unique = []
    keys = set(tested)
    for row in chosen:
        key = tuple(row['point'])
        if key not in keys:
            unique.append(row)
            keys.add(key)
    result = {'rank':rank, 'basis':[list(map(str,p)) for p in basis],
        'rounded_gram':[list(map(int,row)) for row in g.rows()],
        'LLL':[list(map(int,row)) for row in u.rows()], 'fingerprint_columns':columns,
        'height_asymmetry':str(asym), 'extensions_per_anchor':count,
        'full_cosets_scored':len(anchors)*count, 'anchors':records, 'centres':unique}
    checkpoint(wd/'selection.json', result)
    return unique


def run(start):
    v1.guard()
    p = read(D/'protocol.json')
    if p['sources'] != sources() or any(sha(ROOT/k) != v for k,v in p['inputs'].items()):
        raise ArithmeticError('frozen basis-cascade source/input changed')
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from research_runtime.search_state import raw_state
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl
    out = D/'run'
    out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists():
        print('Already terminal');return
    fixture = read(D/'fixture.json')
    model = tuple(map(F,fixture['curve']))
    basis = tuple(tuple(map(F,q)) for q in fixture['points'])
    # Cache finite group tables only; no search results or other-arm data.
    from research_runtime.finite_reduction import ReductionCache
    shared_reductions = ReductionCache(MemoryFactStore())
    mod2.ReductionCache = lambda store: shared_reductions
    cache = Cache(MemoryFactStore())
    state = raw_state(model, basis, cache=cache, prime_bound=1000)
    if state.rank != start or tuple(tuple(map(F,p)) for p in state.basis) != basis:
        raise ArithmeticError('fixed seed certificate failed')
    if not (out/'seed.json').exists():
        checkpoint(out/'seed.json', {'curve':list(map(str,model)), 'points':[list(map(str,p)) for p in basis], 'state':state.record()})
    mapper = load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    stages, tested, total = [], set(), 0
    for epoch in range(p['max_epochs']):
        wd = out/f'epoch-{epoch:02d}'
        wd.mkdir(exist_ok=True)
        basis = tuple(tuple(map(F,p)) for p in state.basis)
        before = len(basis)
        if (wd/'stage.json').exists():
            stage = read(wd/'stage.json')
            cloud = read(wd/stage['audit'])
            mod2.check(wd/stage['audit'])
            for path in sorted(wd.glob('chart-*.json')):
                tested.add(tuple(read(path)['centre']['point']))
                total += 1
            basis2 = tuple(tuple(map(F,q)) for q in cloud['independent_points'])
            state = raw_state(model,basis2,cache=cache,prime_bound=1000)
            stages.append(stage)
            if state.rank >= p['target_rank'] or state.rank == before:
                break
            continue
        centres = read(wd/'selection.json')['centres'] if (wd/'selection.json').exists() else landscape(model,basis,tested,wd,p)
        charts = []
        for j,c in enumerate(centres):
            if total >= p['max_charts']:
                break
            path = wd/f'chart-{j:03d}.json'
            if path.exists():
                chart = read(path)
            else:
                mapping = mapper.mapping(model,basis,c)
                search = PointedQuarticSearch(state=state,centre={'coefficients':c['representative']},coordinate_policy=mapping['coordinate_policy'])
                transcript, points = backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
                if backend.replay(search,mapping,transcript) != points:
                    raise ArithmeticError('exact search witness replay changed')
                chart = {'index':j,'centre':c,'mapping':mapping,'search':transcript}
                checkpoint(path,chart)
            charts.append(chart)
            tested.add(tuple(c['point']));total += 1
            snapshot = wd/f'cloud-{j:03d}.json'
            audit = wd/f'mod2-{j:03d}.json'
            if not snapshot.exists():
                checkpoint(snapshot, {'status':'INCREMENTAL_RETAINED_CLOUD','family':'det1092-basis-cascade-v1-'+os.environ['CASCADE_ARM'],
                    'parameter':'0','curve':list(map(str,model)),'charts':charts,
                    'final_state':state.record(),'rank_lower_bound':before})
            with (wd/f'audit-{j:03d}.log').open('a') as log, contextlib.redirect_stdout(log):
                if not audit.exists():
                    mod2.build(snapshot,audit,1000,sha(snapshot))
                mod2.check(audit)
            cloud = read(audit)
            print('CHART',epoch,j+1,'/',len(centres),'rank',cloud['rank_lower_bound'], 'status',chart['search']['status'],flush=True)
            if cloud['rank_lower_bound'] > before:
                with (wd/f'audit-{j:03d}.log').open('a') as log, contextlib.redirect_stdout(log):
                    if not (wd/'modl.json').exists():
                        modl.build(audit,wd/'modl.json')
                    modl.check(wd/'modl.json')
                enlarged = tuple(tuple(map(F,q)) for q in cloud['independent_points'])
                if enlarged[:before] != basis:
                    raise ArithmeticError('certified extension lost original basis prefix')
                state = raw_state(model,enlarged,cache=cache,prime_bound=1000)
                if state.rank != len(enlarged):
                    raise ArithmeticError('independent subgroup replay mismatch')
                break
        if not charts:
            checkpoint(out/'terminal.json',{'status':'BUDGET_OR_NO_NEW_FINALISTS','final_rank_lower_bound':before,'stages':stages})
            return
        stage = {'epoch':epoch,'before':before,'after':state.rank,'charts':len(charts),
                 'stale_charts_cancelled':len(centres)-len(charts) if state.rank>before else 0,
                 'audit':audit.name,'audit_sha256':sha(audit),'full_cosets_scored':read(wd/'selection.json')['full_cosets_scored']}
        checkpoint(wd/'stage.json',stage)
        stages.append(stage);checkpoint(out/'stages.json',stages)
        print('STAGE',stage,flush=True)
        if state.rank >= p['target_rank'] or state.rank == before or total >= p['max_charts']:
            break
    checkpoint(out/'terminal.json',{'status':'COMPLETE_BOUNDED_CALIBRATION','final_rank_lower_bound':state.rank,
        'stages':stages,'charts':total,'protocol_sha256':sha(D/'protocol.json'),'read_paths':sorted(v1.READS)})


if __name__ == '__main__':
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(7200,7210))
    resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3))
    status='COMPLETED'
    try:
        run(17)
    except BaseException:
        status='EXPERIMENT_FAILURE'
        checkpoint(D/'failure.json',{'starting_commit':read(D/'protocol.json')['starting_commit'],'traceback':traceback.format_exc()})
        raise
    finally:
        own=resource.getrusage(resource.RUSAGE_SELF);child=resource.getrusage(resource.RUSAGE_CHILDREN)
        checkpoint(D/'resources.json', {'starting_commit':read(D/'protocol.json')['starting_commit'],'status':status,'wall_seconds':time.monotonic()-started,'own_cpu_seconds':own.ru_utime+own.ru_stime,'child_cpu_seconds':child.ru_utime+child.ru_stime,'peak_rss_KiB':own.ru_maxrss})
