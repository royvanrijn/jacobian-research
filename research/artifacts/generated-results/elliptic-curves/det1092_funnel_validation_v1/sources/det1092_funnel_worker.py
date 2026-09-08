#!/usr/bin/env sage-python
"""Bounded first-seed construction and prospective M18 -> unchanged V3 adapter.

Run under the controller. Every completed chart is immutable. The seed search
has no V3 landscape; V3 only sees the independently certified 18-point packet.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import sys

import det1092_funnel as f
from v3_warm_support import atomic, bindings, curve_tuple, point_tuple, read, require, sha

ORIGINAL_PARENT = f.ART/'curve302_recovered_mw17_parent_v1.json'


def checked_protocol(folder):
    p = read(folder/'protocol.json')
    bindings(f.ROOT, p['inputs']); bindings(f.ROOT, p['sources'])
    require(sha('/usr/bin/gp') == p['v3']['gp_sha256'], 'GP changed')
    from sage.all import pari
    import sage.version
    require(p['worker_software'] == {'sage':sage.version.version, 'pari':str(pari.version()),
                                   'python':sys.version}, 'Sage worker runtime changed')
    return p


def conic_points(row):
    """Evaluate the frozen conic maps at rational preimages; no point search."""
    from sage.all import QQ, PolynomialRing, EllipticCurve
    from det1092_v3_contract import evaluate
    cover = read(f.COVER); bindings(f.ROOT, cover['inputs'])
    chart = read(f.CHART); s = F(row['parameter'])
    aa, bb, cc, dd = map(int, chart['parameter_matrix'])
    if cc*s+dd == 0:
        return []
    tau = F((aa*s+bb)/(cc*s+dd))
    U = PolynomialRing(QQ, 'u'); K = U.fraction_field()
    def dec(v):
        return K(U(v['numerator']))/U(v['denominator'])
    T = dec(cover['parametrization']['t_of_u'])
    g = T.numerator()-QQ(tau)*T.denominator()
    require(g and g.degree() <= 2, 'conic preimage equation changed')
    roots = sorted(g.roots(QQ, multiplicities=False))
    parent = read(ORIGINAL_PARENT)
    E = EllipticCurve(QQ, [QQ(evaluate(v,tau)) for v in parent['a_invariants']])
    h = QQ(cc*s.numerator+dd*s.denominator); w = QQ(chart['weierstrass_u'])
    points = []
    for u in roots:
        if not T.denominator()(u):
            continue
        maps = [dec(cover['lift'][key]) for key in ('x_of_u','y_of_u')]
        if any(not v.denominator()(u) for v in maps):
            continue
        x, y = [v(u) for v in maps]; E([x,y])
        X = (x+E.b2()/12)*h**4/w**2
        Y = (y+(E.a1()*x+E.a3())/2)*h**6/w**3
        point = F(str(X)), F(str(Y))
        require(point[1]**2 == point[0]**3+F(row['model'][3])*point[0]+F(row['model'][4]), 'conic/model transport failed')
        points.append(point)
    return sorted(set(points))


def first_seed(model, base, proof, points):
    """Stop admission at the first certified extra point; no strict-local gate."""
    from v3_warm_engine import certified_state
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
    from research_runtime.memory_store import MemoryFactStore
    from mod2_reduction_independence import _primes_up_to
    from certify_compact_r17_candidates import checked_rank
    state = certified_state(model, base, proof)
    cache = QuotientOnlyReductionCache(MemoryFactStore())
    for point in points:
        trial = state.adjoin(point, cache=cache, extra_primes=_primes_up_to(1000))
        if trial.rank == 18:
            basis = (*base, point)
            # Standalone group enumeration independently rechecks the certificate.
            exact = json.loads(f.packed(checked_rank(model, basis, trial.reductions.primes, proof['no_rational_2_torsion_prime'])))
            certified_state(model, basis, exact)
            return point, exact
    return None


def seal_seed(folder, row, model, base, point, proof, evidence):
    from certify_compact_r17_candidates import checked_rank
    from v3_warm_engine import certified_state
    points = (*base, point)
    require(len(base) == 17 and len(points) == 18, 'one-point seed boundary violated')
    independent = checked_rank(model, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
    require(json.loads(f.packed(independent)) == proof, 'independent seed certificate differs')
    certified_state(model, points, proof)
    packet = dict(schema='det1092-funnel-m18.v1', status='CERTIFIED_M18',
                  parameter=row['parameter'], original_parameter=row.get('original_parameter'),
                  curve=list(map(str,model)), points=[list(map(str,p)) for p in points],
                  initial_rank=18, generic_rank=17, evidence=evidence,
                  proof=proof, intake_sha256=sha(folder/'intake.json'),
                  claim_boundary='Independent subgroup rank18. No full-rank upper bound or novelty claim.')
    atomic(folder/'m18.json', packet, immutable=True)
    atomic(folder/'result.json', dict(status='CERTIFIED_M18', packet_sha256=sha(folder/'m18.json'),
                                    parameter=row['parameter'], rank_lower_bound=18), immutable=True)
    return packet


def prepare(folder, row):
    from det1092_v3_contract import specialize
    from certify_compact_r17_candidates import checked_rank
    folder.mkdir(parents=True, exist_ok=True)
    atomic(folder/'intake.json', row, immutable=True)
    try:
        model, base = specialize(read(f.PARENT), row)
    except ValueError as error:
        # A pole/failed finite specialization is not a proved rank loss.
        atomic(folder/'result.json', dict(status='INHERITED_RANK_UNRESOLVED', reason=str(error)), immutable=True)
        return None
    try:
        proof = json.loads(f.packed(checked_rank(model, base)))
    except (ArithmeticError, ValueError) as error:
        atomic(folder/'result.json', dict(status='INHERITED_RANK_UNRESOLVED', reason=str(error)), immutable=True)
        return None
    seed = dict(curve=list(map(str,model)), points=[list(map(str,p)) for p in base], proof=proof)
    atomic(folder/'m17.json', seed, immutable=True)
    return model, base, proof


def seed_centres(folder, row, model, base, p):
    """The retained 2048-parity/49-centre geometry, before any quartic search."""
    from sage.all import matrix, ZZ, pari
    from v3_warm_engine import load
    geometry = load('funnel_generic17_geometry',f.CAS/'prospective_half_lattice_v3.sage')
    gram, asymmetry = geometry.canonical_height_gram(model,base)
    G = matrix(ZZ,geometry.rounded_gram(gram,1000000))
    U = matrix(ZZ,pari(G).qflllgram()).transpose()
    require(abs(U.det()) == 1,'nonunimodular seed metric change')
    inverse = U.inverse(); reduced = U*G*U.transpose()
    oracle = geometry.CosetOracle(reduced.rows())
    masks, seen, counter = [],set(),0
    while len(masks) < p['seed']['parity_samples']:
        require(counter < 100*p['seed']['parity_samples'],'parity draw cap')
        h = f.digest((p['domain']+':seed:'+row['id']+':'+str(counter)).encode()); counter += 1
        mask = int(h,16) % (2**17-1)+1
        if mask not in seen:
            masks.append(mask); seen.add(mask)
    samples = []
    for mask in masks:
        residue = matrix(ZZ,1,17,[(mask >> j)&1 for j in range(17)])
        target = [int(v)%2 for v in (residue*inverse).row(0)]
        norm,rep,error = oracle.solve(target)
        word = list(map(int,(matrix(ZZ,1,17,rep)*U).row(0)))
        require(all((word[j]-(mask >> j)) % 2 == 0 for j in range(17)), 'seed parity transport differs')
        require(sum(word[j]*G[j,k]*word[k] for j in range(17) for k in range(17)) == norm,'seed norm differs')
        samples.append(dict(parity=mask,representative=word,metric_norm=int(norm),cvp_error=str(error)))
    centres = sorted(samples,key=lambda c:(-c['metric_norm'],c['parity']))[:p['seed']['max_charts']]
    record = dict(curve=list(map(str,model)),points=[list(map(str,q)) for q in base],
                  rounded_gram=[list(map(int,r)) for r in G.rows()],
                  change_of_basis=[list(map(int,r)) for r in U.rows()],
                  maximum_gram_asymmetry=str(asymmetry),samples=samples,centres=centres,
                  scope='Numerical scheduling only; exact parity and rounded-norm checks, no CVP optimality claim.')
    atomic(folder/'seed-selection.json',record,immutable=True)
    return centres


def seed_search(run, case, replay_only=False):
    from v3_warm_engine import certified_state, load
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    p = checked_protocol(run); selection = read(run/'selection.json')
    require(selection['protocol_digest'] == f.digest(f.packed(p)), 'selection protocol binding differs')
    row = next((r for r in selection['seed_inputs'] if r['id'] == case), None)
    require(row is not None, 'case not in frozen seed selection')
    folder = run/'seeds'/case
    prepared = prepare(folder, row)
    if prepared is None:
        return
    model, base, proof = prepared
    state = certified_state(model, base, proof)
    # Read persisted conic witnesses only after recomputing them from equations.
    conic = conic_points(row)
    atomic(folder/'conic-points.json', [list(map(str,p)) for p in conic], immutable=True)
    found = first_seed(model, base, proof, conic)
    if found:
        require(not list(folder.glob('chart-*.json')), 'seed charts executed after a conic seed')
        return seal_seed(folder, row, model, base, *found, evidence={'kind':'frozen_conic', 'sha256':sha(folder/'conic-points.json')})
    mapper = load('funnel_seed_mapper', f.CAS/'factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    centres = seed_centres(folder,row,model,base,p)
    censored = False
    for i,centre in enumerate(centres):
        mapping = mapper.mapping(model, base, centre)
        search = PointedQuarticSearch(state=state, centre={'coefficients':centre['representative']}, coordinate_policy=mapping['coordinate_policy'])
        path = folder/f'chart-{i:03d}.json'
        if path.exists():
            chart = read(path)
            require(chart['mapping'] == mapping and chart['index'] == i, 'seed map changed on resume')
        else:
            require(not replay_only, 'replay cannot execute missing chart')
            transcript, _ = backend.execute(search, mapping, p['seed']['height'], p['seed']['seconds_per_chart'], p['v3']['gp_sha256'])
            chart = dict(index=i, mapping=mapping, search=transcript)
            atomic(path, chart, immutable=True)
        require(chart['search']['height_bound'] == p['seed']['height'] and
                chart['search']['timeout_seconds'] == p['seed']['seconds_per_chart'] and
                chart['search']['gp_binary_sha256'] == p['v3']['gp_sha256'], 'seed budget/backend binding changed')
        points = backend.replay(search, mapping, chart['search'])
        found = first_seed(model, base, proof, points)
        print('SEED_CHART', case, i+1, 'witnesses', len(points), 'certified_extra', found is not None, flush=True)
        if found:
            require(not (folder/f'chart-{i+1:03d}.json').exists(), 'stale charts after first seed')
            return seal_seed(folder, row, model, base, *found, evidence={'kind':'bounded_chart', 'index':i, 'sha256':sha(path)})
        censored |= chart['search']['status'] != 'bounded_search_complete'
    atomic(folder/'result.json', dict(status='CENSORED' if censored else 'BOUNDED_NO_CERTIFIED_SEED',
                                    parameter=row['parameter'], charts=p['seed']['max_charts'], rank_lower_bound=17), immutable=True)


def amplifier_context(run, case):
    from det1092_v3_contract import specialize
    from v3_warm_engine import Context, certified_state, load
    from certify_compact_r17_candidates import checked_rank
    p = checked_protocol(run)
    source = run/'seeds'/case
    packet = read(source/'m18.json'); row = read(source/'intake.json')
    require(packet['status'] == 'CERTIFIED_M18' and packet['intake_sha256'] == sha(source/'intake.json'), 'unbound M18 seed')
    model, base = specialize(read(f.PARENT), row)
    points = point_tuple(packet['points']); proof = packet['proof']
    require(curve_tuple(packet['curve']) == model and points[:17] == base and len(points) == 18, 'M18 lost exact generic prefix')
    checked_rank(model, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
    state = certified_state(model, points, proof)
    folder = run/'amplifiers'/case; folder.mkdir(parents=True, exist_ok=True)
    atomic(folder/'seed-input.json', {k:packet[k] for k in ('parameter','curve','points','initial_rank','generic_rank')}, immutable=True)
    atomic(folder/'seed-proof.json', proof, immutable=True)
    if (folder/'orbits.tsv').exists():
        require(sha(folder/'orbits.tsv') == sha(f.ORBITS), 'orbit copy changed')
    else:
        from det1092_v3_worker import copy_immutable
        copy_immutable(f.ORBITS, folder/'orbits.tsv')
    engine = load('funnel_unchanged_v3', f.CAS/'adaptive_visibility_cascade_v3.sage')
    require(engine.sources() == p['v3']['sources'], 'frozen V3 numerical sources changed')
    policy = dict(p['v3'])
    policy.update(schema='det1092-prospective-m18-amplifier.v1',
                  scope='Prospective V3 from exactly M17 plus first certified P18; finite search, no upper bound.',
                  seed_packet_sha256=sha(source/'m18.json'))
    atomic(folder/'protocol.json', policy, immutable=True)
    engine.D = engine.v1.D = folder
    engine.ORBITS = engine.v1.ORBITS = folder/'orbits.tsv'
    return Context(case, folder, f.ROOT, policy, engine, model, state, p['sources'])


def amplify(run, case, preflight_only=False, replay_only=False):
    from v3_warm_engine import preflight
    from det1092_v3_worker import run_search
    ctx = amplifier_context(run, case)
    preflight(ctx)
    if preflight_only:
        atomic(ctx.folder/'preflight.json', dict(status='PASS_M18_EXACT_V3_MAP', initial_rank=18,
                    packet_sha256=ctx.policy['seed_packet_sha256'], protocol_sha256=sha(ctx.folder/'protocol.json')), immutable=True)
        return
    # All transitive bindings have been checked before the artifact read guard.
    ctx.engine.v1.guard()
    if not replay_only:
        run_search(ctx)
    from det1092_funnel_replay import replay
    replay(ctx)


def conic_control(run):
    """Formula-derived u=0 control, outside population selection and its yield."""
    from sage.all import QQ, PolynomialRing
    p = checked_protocol(run)
    require(p['profile'] == 'smoke', 'constructive control only in smoke namespace')
    cover = read(f.COVER); U = PolynomialRing(QQ, 'u')
    T = cover['parametrization']['t_of_u']
    tau = F(str(U(T['numerator'])(0)/U(T['denominator'])(0)))
    aa, bb, cc, dd = map(int, read(f.CHART)['parameter_matrix'])
    s = (bb-dd*tau)/(cc*tau-aa)
    arithmetic = f.Arithmetic(p)
    row = arithmetic.candidate(arithmetic.record(-1, s.numerator, s.denominator), 'constructive_positive_control')
    row['id'] = 'conic-control'
    folder = run/'seeds'/row['id']
    model, base, proof = prepare(folder,row)
    require(first_seed(model,base,proof,[base[0]]) is None, 'known generic point falsely promoted')
    try:
        first_seed(model,base,proof,[(base[0][0],base[0][1]+1)])
    except ValueError:
        pass
    else:
        raise ValueError('off-curve seed accepted')
    points = conic_points(row); found = first_seed(model,base,proof,points)
    require(found is not None, 'formula conic control did not certify M18')
    atomic(folder/'conic-points.json', [list(map(str,q)) for q in points], immutable=True)
    seal_seed(folder,row,model,base,*found,evidence={'kind':'separate_conic_control','sha256':sha(folder/'conic-points.json')})
    amplify(run, row['id'], preflight_only=True)
    print('PASS_CONSTRUCTIVE_M18_HANDOFF_CONTROL', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['seed','replay-seed','amplify','replay-amplifier','preflight','conic-control','tables'])
    ap.add_argument('--directory', type=Path, required=True); ap.add_argument('--case')
    a = ap.parse_args(); folder = a.directory.resolve()
    if a.action == 'conic-control':
        conic_control(folder)
    elif a.action == 'tables':
        from sage.all import pari
        p = checked_protocol(folder); arithmetic = f.Arithmetic(p); checked = 0
        for prime, table in zip(arithmetic.primes, arithmetic.tables):
            for residue, (trace, roots, _) in enumerate(table):
                av, bv = arithmetic.model(residue,1) if residue < prime else arithmetic.model(1,0)
                if trace is not None:
                    require(trace == prime+1-int(pari.ellinit([0,0,0,av % prime,bv % prime],prime).ellcard()), 'independent table trace differs')
                    from sage.all import GF, PolynomialRing
                    R = PolynomialRing(GF(prime), 'x'); x = R.gen()
                    require(roots == len((x**3+av*x+bv).roots()), 'splitting differs')
                checked += 1
        atomic(folder/'tables-replay.json', dict(status='PASS_INDEPENDENT_PARI_TABLES', entries=checked,
                                                table_digest=f.digest(f.packed(arithmetic.tables))), immutable=True)
    else:
        require(a.case, '--case required')
        if a.action in ('seed','replay-seed'):
            seed_search(folder,a.case,replay_only=a.action == 'replay-seed')
        else:
            amplify(folder,a.case,preflight_only=a.action == 'preflight',replay_only=a.action == 'replay-amplifier')


if __name__ == '__main__':
    main()
