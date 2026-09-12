#!/usr/bin/env sage-python
"""Exact fixed-carrier trace and branch-overlap analysis, no new atlas search."""
import importlib.machinery
import importlib.util
from pathlib import Path
from sage.all import QQ

path = Path(__file__).with_name('principal28_ancestry.sage')
loader = importlib.machinery.SourceFileLoader('principal28', str(path))
spec = importlib.util.spec_from_loader(loader.name, loader)
a = importlib.util.module_from_spec(spec)
loader.exec_module(a)
g, R, OUT = a.g, a.R, a.OUT


def main():
    geometry = g.read(OUT/'geometry.json')
    replay = g.read(OUT/'verified.json')
    assert replay['carrier_identities_checked'] == 900
    for rel, digest in replay['bindings'].items():
        assert g.sha(a.ROOT/rel) == digest
    g.save(OUT/'branch-analysis-plan.json', {'sources': {str(p.relative_to(a.ROOT)): g.sha(p) for p in
        [Path(__file__), OUT/'geometry.json', OUT/'verified.json']}, 'seconds': 60,
        'selection': 'unique minimum-genus cover in each already frozen target row; no new curves',
        'tests': 'exact trace identities and pairwise branch gcds within each fibre'})
    results = []
    for cid in ['302', '11952']:
        rows = [r for r in geometry['rows'] if r['fibre'] == cid]
        best = []
        traces = []
        A, B, sections, t0 = a.parent(cid)
        for row in rows:
            covers = [c for c in row['covers'] if c['normalization_genus'] == row['best_genus_in_atlas']]
            assert len(covers) == 1 and covers[0]['normalization_genus'] == 1
            c = covers[0]
            q = R(c['branch_monic_polynomial'])
            assert q.degree() == 4 and q.is_monic() and q.gcd(q.derivative()) == 1
            assert not c['branch_at_infinity']
            x0, x1, y0, y1 = [g.dec(c['maps'][k]) for k in ['x0', 'x1', 'y0', 'y1']]
            assert x1
            m = y1/x1
            tx = m*m-2*x0
            ty = m*(x0-tx)-y0
            assert ty*ty == tx**3+A*tx+B
            traces.append((tx, ty))
            best.append({'target': row['target'], 'cover_index': row['covers'].index(c),
                         'pencil': c['tag'], 'branch_monic_polynomial': c['branch_monic_polynomial'],
                         'constant_twist_representative': c['constant_twist_representative'],
                         'trace_section': [g.enc(tx), g.enc(ty)]})
        pairs = []
        for i, x in enumerate(best):
            for y in best[i+1:]:
                d = R(x['branch_monic_polynomial']).gcd(R(y['branch_monic_polynomial']))
                pairs.append({'targets': [x['target'], y['target']], 'gcd': list(map(str, d.list()))})
        disjoint = all(p['gcd'] == ['1'] for p in pairs)
        n = len(best)
        record = {'fibre': cid, 'selected_carriers': best, 'pairwise_branch_gcds': pairs,
                  'same_exact_generic_trace': all(T == traces[0] for T in traces),
                  'distinct_pencil_parameters': len({c['pencil']['parameter'] for c in best}),
                  'pairwise_disjoint_geometric_branch_divisors': disjoint,
                  'quadratic_squareclass_rank_over_algebraic_closure_of_constants': n if disjoint else None,
                  'compositum_degree': str(2**n) if disjoint else None,
                  'compositum_genus': str(1+2**n*(n-1)) if disjoint else None,
                  'scope': 'These selected minimum-genus carriers only. Other ancestry or shared-cover constructions are not excluded.'}
        results.append(record)
        print('BRANCH', cid, 'pairs', len(pairs), 'disjoint', disjoint,
              'same trace', record['same_exact_generic_trace'], 'compositum degree', record['compositum_degree'], flush=True)
    g.save(OUT/'branch-analysis.json', {'status': 'PASS_EXACT_SELECTED_CARRIER_BRANCH_ANALYSIS', 'rows': results,
          'proof': 'Distinct squarefree quartics with pairwise gcd1 have disjoint branch supports. Valuation at a branch of each quartic proves the radicals independent even over Qbar(t). The 2^n-cover has4n branch points with inertia2; Riemann-Hurwitz gives genus1+2^n(n-1).',
          'boundaries': 'This is a minimum common-cover degree for the emitted quadratic extensions, not for all multisection ancestries of the targets; a pencil is not a single quadratic field.',
          'bindings': {str(p.relative_to(a.ROOT)): g.sha(p) for p in [Path(__file__), OUT/'branch-analysis-plan.json', OUT/'geometry.json']}})


if __name__ == '__main__':
    main()
