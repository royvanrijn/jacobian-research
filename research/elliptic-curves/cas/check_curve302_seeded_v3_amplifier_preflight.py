#!/usr/bin/env python3
"""Foreground, zero-chart preflight for the curve-302 seeded V3 amplifier.

Run under Sage Python. This reconstructs/reuses both immutable rank-18 seeds,
replays their finite certificates, loads the frozen V3 numerical context, and
verifies the exact ordered basis. It deliberately does not build landscapes or
execute point-search charts.
"""
from run_curve302_seeded_v3_amplifier import SEEDS, construct_seed, context
from v3_warm_support import assert_basis


def main():
    for seed_id in SEEDS:
        construct_seed(seed_id)
        ctx=context(seed_id)
        assert_basis(ctx.state,ctx.model,ctx.state.basis,18)
        print(f'SEEDED_V3_PREFLIGHT_PASS|seed={seed_id}|rank={ctx.state.rank}',flush=True)
    print('SEEDED_V3_PREFLIGHT_COMPLETE|cases=2|charts=0|status=PASS',flush=True)


if __name__=='__main__':
    main()
