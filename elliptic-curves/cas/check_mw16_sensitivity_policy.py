#!/usr/bin/env python3
"""Audit declared bounds/centre coverage, then run the exact evidence checker.

The original verifier and its executed source snapshots remain unchanged.
This additional check binds each completed chart to the frozen policy label.
"""
from itertools import product
import verify_mw16_sensitivity as exact


def check_policy(payload):
    budget=payload['declared_budget']
    if [r['parent_id'] for r in payload['results']]!=budget['parent_ids']:
        raise ArithmeticError('candidate order differs from the declared policy')
    expected_settings=set(product(budget['centres'],budget['specifications'],budget['heights']))
    for row in payload['results']:
        actual=[(s['centre'],s['specification'],s['height']) for s in row['settings']]
        expected=set() if row.get('skipped_adaptive_reason') else expected_settings
        if len(actual)!=len(set(actual)) or set(actual)!=expected:
            raise ArithmeticError('setting coverage differs from the declared policy')
        deepest=set(row['deepest_masks'])
        if len(deepest)!=len(row['deepest_masks']):
            raise ArithmeticError('repeated deepest class')
        for setting in row['settings']:
            seen=[]
            for chart in setting['charts']:
                job=chart['centre_construction']
                if (chart['specification'],chart['height_bound'],job['centre'])!=(setting['specification'],setting['height'],setting['centre']):
                    raise ArithmeticError('chart bound or coordinate differs from its setting label')
                if budget['mode']=='adaptive':
                    if chart['mask']!=job['old_mask']+(job['quotient_word']<<16):
                        raise ArithmeticError('adaptive class differs from the declared quotient word')
                    seen.append((job['old_mask'],job['quotient_word']))
                else:
                    seen.append(chart['mask'])
            if budget['mode']=='adaptive':
                bits=len(row['search_basis'])-16
                if not 1<=bits<=budget['adaptive_bits']:
                    raise ArithmeticError('adaptive slice exceeds its declared bit cap')
                expected_charts=set(product(deepest,range(1,1<<bits)))
            else:
                expected_charts=deepest
            if len(seen)!=len(set(seen)) or set(seen)!=expected_charts:
                raise ArithmeticError('deepest-class or quotient-word coverage changed')


def main():
    verify=exact.verify
    def with_policy(payload,documents):
        check_policy(payload)
        return verify(payload,documents)
    exact.verify=with_policy
    exact.main()


if __name__=='__main__': main()
