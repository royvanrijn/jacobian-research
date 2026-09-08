#!/usr/bin/env python3
"""Replay the frozen control, including measured floating-point timing fields.

The original CLI wrongly passed those descriptive floats to the exact-only
protocol digest. Preserve its source and compare the complete JSON instead.
"""
import json
import audit_factor_free_control as audit

if __name__=='__main__':
    actual=audit.expected();saved=audit.cert.read(audit.OUT)
    if json.dumps(actual,sort_keys=True)!=json.dumps(saved,sort_keys=True):
        raise ArithmeticError('retained control differs')
    print(actual['status'],'rank',actual['rank_lower_bound'],'points',len(actual['points']))
