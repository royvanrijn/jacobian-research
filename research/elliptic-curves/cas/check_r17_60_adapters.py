#!/usr/bin/env python3
"""Zero-search regression of new adapters on one retained, replayed seed."""
import os
from pathlib import Path
import r17_60_arithmetic as arithmetic
from v3_warm_support import atomic, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
fixture = ROOT/'artifacts/local/elliptic-curves/r17-60-adapter-check-v1'
source = ROOT/'artifacts/generated-results/elliptic-curves/fresh6_lowheight_first_m18_v1/result.json'
row = next(r for r in read(source)['records'] if r['family']=='074d9')
case = 'retained-adapter-control'
origin = (ROOT/row['evidence']['terminal.json']['path']).parent
folder = fixture/'cases'/case
folder.mkdir(parents=True,exist_ok=True)
link = folder/'seed-search'
if not link.exists():
    link.symlink_to(origin, target_is_directory=True)
require(link.resolve()==origin.resolve(),'retained control path differs')
atomic(fixture/'roster.json',{'rows':[{'id':case,'family':row['family'],'parameter':row['parameter']}]},immutable=True)
atomic(fixture/'protocol.json',{'status':'RETAINED_INPUT_ADAPTER_CONTROL','point_searches':0,
                              'source_sha256':sha(source)},immutable=True)
arithmetic.D = fixture
arithmetic.reconcile_seed(case)
arithmetic.reconcile_seed(case)
arithmetic.complement(case)
arithmetic.complement(case)
packet = read(folder/'seed-reconciled.json')
bank = read(folder/'complement-preparation/anchor-bank.json')
require(packet['rank_lower_bound']>=18 and len(bank['rows'])==16,'adapter output differs')
print('R17_60_ADAPTER_CHECK_PASS|retained_seed=074d9|rank='+str(packet['rank_lower_bound'])+'|parents=16|new_charts=0',flush=True)
