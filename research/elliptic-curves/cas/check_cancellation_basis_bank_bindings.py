#!/usr/bin/env sage -python
"""Successor commissioning after explicit object/file binding was added."""
import copy
from pathlib import Path
import resource
import shutil
import time
import traceback

from cancellation_basis_early import guard, RAW as RETAINED
from cancellation_basis_epoch import read, sha, need
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL
from verify_cancellation_basis_bank import verify_selection

OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1'
RAW = LOCAL/'cancellation-basis-bank-certificate-v1'
CAS = Path(__file__).resolve().parent


def run():
    resource.setrlimit(resource.RLIMIT_CPU, (20,25)); start = time.process_time()
    plan, cases = guard(); case = next(c for c in cases if c['id'] == 'fcf6ce8e6a518eb715d2')
    sources = source_closure([Path(__file__)])
    recipe = {'status':'FROZEN_BOUND_SELECTION_SUCCESSOR','source_sha256':sources,
        'prior_development_sha256':sha(OUT/'development.json'),
        'prior_corruption_checks_sha256':sha(OUT/'corruption-checks.json'),
        'case':case['id'],'epochs':[1,2],'cpu_soft_seconds':20,'cpu_hard_seconds':25,'point_search_calls':0,
        'change':'Bind the supplied selection to the exact bytes read before verification; preserve their hash in the returned certificate. Earlier source snapshots and commissioning receipts remain unchanged.'}
    new_write(RAW/'binding-checks/recipe.json',recipe)
    for name in sources:
        target = RAW/'binding-checks/sources'/name; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/name,target)
    rows = []
    try:
        for epoch in recipe['epochs']:
            old = RETAINED/'arms'/case['id']/f'basis_refresh/epoch-{epoch:02d}'
            packet = read(old/'seed.json'); selection = read(old/'producer/selection.json')
            proof = verify_selection(packet,case,plan,old/'producer',selection)
            need(proof['selection_sha256'] == sha(old/'producer/selection.json'),'selection hash differs')
            new_write(RAW/'binding-checks'/f'epoch-{epoch}-proof.json',proof)
            rows.append({'epoch':epoch,'rank':proof['rank'],'status':proof['status'],'cpu_seconds':proof['cpu_seconds']})
            changed = copy.deepcopy(selection); changed['centres'].reverse()
            try:
                verify_selection(packet,case,plan,old/'producer',changed)
            except ArithmeticError as error:
                need(str(error) == 'selection object differs from bound file','unexpected binding rejection')
                rows.append({'epoch':epoch,'status':'REJECTED_UNBOUND_SELECTION','error':str(error)})
            else:
                raise ArithmeticError('unbound selection accepted')
        result = {'status':'PASS_BOUND_SELECTION_SUCCESSOR','rows':rows,'cpu_seconds':time.process_time()-start,
            'checker_sha256':sha(CAS/'verify_cancellation_basis_bank.py'),
            'recipe_sha256':sha(RAW/'binding-checks/recipe.json'),'point_search_calls':0}
    except BaseException as error:
        new_write(OUT/'binding-checks.json',{'status':'RETAINED_BINDING_CHECK_FAILURE','rows':rows,
            'error':type(error).__name__+': '+str(error),'traceback':traceback.format_exc(),
            'recipe_sha256':sha(RAW/'binding-checks/recipe.json'),'cpu_seconds':time.process_time()-start})
        raise
    new_write(OUT/'binding-checks.json',result)
    print({k:v for k,v in result.items() if k != 'rows'},flush=True)


if __name__ == '__main__':
    run()
