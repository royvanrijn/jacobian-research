"""Bind live producer tuples through the bank's declared JSON representation.

The earlier checker and failed timing harness remain source-frozen. This
adapter changes serialization only; all arithmetic checks are delegated to
the same independent checker on exactly the saved JSON selection.
"""
from fractions import Fraction as F
import json

from cancellation_basis_epoch import sha, need
from finite_cancellation_corpus import canonical, digest, write
from verify_cancellation_basis_bank import verify_selection as verify_json_selection


def verify_selection(packet, case, plan, folder, selection):
    return verify_json_selection(packet,case,plan,folder,json.loads(canonical(selection)))


def rebuild(packet, case, plan, folder):
    from cancellation_scheduler_cpu import cpu
    from visibility_complement_subset import landscape
    from v3_warm_engine import certified_state
    start = cpu(); model = tuple(map(F,packet['curve'])); basis = [tuple(map(F,p)) for p in packet['points']]
    certified_state(model,basis,packet['proof']); write(folder/'seed.json',packet)
    parent = case['parent_bank']; policy = {**plan['bank_policy'],'generic_rank':parent['dimension'],
        'scaled_shells':parent['shells'],'anchor_count':len(parent['rows'])}
    need(len(basis)-parent['dimension'] <= plan['maximum_extension_dimension'],'extension limit exceeded')
    tick = cpu(); selection = landscape(model,basis,set(),folder/'producer',policy,parent); producer_cpu = cpu()-tick
    tick = cpu(); proof = verify_selection(packet,case,plan,folder/'producer',selection); verifier_cpu = cpu()-tick
    write(folder/'selection-proof.json',proof)
    bank = {'seed':packet,'centres':selection['centres'][:plan['maximum_centres']],'rank':len(basis)}
    need(bank['centres'],'empty bank'); write(folder/'bank.json',bank)
    receipt = {'status':'PASS_INDEPENDENT_BANK_CERTIFICATE','rank':len(basis),
        'bank_sha256':sha(folder/'bank.json'),'seed_sha256':sha(folder/'seed.json'),
        'producer_sha256':sha(folder/'producer/selection.json'),'selection_proof_sha256':sha(folder/'selection-proof.json'),
        'parent_bank_sha256':digest(canonical(parent)),'exported_centres':len(bank['centres']),
        'production_cpu_seconds':producer_cpu,'independent_verification_cpu_seconds':verifier_cpu,
        'cpu_seconds':cpu()-start,'boundary':proof['boundary'],'serialization':'Exact JSON round trip before bound-selection verification; native producer tuples are JSON arrays.'}
    write(folder/'verification.json',receipt)
    return bank,receipt
