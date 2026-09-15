#!/usr/bin/env sage -python
"""Independent retained-output accounting, not another ellipsoid/search run."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import time

from cancellation_basis_bank_cost_v2 import guard, OUT, RAW, ARMS
from cancellation_basis_epoch import read, sha, need
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import canonical, digest, ROOT
from run_complement_seed_v3 import normalized_selection


def verify(write_result=True):
    import numpy as np
    started = time.process_time(); plan, inputs = guard(); ph = sha(OUT/'protocol.json')
    summary = read(OUT/'summary.json'); supervision = read(OUT/'supervision.json')
    need(summary['protocol_sha256'] == ph == supervision['protocol_sha256'] and supervision['status'] == 'COMPLETE', 'comparison incomplete')
    need(len(supervision['records']) == len(summary['rows']) == 68, 'worker count differs')
    expected = {(t['id'],a):t for t in inputs for a in ARMS}; seen = set(); totals = [0.,0.]
    by_family = defaultdict(dict); strata = defaultdict(lambda:[0.,0.]); counters = Counter(); certificates = []
    minimum_vectors = 0; direct_norms = 0
    for supervised,reported in zip(supervision['records'],summary['rows']):
        key = (supervised['id'],supervised['arm']); need(key in expected and key not in seen,'duplicate or unexpected worker'); seen.add(key)
        item = expected[key]; arm = supervised['arm']; dest = RAW/'arms'/item['id']/arm
        result = read(dest/'result.json'); receipt = read(dest/'bank/verification.json')
        need(read(dest/'supervisor.json') == supervised and
             all(reported[k] == v for k,v in supervised.items()),'supervisor/report mismatch')
        need(supervised['status'] == 'COMPLETE' and supervised['returncode'] == 0 and
             supervised['protocol_sha256'] == result['protocol_sha256'] == ph and
             result['status'] == 'PASS_IDENTICAL_CERTIFIED_BANK' and result['point_search_calls'] == 0,'unsuccessful worker')
        need(supervised['result_sha256'] == sha(dest/'result.json') and supervised['worker_log_sha256'] == sha(dest/'worker.log'),'worker file binding differs')
        need(result['verification_sha256'] == sha(dest/'bank/verification.json') and
             result['native_rank_sha256'] == sha(dest/'native-rank.json'),'proof binding differs')
        rank = read(dest/'native-rank.json'); bank = read(dest/'bank/bank.json'); old = ROOT/item['retained']
        need(rank['status'] == 'PASS' and rank['certificate']['rank'] == item['rank'] and
             rank['seed_sha256'] == digest(canonical(item['packet'])),'native subgroup receipt differs')
        need(bank == read(old/'bank.json') and bank['seed'] == item['packet'] and bank['rank'] == item['rank'] and
             sha(dest/'bank/bank.json') == sha(old/'bank.json') == result['bank_sha256'] == receipt['bank_sha256'], 'exported bank bytes differ')
        selection = read(dest/'bank/producer/selection.json')
        need(normalized_selection(selection) == normalized_selection(read(old/'producer/selection.json')), 'complete selection differs')
        need(receipt['seed_sha256'] == sha(dest/'bank/seed.json') and
             receipt['producer_sha256'] == sha(dest/'bank/producer/selection.json') and
             receipt['parent_bank_sha256'] == digest(canonical(item['case']['parent_bank'])), 'bank input binding differs')
        for ai,anchor in enumerate(selection['anchors']):
            need(sha(dest/'bank/producer'/f'anchor-{ai:02d}-full.npz') == anchor['full_scores_sha256'] and
                 sha(dest/'bank/producer'/f'anchor-{ai:02d}-maps.json') == anchor['maps_sha256'], 'producer payload differs')
        if arm == ARMS[0]:
            need(receipt['status'] == 'PASS_INDEPENDENT_BANK' and
                 receipt['reference_sha256'] == sha(dest/'bank/reference/selection.json') and
                 normalized_selection(selection) == normalized_selection(read(dest/'bank/reference/selection.json')), 'baseline reference differs')
        else:
            proof_path = dest/'bank/selection-proof.json'; proof = read(proof_path)
            need(receipt['status'] == 'PASS_INDEPENDENT_BANK_CERTIFICATE' and receipt['selection_proof_sha256'] == sha(proof_path) and
                 proof['status'] == 'PASS_INDEPENDENT_BANK_SELECTION_CERTIFICATE' and
                 proof['selection_sha256'] == sha(dest/'bank/producer/selection.json') and
                 proof['seed_sha256'] == digest(canonical(item['packet'])) and
                 proof['parent_bank_sha256'] == digest(canonical(item['case']['parent_bank'])), 'independent certificate differs')
            claimed = {(a['anchor'],a['extension']):a for a in proof['cvp_certificates']}
            need(len(claimed) == len(proof['cvp_certificates']),'duplicate minimum packet')
            n = item['rank']; g = selection['rounded_gram']; u = selection['LLL']; generic = selection['generic_rank']
            matched = set()
            for ai,anchor in enumerate(selection['anchors']):
                for row in anchor['refined']:
                    idx = (ai,row['extension']); need(idx in claimed, 'omitted minimum packet'); matched.add(idx)
                    certified = claimed[idx]; original = row['cvp']
                    need(certified['status'] == 'PASS_EXACT_PARITY_MINIMUM' and certified['norm'] == original['norm'] and
                         sorted(certified['minima']) == sorted(original['minima']) and
                         0 < certified['verification_nodes'] <= original['nodes'] <= plan['bank_plan']['bank_policy']['exact_cvp_node_limit'], 'minimum certificate differs')
                    parity = [v % 2 for v in anchor['anchor']['representative'][:generic]]
                    parity += [(row['extension'] >> i) & 1 for i in range(n-generic)]
                    for vector in certified['minima']:
                        word = [sum(vector[i]*u[i][j] for i in range(n)) for j in range(n)]
                        need([w % 2 for w in word] == parity, 'certificate parity transport differs')
                        norm = sum(word[i]*g[i][j]*word[j] for i in range(n) for j in range(n))
                        need(norm == certified['norm'], 'direct integer norm differs')
                        minimum_vectors += 1; direct_norms += 1
            need(matched == set(claimed),'extra minimum packet')
            counters['independent_minimum_packets'] += len(claimed)
            counters['verification_nodes'] += sum(p['verification_nodes'] for p in claimed.values())
            counters['optimizer_nodes'] += sum(r['cvp']['nodes'] for a in selection['anchors'] for r in a['refined'])
            counters['native_word_checks'] += proof['native_word_checks']; counters['map_recipes'] += proof['map_recipes']
            certificates.append({'id':item['id'],'selection_proof_sha256':sha(proof_path),'cvp_packets':len(claimed)})
        which = ARMS.index(arm); elapsed = supervised['charged_cpu_seconds']; need(0 < elapsed < 35,'worker CPU bound differs')
        totals[which] += elapsed; strata[item['rank']][which] += elapsed
        values = by_family[item['case']['family']].setdefault(item['case']['id'],[0.,0.]); values[which] += elapsed
    need(seen == set(expected),'incomplete roster')
    need(sum(totals) <= plan['campaign_child_cpu_seconds'] and
         abs(sum(totals)-supervision['charged_child_cpu_seconds']) < 1e-9,'campaign CPU differs')
    for i,arm in enumerate(ARMS): need(abs(totals[i]-summary['full_child_cpu_seconds'][arm]) < 1e-9,'reported CPU differs')
    rng = np.random.default_rng(20260914); boot = []
    families = [list(rows.values()) for _,rows in sorted(by_family.items())]
    for _ in range(10000):
        a = b = 0.
        for rows in families:
            for ix in rng.integers(len(rows),size=len(rows)):
                a += rows[int(ix)][0]; b += rows[int(ix)][1]
        boot.append(a/b)
    interval = [float(np.quantile(boot,q)) for q in (.0125,.9875)]
    need(np.allclose(interval,summary['paired_family_whole_curve_bootstrap_central_97_5_percent'],rtol=0,atol=1e-12),'bootstrap differs')
    saving = 1-totals[1]/totals[0]
    gate = saving >= .25 and interval[0] > 1 and all(a/b > 1 for a,b in strata.values())
    need(bool(gate) == summary['engineering_gate'] and abs(saving-summary['candidate_cpu_saving_fraction']) < 1e-12,'engineering gate differs')
    result = {'status':'PASS_RETAINED_BANK_COST_ACCOUNTING','protocol_sha256':ph,'summary_sha256':sha(OUT/'summary.json'),
        'checker_sha256':sha(Path(__file__)),'workers':len(seen),'identical_bank_pairs':34,
        'fresh_native_subgroup_receipts':68,'counts':dict(counters),'minimum_vectors':minimum_vectors,
        'direct_integer_norm_checks':direct_norms,'certificate_bindings':certificates,
        'full_child_cpu_seconds':dict(zip(ARMS,totals)),'engineering_gate':bool(gate),
        'cpu_seconds':time.process_time()-started,'point_search_calls':0,
        'boundary':'Independent byte, seed, native-rank receipt, minimum-packet, exact integer norm/parity, paired timing and bootstrap accounting. The primary workers performed full independent ellipsoid and native group checks. This accounting does not rerun ellipsoids, native rank proofs, maps or point searches; timing claims are finite execution observations.'}
    if write_result: new_write(OUT/'audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k != 'certificate_bindings'},indent=2),flush=True)
    return result


if __name__ == '__main__':
    verify()
