#!/usr/bin/env python3
"""Audit retained banks and freeze a roster before any validation point calls."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import subprocess

from finite_cancellation_corpus import ROOT, OUT as V1, LOCAL, canonical, digest, write

OUT = ROOT / 'artifacts/generated-results/elliptic-curves/finite_cancellation_validation_v2'
FAMILIES = ('074d9', '07ca9', '08234', '08f72', '103b2', '11952')


def j_group(model):
    a1, a2, a3, a4, a6 = map(F, model)
    b2 = a1*a1 + 4*a2
    b4 = a1*a3 + 2*a4
    b6 = a3*a3 + 4*a6
    c4 = b2*b2 - 24*b4
    c6 = -b2**3 + 36*b2*b4 - 216*b6
    return digest(str(1728*c4**3/(c4**3-c6**2)).encode())


def main():
    if (OUT/'roster.json').exists():
        raise FileExistsError('Preserve the frozen roster; no outcome-based replacement.')
    # This design is written before examining retained bank eligibility.
    design = {
        'status': 'DESIGN_BEFORE_BANK_AUDIT', 'maximum_cases': 48,
        'strata': {'shallow': [18, 20], 'deep': [24, 27]},
        'per_family_per_stratum': 4, 'families': list(FAMILIES),
        'selection': 'All currently retained local */landscape/selection.json banks joined by exact j to the V1 corpus. Require a matching certified seed, at least 16 centres and a separately documented larger subgroup. Exclude all 12 V1 CPU curves and Curve302. For each j use the lexically first shallow bank; for deep use the highest starting rank then lexical bank. Select deep before shallow, lowest j hashes per family, at most four per family and stratum, with globally distinct j. Missing slots remain missing.',
        'maximum_centres': 48, 'height': 125000,
        'search_cpu_seconds_per_arm': 40, 'hard_process_cpu_seconds': 60,
        'point_wall_seconds': 5, 'arm_wall_seconds': 90,
        'arms': ['factor_free', 'real_q', 'real_q_g'],
        'fit': 'Fixed ridge penalty 10. Entire validation j groups, the original 12 CPU groups and Curve302 are excluded. No CPU outcomes or target heights select banks.',
        'feature_change': 'Both fitted arms share real features and a genuinely q-only residue tree. The full arm additionally computes the original two gcd expectations and their unresolved/missing masses. No outcome-driven tuning.',
        'primary': 'Certified next-direction recoveries divided by all charged process-tree CPU, including failed/unsuccessful arms. Literal withheld representatives counted separately.',
        'gate': 'For each fitted arm vs factor-free, require at least as many certified successes, at least 10 percent better aggregate recoveries per CPU, and a stratified paired case-bootstrap 95 percent interval entirely above 1. A cancellation-specific gate compares real_q_g with real_q by the same rule. Zero-success denominators make the rate gate UNKNOWN. Report the original-style 10 percent aggregate-time gate separately.',
        'uncertainty': '10000 paired bootstrap samples, sampling cases within each observed stratum and family, seed 20260914. Fixed sample/caps; no retries to improve timing or extra cases after outcomes.',
        'boundary': 'New retrospective controls, not new curves. The whole corpus has already informed development. This is separation from fitting and from new CPU outcomes, not a pristine external population. Cached centre orders are common inputs; no optimality/CVP theorem is claimed.'}
    write(OUT/'design.json', design)
    corpus = json.loads(gzip.decompress((V1/'corpus.json.gz').read_bytes()))
    old = json.loads((V1/'cpu/inputs.json').read_text())
    excluded = {c['j_group'] for c in old} | {c['j_group'] for c in corpus if c['family']=='Curve302-development'}
    endpoints = defaultdict(list)
    for c in corpus:
        if c['family'] in FAMILIES and c['j_group'] not in excluded:
            endpoints[c['j_group']].append(c)
    paths = sorted(p for p in subprocess.check_output(
        ['rg', '--files', str(LOCAL), '-g', 'selection.json'], text=True).splitlines()
        if p.endswith('/landscape/selection.json'))
    counts = Counter(); candidates = []; failures = []
    for raw in paths:
        lp = Path(raw); ep = lp.parent.parent
        sp = ep/'seed.json'
        if not sp.exists():
            sp = ep.parent/'seed.json'
        if not sp.exists():
            counts['missing_seed'] += 1; continue
        try:
            seed = json.loads(sp.read_text()); rank = len(seed.get('points', []))
            stratum = next((s for s,(lo,hi) in design['strata'].items() if lo<=rank<=hi), None)
            if stratum is None:
                counts['outside_rank_strata'] += 1; continue
            group = j_group(seed['curve'])
            cs = [c for c in endpoints[group] if c['rank_lower_bound']>rank]
            if not cs:
                counts['no_eligible_larger_endpoint_or_excluded'] += 1; continue
            selection = json.loads(lp.read_text())
            if selection.get('basis') != seed['points']:
                counts['basis_mismatch'] += 1; continue
            proof = seed.get('proof', {})
            if proof.get('rank_lower_bound')!=rank or not proof.get('signatures') or not proof.get('no_rational_2_torsion_prime'):
                counts['missing_seed_proof'] += 1; continue
            centres = [r['representative'] for r in selection.get('centres', [])[:48]]
            if len(centres)<16 or any(len(w)!=rank or any(type(v)!=int for v in w) or not any(w) for w in centres):
                counts['invalid_or_short_centre_bank'] += 1; continue
            c = min(cs, key=lambda r:(-r['rank_lower_bound'],r['family'],r['id']))
            row = {'j_group':group, 'family':c['family'], 'stratum':stratum, 'initial_rank':rank,
                   'corpus_id':c['id'], 'endpoint_rank':c['rank_lower_bound'],
                   'source_seed':str(sp.relative_to(ROOT)), 'source_landscape':str(lp.relative_to(ROOT)),
                   'seed_sha256':digest(sp.read_bytes()), 'landscape_sha256':digest(lp.read_bytes()),
                   'centres':len(centres)}
            candidates.append(row); counts['eligible_banks'] += 1
        except (KeyError, ValueError, TypeError, ZeroDivisionError) as e:
            failures.append({'path':str(lp.relative_to(ROOT)), 'error':type(e).__name__+': '+str(e)})
    chosen = []; used = set()
    for stratum in ('deep','shallow'):
        banks = {}
        key = (lambda r:(-r['initial_rank'],r['source_landscape'])) if stratum=='deep' else (lambda r:r['source_landscape'])
        for row in sorted((r for r in candidates if r['stratum']==stratum), key=key):
            banks.setdefault(row['j_group'], row)
        for family in FAMILIES:
            eligible = sorted((r for g,r in banks.items() if g not in used and r['family']==family), key=lambda r:r['j_group'])
            for row in eligible[:4]:
                chosen.append(row); used.add(row['j_group'])
    assert len(used)==len(chosen) and not used & excluded
    chosen.sort(key=lambda r:(r['stratum'],r['family'],r['j_group']))
    inputs=[]; oracle=[]
    for row in chosen:
        seed=json.loads((ROOT/row['source_seed']).read_text())
        selection=json.loads((ROOT/row['source_landscape']).read_text())
        ident=digest(canonical([row['j_group'],row['source_seed'],row['source_landscape']]))[:20]
        row['id']=ident
        inputs.append({'id':ident,'j_group':row['j_group'],'family':row['family'],'stratum':row['stratum'],
                       'seed':{k:seed[k] for k in ('curve','points','proof')},
                       'centres':[c['representative'] for c in selection['centres'][:48]]})
        c=next(c for c in corpus if c['id']==row['corpus_id'])
        oracle.append({'id':ident,'corpus_id':c['id'],'curve':c['curve'],'points':c['generic_points']+c['targets'],
                       'endpoint_rank':c['rank_lower_bound'],'source':c['source'],'proof_sha256':c['proof_sha256']})
    write(OUT/'inputs.json',inputs)
    with gzip.GzipFile(str(OUT/'oracle.json.gz'),'wb',mtime=0) as f:f.write(canonical(oracle))
    with gzip.GzipFile(str(OUT/'eligible_banks.json.gz'),'wb',mtime=0) as f:f.write(canonical(candidates))
    result={'status':'FROZEN_ROSTER_BEFORE_NEW_FIT_OR_CPU', 'design_sha256':digest((OUT/'design.json').read_bytes()),
            'landscape_files_seen':len(paths),'counts':dict(counts),'read_failures':failures,
            'candidate_j_groups':len({r['j_group'] for r in candidates}),'eligible_by_stratum_rank':dict(Counter(f"{r['stratum']}:M{r['initial_rank']}" for r in candidates)),
            'chosen':chosen, 'excluded_prior_j_groups':sorted(excluded),
            'excluded_training_j_groups':sorted(used|excluded),
            'inputs_sha256':digest((OUT/'inputs.json').read_bytes()),'oracle_sha256':digest((OUT/'oracle.json.gz').read_bytes()),
            'eligible_banks_sha256':digest((OUT/'eligible_banks.json.gz').read_bytes()),
            'corpus_sha256':digest((V1/'corpus.json.gz').read_bytes()),'source_sha256':digest(Path(__file__).read_bytes())}
    write(OUT/'roster.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('chosen','excluded_training_j_groups','excluded_prior_j_groups')},indent=2))
    print('SELECTED',dict(Counter((r['stratum'],r['initial_rank']) for r in chosen)),flush=True)


if __name__=='__main__':main()
