#!/usr/bin/env python3
"""Render the elliptic-curve inventory and ICARM-like curve data from certificates."""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import certify_compact_r17_candidates as cert
from local_conductor_database import load_conductor_inventory

ROOT = Path(__file__).resolve().parents[2]
REPO = next(p for p in ROOT.parents if (p/'.git').exists()) if not (ROOT/'.git').exists() else ROOT
PREFIX = str(ROOT.relative_to(REPO))+'/' if ROOT != REPO else ''
METRICS = ROOT/'artifacts/generated-results/elliptic-curves/inventory201_table_metrics_v1.json'
OUT = ROOT/'elliptic-curves/data/research_curves'
BEGIN = '<!-- BEGIN GENERATED ELLIPTIC CURVE TABLE -->'
END = '<!-- END GENERATED ELLIPTIC CURVE TABLE -->'


def encoded(value):
    return json.dumps(value,indent=2,sort_keys=True)+'\n'


def put(path,text,check):
    if check:
        if path.read_text() != text:
            raise ArithmeticError('generated curve view differs: '+str(path))
    else:
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(text)


def log_integer(n):
    return math.log(int(n))


def run(check=False):
    metrics = json.loads(METRICS.read_text())
    byid = {r['id']:r for r in metrics['rows']}
    rows = []
    for r in load_conductor_inventory():
        m = byid[r['id']]
        info = r['conductor_information']
        if not cert.isomorphic(r['curve'],m['ainvs']):
            raise ArithmeticError('table model identity differs')
        source = json.loads((ROOT/info['certificate']).read_text()) if info.get('certificate') else {}
        known_primes = [q['prime'] for q in source.get('local_data',[]) if q['conductor_exponent']]
        rank_source = r['source_certificate']
        if not (ROOT/rank_source).is_file():
            rank_source = 'artifacts/generated-results/elliptic-curves/'+rank_source
        if not (ROOT/rank_source).is_file():
            raise FileNotFoundError('rank source certificate: '+rank_source)
        conductor = info['conductor'] if info['status']=='EXACT' else None
        if conductor is not None and not info['bad_primes']:
            raise ArithmeticError('exact high-rank row has no complete prime list')
        row = {'id':r['id'],'ainvs':m['ainvs'],'rank_lower_bound':r['rank_lower_bound'],
            'local_search_rank_lower_bound':r['local_search_rank_lower_bound'],
            'rank_provenance':r['rank_provenance'],'family':r['family'],'parameter':r['parameter'],
            'icarm_ids':r['current_catalogue_matches'],'conductor':conductor,
            'conductor_status':info['status'],'log_conductor':log_integer(conductor) if conductor else None,
            'conductor_divisor':info.get('conductor_divisor',conductor),
            'conductor_upper_bound':info.get('conductor_upper_bound',conductor),
            'bad_primes':info.get('bad_primes'),'known_bad_primes':info.get('bad_primes') or known_primes,
            'discriminant':m['discriminant'],'log_abs_discriminant':m['log_abs_discriminant'],
            'naive_height':m['naive_height'],'faltings_height':m['faltings_height'],
            'points':m['points'],'original_to_minimal_isomorphism':m['original_to_minimal_isomorphism'],
            'conductor_certificate':info.get('certificate'),'rank_source_certificate':rank_source,
            'rank_certificate':r['rank_certificate'],
            'rank_certificate_model':'Frozen original inventory model; points in this export are transported to the displayed minimal model by the saved isomorphism.'}
        rows.append(row)
    rows.sort(key=lambda r:(-r['rank_lower_bound'],r['conductor'] is None,
                            int(r['conductor'] or r['conductor_upper_bound'] or 10**1000),r['id']))
    counts = Counter(r['conductor_status'] for r in rows)
    data = {'schema':'elliptic-curves.research-table.v1','count':len(rows),'curves':rows,
        'conductor_status_counts':dict(counts),'logarithm':'natural','artifact_path_base':str(ROOT.relative_to(REPO)),
        'height_convention':metrics['definitions'],
        'provenance':{'metrics':str(METRICS.relative_to(ROOT)),
            'metrics_sha256':hashlib.sha256(METRICS.read_bytes()).hexdigest(),
            'conductor_manifest':'elliptic-curves/data/conductor_screen_current.json',
            'conductor_manifest_sha256':hashlib.sha256((ROOT/'elliptic-curves/data/conductor_screen_current.json').read_bytes()).hexdigest()},
        'claim_boundary':'Rank lower bounds, not exact ranks. Logarithms/heights are rounded approximations. Exact conductor and complete bad_primes are null until certified; partial information is kept in separate fields. Includes already-public matches and one public-point rank28 reproduction.'}
    put(OUT/'database.json',encoded(data),check)
    buf = io.StringIO(newline='')
    fields = ['id','ainvs','rank_lower_bound','conductor','conductor_status','log_conductor',
              'naive_height','faltings_height','discriminant','log_abs_discriminant','bad_primes','conductor_upper_bound']
    writer = csv.DictWriter(buf,fieldnames=fields,lineterminator='\n');writer.writeheader()
    for r in rows:
        writer.writerow({k:json.dumps(r[k]) if isinstance(r[k],list) else r[k] for k in fields})
    put(OUT/'database.csv',buf.getvalue(),check)
    for r in rows:
        put(OUT/(r['id']+'.json'),encoded(r),check)
        public = ', '.join(f'[ICARM#{i}](https://elliptic-rank.icarm.cloud/curve/{i})' for i in r['icarm_ids']) or 'Unmatched in the current ICARM snapshot'
        text = [f'# {r["id"]}','',f'Rank **≥ {r["rank_lower_bound"]}**. Family `{r["family"]}` at `{r["parameter"]}`. {public}.','',
            f'[Full data and transported points]({r["id"]}.json) · [Inventory](../../INVENTORY.md)',
            '', 'Minimal a-invariants:', '', '```text',', '.join(r['ainvs']),'```','',
            f'Conductor status: **{r["conductor_status"]}**.','']
        if r['conductor']:
            text += ['Exact conductor:','','```text',r['conductor'],'```','','Complete bad primes:','','```text',', '.join(r['bad_primes']),'```','']
        else:
            text += ['The exact conductor and complete bad-prime list remain unknown.','']
            if r['conductor_divisor']:
                text += ['Certified conductor divisor:','','```text',r['conductor_divisor'],'```','','Certified conductor upper bound:','','```text',r['conductor_upper_bound'],'```','']
            if r['known_bad_primes']:
                text += ['Proved bad primes so far (incomplete; not a submission-ready list):','','```text',', '.join(r['known_bad_primes']),'```','']
        text += ['Exact minimal discriminant:','','```text',r['discriminant'],'```','',
            '| log N | Naive height | Faltings height | log abs(Δ) |','|---:|---:|---:|---:|',
            f'| {r["log_conductor"]:.4f} | {r["naive_height"]:.4f} | {r["faltings_height"]:.4f} | {r["log_abs_discriminant"]:.4f} |' if r['conductor'] else f'| — | {r["naive_height"]:.4f} | {r["faltings_height"]:.4f} | {r["log_abs_discriminant"]:.4f} |','',
            'Natural logarithms; numerical columns are rounded. The Faltings column uses ICARM’s minimal-model period-area convention.','',
            f'Point count: {len(r["points"])}. Local search bound: ≥{r["local_search_rank_lower_bound"]}. Rank provenance: `{r["rank_provenance"]}`.','',
            f'[Rank source](../../../{r["rank_source_certificate"]}) · [Minimal model and transport checks](../../../{METRICS.relative_to(ROOT)})']
        if r['conductor_certificate']:
            text += [f'· [Conductor certificate](../../../{r["conductor_certificate"]})']
        put(OUT/(r['id']+'.md'),'\n'.join(text)+'\n',check)
    table = [BEGIN,'## Elliptic curve inventory','',
        f'**201 research curves · {counts["EXACT"]} exact conductors · {counts["UNKNOWN"]} unresolved**'+(f' · {counts["REPORTED"]} reported only' if counts['REPORTED'] else '')+'.',
        'Includes seven ICARM matches; the rank-28 row is a public-point reproduction. Rank values are proved lower bounds.','',
        'Columns and height conventions follow [ICARM’s table](https://elliptic-rank.icarm.cloud/curves). Logs are natural and shown to two decimals. A dash means the exact conductor is unknown; bounds and partial primes are available on the linked curve page. Coefficients are clipped here; each page contains the complete equation and data.',
        '', '[Download JSON](elliptic-curves/data/research_curves/database.json) · [Download CSV](elliptic-curves/data/research_curves/database.csv) · [Arithmetic and replay notes](elliptic-curves/notes/INVENTORY201_TABLE_AND_CONDUCTORS_2026-09-07.md)',
        '', '| Curve | a-invariants | Rank | log N | Naive height | Faltings height | log abs(Δ) |','|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        name = f'[{r["id"]}](elliptic-curves/data/research_curves/{r["id"]}.md)'
        if r['icarm_ids']:
            name += ' '+', '.join(f'[#{i}](https://elliptic-rank.icarm.cloud/curve/{i})' for i in r['icarm_ids'])
        ainvs = '['+', '.join(a if len(a)<=14 else a[:14]+'…' for a in r['ainvs'])+']'
        ln = f'{r["log_conductor"]:.2f}' if r['conductor'] else '—'
        table.append(f'| {name} | `{ainvs}` | ≥ {r["rank_lower_bound"]} | {ln} | {r["naive_height"]:.2f} | {r["faltings_height"]:.2f} | {r["log_abs_discriminant"]:.2f} |')
    table += ['',END]
    section = '\n'.join(table)
    inventory = ROOT/'elliptic-curves/INVENTORY.md'
    inventory_section = section.replace('](elliptic-curves/','](')
    inventory_old = inventory.read_text() if inventory.exists() else '# Elliptic-curve inventory\n'
    if BEGIN in inventory_old:
        a = inventory_old.index(BEGIN);b = inventory_old.index(END,a)+len(END)
        inventory_new = inventory_old[:a]+inventory_section+inventory_old[b:]
    else:
        inventory_new = inventory_old.rstrip()+'\n\n'+inventory_section+'\n'
    put(inventory,inventory_new,check)
    print('INVENTORY TABLE PASS',len(rows),'curves;',dict(counts))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
