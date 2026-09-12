#!/usr/bin/env python3
"""Render a compact main README summary from the existing curve inventory.

The legacy selection helpers remain as regressions for the former highlighted
table. The full table and all curve pages live in elliptic-curves/INVENTORY.md.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
REPO = next(p for p in ROOT.parents if (p/'.git').exists()) if not (ROOT/'.git').exists() else ROOT
BEGIN = '<!-- BEGIN GENERATED ELLIPTIC CURVE TABLE -->'
END = '<!-- END GENERATED ELLIPTIC CURVE TABLE -->'

# Editorial exceptions, not claims of high rank or small conductor. One exported
# representative per construction, plus the smaller equation-derived control.
STRUCTURAL_EXCEPTIONS = {
    'det1092-orbit8044-000000': 'parametrized seed factory',
    'det1092-progression-n0': 'uniform arithmetic seed progression',
    'det1092-bifibration-47755-u0': 'alternate-fibration conic seed',
    'det1092-small-conic': 'small equation-derived conic control',
}


def conductor_benchmarks(catalogue):
    """Pinned comparison only; public reported ranks are not local proofs."""
    return {rank: min((int(r['conductor']) for r in catalogue
                      if r['rank_lower_bound'] >= rank and r.get('conductor')),
                     default=None) for rank in range(0, 22)}


def retention_reason(row, benchmarks):
    rank = row['local_search_rank_lower_bound']
    if rank >= 22:
        return 'certified lower bound at least 22'
    if row['id'] in STRUCTURAL_EXCEPTIONS:
        return STRUCTURAL_EXCEPTIONS[row['id']]
    benchmark = benchmarks.get(rank)
    if (row.get('conductor_status') == 'EXACT' and row.get('conductor')
            and benchmark is not None and int(row['conductor']) <= benchmark):
        return 'exact conductor at or below pinned rank-specific benchmark'
    return None


def select_section(section, rows, benchmarks):
    """Filter the presentation only; never change the canonical census."""
    kept = {r['id'] for r in rows if retention_reason(r, benchmarks)}
    removed = {r['id'] for r in rows} - kept
    seen = set()
    lines = []
    for line in section.splitlines():
        match = re.match(r'\| \[([^]]+)\]\(data/research_curves/([^/)]+)\.md\)', line)
        if match:
            identifier = match.group(2)
            if identifier in removed:
                seen.add(identifier)
                continue
        lines.append(line)
    if seen != removed:
        raise ArithmeticError('unmatched hidden curve row: '+str(removed-seen))
    return '\n'.join(lines), kept


def short_introduction(database):
    counts = database['conductor_status_counts']
    return (BEGIN+'\n## Elliptic curve inventory\n\n'
        f'**{database["count"]} research curves · {counts.get("EXACT", 0)} exact conductors · {counts.get("UNKNOWN", 0)} unresolved**'
        +(f' · {counts["REPORTED"]} reported only' if counts.get('REPORTED') else '')+'.\n\n'
        '[Full inventory](INVENTORY.md) · [JSON](data/research_curves/database.json) · '
        '[CSV](data/research_curves/database.csv) · '
        '[Methods and selection](notes/INVENTORY_REFRESH_2026-09-09.md)\n\n'
        'Ranks are certified lower bounds. Full equations, points, conductor bounds and '
        'provenance are retained in the linked inventory.\n\n')


def run(check=False, from_inventory=False):
    argv = [sys.executable,str(ROOT/'elliptic-curves/cas/render_inventory201_readme.py')]
    if check:
        argv.append('--check')
    if not from_inventory:
        subprocess.run(argv,check=True)
    inventory = ROOT/'elliptic-curves/INVENTORY.md'
    database = json.loads((inventory.parent/'data/research_curves/database.json').read_text())
    rows = database['curves']
    if database['count'] != len(rows):
        raise ArithmeticError('inventory count differs from exported rows')
    if dict(Counter(r['conductor_status'] for r in rows)) != database['conductor_status_counts']:
        raise ArithmeticError('conductor status totals differ from exported rows')
    section = short_introduction(database)+END
    def link(match):
        target = match.group(1)
        if '://' in target or target.startswith('#'):
            return match.group(0)
        path,separator,anchor = target.partition('#')
        relative = (inventory.parent/path).resolve().relative_to(REPO)
        return ']('+str(relative)+(separator+anchor if separator else '')+')'
    section = re.sub(r'\]\(([^)]+)\)',link,section)
    path = REPO/'README.md'
    previous = path.read_text()
    if BEGIN in previous:
        a,b = previous.index(BEGIN),previous.index(END)+len(END)
        expected = previous[:a]+section+previous[b:]
    else:
        expected = previous.rstrip()+'\n\n'+section+'\n'
    if check:
        if previous != expected:
            raise ArithmeticError('main README summary differs from canonical inventory')
    else:
        path.write_text(expected)
    print(f'MAIN README SUMMARY PASS: {len(rows)} curves; full inventory retained')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--from-inventory',action='store_true',
                        help='Presentation-only render from the already generated inventory; do not rebuild arithmetic data.')
    args=parser.parse_args()
    run(args.check,args.from_inventory)
