#!/usr/bin/env python3
"""Publish the canonical inventory table as an expandable main README section."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from refresh_icarm_local_database import load_catalogue
from curve_table_highlights import highlight_rank_minima

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


def run(check=False):
    argv = [sys.executable,str(ROOT/'elliptic-curves/cas/render_inventory201_readme.py')]
    if check:
        argv.append('--check')
    subprocess.run(argv,check=True)
    inventory = ROOT/'elliptic-curves/INVENTORY.md'
    source = inventory.read_text()
    section = source[source.index(BEGIN):source.index(END)+len(END)]
    database = json.loads((inventory.parent/'data/research_curves/database.json').read_text())
    rows = database['curves']
    section, kept = select_section(section, rows, conductor_benchmarks(load_catalogue()['curves']))
    # Recompute after curation: minima in the displayed subset can differ from
    # the complete inventory. The helper removes previous numeric bold first.
    section = highlight_rank_minima(section)
    pending = sum(r['id'] not in kept and r['conductor_status'] != 'EXACT' for r in rows)
    explanation = (f'**Main-table selection: {len(kept)} of {len(rows)} curves.** '
        'Keep every certified lower bound ≥22, documented structural exceptions, '
        'and smaller-rank curves whose exact conductor meets or beats the '
        'rank-specific minimum in the [pinned ICARM snapshot](data/icarm_current.json). '
        'This is an editorial filter, not a live record claim. '
        f'The [complete inventory](INVENTORY.md) retains every curve, point and certificate; '
        f'{pending} hidden curves still have unresolved conductors and are not classified as high-conductor.\n\n'
        'Below-22 structural examples retained: '+ '; '.join(
            f'[{identifier}](data/research_curves/{identifier}.md) ({reason})'
            for identifier, reason in STRUCTURAL_EXCEPTIONS.items() if identifier in kept)+'.\n\n')
    section = section.replace('| Curve |', explanation+'| Curve |', 1)
    def link(match):
        target = match.group(1)
        if '://' in target or target.startswith('#'):
            return match.group(0)
        path,separator,anchor = target.partition('#')
        relative = (inventory.parent/path).resolve().relative_to(REPO)
        return ']('+str(relative)+(separator+anchor if separator else '')+')'
    section = re.sub(r'\]\(([^)]+)\)',link,section)
    count = len(kept)
    section = section.replace('| Curve |',f'<details>\n<summary>Show {count} highlighted curves</summary>\n\n| Curve |',1)
    section = section.replace('\n'+END,'\n</details>\n\n'+END,1)
    path = REPO/'README.md'
    previous = path.read_text()
    if BEGIN in previous:
        a,b = previous.index(BEGIN),previous.index(END)+len(END)
        expected = previous[:a]+section+previous[b:]
    else:
        expected = previous.rstrip()+'\n\n'+section+'\n'
    if check:
        if previous != expected:
            raise ArithmeticError('main README table differs from canonical inventory')
    else:
        path.write_text(expected)
    print(f'MAIN README TABLE PASS: {count} highlighted / {len(rows)} archived rows')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
