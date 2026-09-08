#!/usr/bin/env python3
"""Publish the canonical inventory table as an expandable main README section."""
import argparse
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
REPO = next(p for p in ROOT.parents if (p/'.git').exists()) if not (ROOT/'.git').exists() else ROOT
BEGIN = '<!-- BEGIN GENERATED ELLIPTIC CURVE TABLE -->'
END = '<!-- END GENERATED ELLIPTIC CURVE TABLE -->'


def run(check=False):
    argv = [sys.executable,str(ROOT/'elliptic-curves/cas/render_inventory201_readme.py')]
    if check:
        argv.append('--check')
    subprocess.run(argv,check=True)
    inventory = ROOT/'elliptic-curves/INVENTORY.md'
    source = inventory.read_text()
    section = source[source.index(BEGIN):source.index(END)+len(END)]
    def link(match):
        target = match.group(1)
        if '://' in target or target.startswith('#'):
            return match.group(0)
        path,separator,anchor = target.partition('#')
        relative = (inventory.parent/path).resolve().relative_to(REPO)
        return ']('+str(relative)+(separator+anchor if separator else '')+')'
    section = re.sub(r'\]\(([^)]+)\)',link,section)
    section = section.replace('| Curve |','<details>\n<summary>Show all 201 curves</summary>\n\n| Curve |',1)
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
    print('MAIN README TABLE PASS: 201 rows, synchronized with inventory')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
