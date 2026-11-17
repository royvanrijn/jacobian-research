#!/usr/bin/env python3
"""Exact two-source comparison for the24 A6 affine certificates.

Default: retained snapshot consistency, with no external workbook read.
--workbook PATH: hash and re-extract the original XLSX read-only via ZIP XML.
--output PATH: retain the actual mode and pins of this invocation.
No canonical receipt is overwritten unless --output explicitly names it.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import platform
import xml.etree.ElementTree as ET
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
AUDIT = HERE/'a6_priority_audit.json'
NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
BOOK_SHA = '514b12e118fa681520da5702a5ed4307ce7cf93b6b9f7b5e53f265ca7b01647b'
BOOK_BLOB = '3882c6da34c9057f8f13da31afb9d8973b7158cb'
OPEN_SHA = '3fae47725d8808d58a6b1a24e36ade4e19121092c92714020ae0e1d0070b1edc'
OPEN_BLOB = '01c7c332a0bf1565d7558aec4c769ff36cbd31a9'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def blob(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()

def snapshot_check(audit,inputs,open_bytes):
    assert audit['positive_inputs'] == 'a6_affine_inputs.json'
    names = [case['name'] for case in inputs['cases']]
    assert len(names) == len(set(names)) == 24
    assert audit['workbook_sha256'] == BOOK_SHA
    assert audit['workbook_git_blob'] == BOOK_BLOB
    assert audit['workbook_bytes'] == 23554211
    source = audit['still_open_source']
    assert source['path'] == 'sources/bhattacharjee-still-open.json'
    assert source['sha256'] == sha(open_bytes) == OPEN_SHA
    assert source['git_blob'] == blob(open_bytes) == OPEN_BLOB
    open_names = json.loads(open_bytes)
    assert len(open_names) == len(set(open_names)) == 293
    assert set(names) <= set(open_names), 'name missing from retained293 list'
    assert [row['name'] for row in audit['rows']] == names, 'rows differ from all24 inputs'
    assert len({row['row'] for row in audit['rows']}) == 24
    assert set(audit['column_headers']) == set('ABCDEFGHIJ')
    for row in audit['rows']:
        number = str(row['row'])
        cells = row['cached_cells']
        assert set(cells) == {c+number for c in 'ABCDEFGHIJ'}
        assert cells['A'+number] == row['name'].replace('_','')
        assert cells['C'+number] == '4', 'recorded Wirtinger number differs'
        assert all(cells[c+number] == '0' for c in 'DEGI'), 'recorded Coxeter hit'
        assert all(cells[c+number] == 'N/A' for c in 'FHJ')
        assert row['formula_cells'] == []
    return names

def workbook_check(path,audit):
    data = path.read_bytes()
    assert len(data) == audit['workbook_bytes']
    assert sha(data) == BOOK_SHA and blob(data) == BOOK_BLOB
    wanted = {row['row']:row for row in audit['rows']}
    found = set()
    with ZipFile(path) as archive:
        shared = ET.fromstring(archive.read('xl/sharedStrings.xml'))
        strings = [''.join(item.itertext()) for item in shared]
        for _,row in ET.iterparse(archive.open('xl/worksheets/sheet1.xml'),events=('end',)):
            if row.tag.rsplit('}',1)[-1] != 'row':
                continue
            number = int(row.attrib['r'])
            if number == 1 or number in wanted:
                cells,formulas = {},[]
                for cell in row:
                    value = cell.find('s:v',NS)
                    if value is not None:
                        cells[cell.attrib['r']] = (strings[int(value.text)]
                            if cell.attrib.get('t') == 's' else value.text)
                    if cell.find('s:f',NS) is not None:
                        formulas.append(cell.attrib['r'])
                if number == 1:
                    assert cells == {c+'1':v for c,v in audit['column_headers'].items()}
                else:
                    assert cells == wanted[number]['cached_cells']
                    assert formulas == wanted[number]['formula_cells']
                    found.add(number)
            row.clear()
            if found == set(wanted):
                break
    assert found == set(wanted), 'original workbook rows incomplete'

def rejection(call):
    try:
        call()
    except (AssertionError,KeyError,IndexError,ValueError):
        return True
    raise AssertionError('corrupt comparison accepted')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    audit = json.loads(AUDIT.read_text())
    input_path = HERE/audit['positive_inputs']
    inputs = json.loads(input_path.read_text())
    open_path = HERE/audit['still_open_source']['path']
    open_bytes = open_path.read_bytes()
    names = snapshot_check(audit,inputs,open_bytes)
    bad = deepcopy(audit)
    row = bad['rows'][0]
    row['cached_cells']['D'+str(row['row'])] = '1'
    controls = {'changed_recorded_hit':rejection(lambda:snapshot_check(bad,inputs,open_bytes))}
    bad = deepcopy(audit)
    bad['rows'][0]['name'] = '13a_1'
    controls['changed_name'] = rejection(lambda:snapshot_check(bad,inputs,open_bytes))
    if args.workbook:
        workbook_check(args.workbook,audit)
    paths = [Path(__file__).resolve(),AUDIT,input_path,open_path]
    report = {'status':'PASS','python_version':platform.python_version(),
              'mode':('original pinned workbook re-extraction plus snapshot checks'
                      if args.workbook else 'retained snapshot consistency only'),
              'workbook_reextracted':bool(args.workbook),
              'names_checked':names,'positive_count':24,'still_open_source_count':293,
              'zero_recorded_maximal_finite_coxeter_hits':24,
              'negative_controls':controls,'scope':audit['scope'],
              'input_sha256':{p.relative_to(HERE).as_posix():sha(p.read_bytes()) for p in paths}}
    if args.workbook:
        report['workbook_sha256'] = BOOK_SHA
        report['workbook_git_blob'] = BOOK_BLOB
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
