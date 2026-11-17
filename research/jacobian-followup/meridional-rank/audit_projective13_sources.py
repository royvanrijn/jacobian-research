#!/usr/bin/env python3
"""Retain pinned geometry-source cues without inferring family membership."""

import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

from audit_knotinfo_snapshots import RELEASES

HERE=Path(__file__).resolve().parent


def main():
    result=json.loads((HERE/'projective13_cover_results.json').read_text())
    names=sorted(result['base_homology']);assert len(names)==8
    release=next(r for r in RELEASES if r['version']=='2026.11.2')
    archive=HERE/'tmp/sources'/release['url'].rsplit('/',1)[1]
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==release['sha256']
    columns=['name','conway_polynomial','three_genus','conway_notation','montesinos_notation','symmetry_type']
    csv.field_size_limit(10**9)
    with zipfile.ZipFile(archive) as z:
        member=next(m for m in z.namelist() if m.endswith('/knotinfo_data_complete.csv'))
        raw=z.read(member);reader=csv.reader(io.StringIO(raw.decode('utf-8-sig')),delimiter='|')
        header=next(reader);indices={k:header.index(k) for k in columns}
        records=[{k:row[i] for k,i in indices.items()} for row in reader if row[indices['name']] in names]
    assert sorted(r['name'] for r in records)==names
    retained={'release':release,'csv_member':member,'csv_sha256':hashlib.sha256(raw).hexdigest(),
        'selection':'Exactly8 source knots with new SL3(F3) quotients; archived release, not a live classification',
        'columns_retained':columns,'records':records}
    target=HERE/'sources/knotinfo/2026.11.2_projective13_geometry_subset.json'
    if target.exists():assert json.loads(target.read_text())==retained
    else:target.write_text(json.dumps(retained,indent=2)+'\n')
    groups=[]
    for D in sorted({r['H1_torsion_factors'][0] for r in result['base_homology'].values()}):
        selected=[r for r in records if result['base_homology'][r['name']]['H1_torsion_factors']==[D]]
        assert len(selected)==4 and len({r['conway_polynomial'] for r in selected})==1
        groups.append({'double_cover_H1_order':D,'names':sorted(r['name'] for r in selected),
            'common_Conway_polynomial':selected[0]['conway_polynomial'],
            'three_genus_values':sorted({r['three_genus'] for r in selected})})
    catalogue=json.loads((HERE/'sources.json').read_text())
    entry=next(r for r in catalogue['sources'] if r['id']=='KNOTINFO-SNAPSHOTS')
    supplement={'snapshot':'sources/knotinfo/2026.11.2_projective13_geometry_subset.json',
        'columns':columns,'assurance':'Pinned archived polynomial/genus/notation/symmetry fields read. Two4knot blocks have equal Conway polynomials and the new homology types. Blank notation is not a family exclusion; neither mutation nor common branched covers nor bridge value follows from these invariants.'}
    if 'projective13_supplementary_read' in entry:assert entry['projective13_supplementary_read']==supplement
    else:
        entry['projective13_supplementary_read']=supplement
        (HERE/'sources.json').write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    assert len(catalogue['sources'])==len({r['id'] for r in catalogue['sources']})>=40
    files=['audit_projective13_sources.py','audit_knotinfo_snapshots.py','sources.json',
        'sources/knotinfo/2026.11.2_projective13_geometry_subset.json','projective13_cover_verification.json']
    audit={'status':'PASS pinned projective13 structural source cues','catalogue_distinct_ids':len(catalogue['sources']),
        'source_release':'2026.11.2','groups':groups,'common_cover_or_mutation_proved':False,
        'next_structural_question':'Do the4knot blocks admit explicit Conway mutations or a common branched-cover presentation? Prove this geometrically/algebraically before reducing the actual rank problem; matching invariants do not suffice.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'projective13_source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({'status':audit['status'],'group_orders':[r['double_cover_H1_order'] for r in groups],'common_cover_or_mutation_proved':False}))


if __name__=='__main__':main()
