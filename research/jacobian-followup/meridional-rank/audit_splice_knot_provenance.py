#!/usr/bin/env python3
"""Exact independent HT-census/archived-table identity for the four cases."""
import csv
import hashlib
import importlib.metadata as metadata
import io
import json
from pathlib import Path
import zipfile
import snappy
import snappy_15_knots
from spherogram import Link
from audit_diagram_identity import normalize,pd_walk

HERE=Path(__file__).resolve().parent


def main():
    targets=list(json.loads((HERE/'three_sheet_splice_plan.json').read_text())['model_assignments'])
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    snapshots={version:{r['name']:r for r in json.loads((HERE/f'sources/knotinfo/{version}_selected.json').read_text())['records']} for version in ['2026.10.5','2026.11.2']}
    archive=next((HERE/'tmp/sources').glob('*2026.11.2*.whl'))
    sha=hashlib.sha256(archive.read_bytes()).hexdigest();assert sha=='5a91e8b77d24fb7d90a4a7048074acb97ca81cde205873983344a64458ad88a3'
    with zipfile.ZipFile(archive) as z:
        member=next(p for p in z.namelist() if p.endswith('/knotinfo_data_complete.csv'))
        reader=csv.DictReader(io.TextIOWrapper(z.open(member),encoding='utf-8-sig',newline=''),delimiter='|')
        fields=['name','classical_conway_name','conway_notation','montesinos_notation','two_bridge_notation','bridge_index','determinant','three_genus']
        family_data={r['name']:{k:r[k] for k in fields} for r in reader if r['name'] in targets}
    pinned=json.loads((HERE/'sources/current13/data/still_open.json').read_text())
    current=json.loads((HERE/'sources/current13/still_open_2026-11-17.json').read_text())
    deposited={r['name'] for r in json.loads((HERE/'sources/current13/data/certificates.json').read_text())}
    versions={name:metadata.version(name) for name in ['snappy','spherogram','snappy_15_knots']}
    assert versions=={'snappy':'3.3.2','spherogram':'2.4.1','snappy_15_knots':'1.2.1'}
    results=[]
    for name in targets:
        source=panel[name]['pd'];link=Link('K'+name.replace('_',''));pd=[list(c) for c in link.PD_code()]
        assert len(link.link_components)==1 and len(pd)==13
        assert all(snapshot[name]['pd_notation']==source and snapshot[name]['bridge_index']=='[3,4]' for snapshot in snapshots.values())
        match,change=normalize(pd_walk(pd),True);source_match,source_change=normalize(pd_walk(source),True)
        assert match==source_match
        literal=[[e+1 for e in c] for c in pd]==source
        assert name in pinned and name in current and name not in deposited
        results.append({'name':name,'census_lookup':'K'+name.replace('_',''),'census_pd':pd,'census_DT':link.DT_code(),
                        'literal_PD_matches_after_adding_one':literal,'exact_full_walk_match':True,
                        'source_walk_change':source_change,'census_walk_change':change,
                        'archived_KnotInfo_releases':['2026.10.5','2026.11.2'],'archived_bridge_interval':[3,4],
                        'in_pinned_and_fetched_project_open_list':True,'absent_from243_deposited_certificates':True,
                        'archived_family_metadata':family_data[name]})
    paths=['audit_splice_knot_provenance.py','three_sheet_splice_plan.json','joint_a5_panel.json','audit_diagram_identity.py',
           'sources/knotinfo/2026.10.5_selected.json','sources/knotinfo/2026.11.2_selected.json',
           'sources/current13/data/still_open.json','sources/current13/data/certificates.json','sources/current13/still_open_2026-11-17.json']
    report={'status':'PASS exact named PD provenance','versions':versions,'results':results,
            'KnotInfo_wheel_sha256':sha,'KnotInfo_csv_member':member,
            'project_status_fetch':{'date':'2026-11-17','url':'https://github.com/anaxia-research/conjecture-to-choose/blob/main/data/still_open.json','connector_reported_modified_date':'2026-10-20T14:48:34Z','names':len(current)},
            'priority_boundary':'Four names remain in the fetched project open list and archived[3,4] rows. This is not an exhaustive prior-family, live KnotInfo or worldwide novelty/priority audit; empty family metadata is not a proof of non-membership.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'splice_knot_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'targets':targets,'literal_PD_matches':[r['literal_PD_matches_after_adding_one'] for r in results],'family_metadata':family_data},indent=2))


if __name__=='__main__':main()
