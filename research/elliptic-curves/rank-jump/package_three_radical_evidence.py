#!/usr/bin/env python3
"""Deterministic compressed preservation of the two exact norm certificates."""
import argparse
import json
from pathlib import Path
import zipfile
import retrospective as r
import three_radical_incidence as first
import three_radical_private_ramification as second

BUNDLE=r.OUT/'rank_jump_three_radical_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_three_radical_evidence_v1.json'
PATHS=(first.OUTPUT,second.OUTPUT)


def package():
    members=[]
    with zipfile.ZipFile(BUNDLE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for path in PATHS:
            raw=path.read_bytes();info=zipfile.ZipInfo(path.name,(2026,9,7,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,raw,compresslevel=9)
            members.append({'path':str(path.relative_to(r.ROOT)),'member':path.name,
                            'sha256':r.digest(raw),'bytes':len(raw)})
    a,b=map(r.read,PATHS);assert a['status']==b['status']=='PASS'
    rows=[]
    for x,y in zip(a['rows'],b['rows']):
        assert x['token']==y['token']
        rows.append({'token':x['token'],'generic_dimension':x['generic_dimension'],**x['summary'],
            'ramification_rank_lower_bound':y['ramification_rank_lower_bound'],
            'candidate_span_unramified_dimension_upper_bound':y['candidate_span_unramified_dimension_upper_bound'],
            'new_Selmer_dimension_after_adding_generic_subgroup':y['new_Selmer_dimension_after_adding_generic_subgroup']})
    r.write_new(MANIFEST,{'schema':'rank-jump.three-radical-evidence.v1','status':'PASS',
        'bundle':str(BUNDLE.relative_to(r.ROOT)),'bundle_sha256':r.digest(BUNDLE.read_bytes()),
        'members':members,'rows':rows,'packager_sha256':r.digest(Path(__file__).read_bytes()),
        'boundary':'Unlabelled frozen generic-only panel. Norm and private-factor certificates are preserved byte for byte in the ZIP, not discarded.'})


def read_member(path):
    with zipfile.ZipFile(BUNDLE) as z:return z.read(path.name)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['package','unpack']);args=p.parse_args()
    if args.mode=='package':package()
    else:
        manifest=r.read(MANIFEST);assert r.digest(BUNDLE.read_bytes())==manifest['bundle_sha256']
        for entry,path in zip(manifest['members'],PATHS):
            raw=read_member(path);assert r.digest(raw)==entry['sha256']
            if path.exists():assert path.read_bytes()==raw
            else:
                with path.open('xb') as f:f.write(raw)
    print('PASS',args.mode)
