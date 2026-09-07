#!/usr/bin/env python3
"""Keep all branch polynomials and the failed interface run in one archive."""
import argparse
from pathlib import Path
import zipfile
import retrospective as r
import additive_branch_geometry as first
import complete_additive_branch_geometry as second
import additive_node_branch as node

BUNDLE=r.OUT/'rank_jump_additive_branch_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_additive_branch_evidence_v1.json'
PATHS=(first.OUTPUT,second.OUTPUT,node.OUTPUT)


def package():
    entries=[]
    with zipfile.ZipFile(BUNDLE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in PATHS:
            raw=p.read_bytes();info=zipfile.ZipInfo(p.name,(2026,9,7,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,raw,compresslevel=9);entries.append({'path':str(p.relative_to(r.ROOT)),'member':p.name,'sha256':r.digest(raw),'bytes':len(raw)})
    n=r.read(node.OUTPUT);assert n['status']=='PASS';rows=[]
    for row in r.read(second.OUTPUT)['rows']:
        if row['family']=='a1-fibration-01':
            rows.append({'family':row['family'],'status':'PASS_AFTER_NODE_FACTOR',
                'norm_count':len(n['rows']),'branch_degree_distribution':n['branch_degree_distribution'],
                'proved_coprime_pairs':n['pair_count'],'simultaneous_norm_cover_genus_ranges_first_four':n['simultaneous_norm_cover_genus_ranges_first_four']})
        else:assert row['status']=='PASS';rows.append({'family':row['family'],'status':'PASS',**row['summary']})
    r.write_new(MANIFEST,{'schema':'rank-jump.additive-branch-evidence.v1','status':'PASS',
        'bundle':str(BUNDLE.relative_to(r.ROOT)),'bundle_sha256':r.digest(BUNDLE.read_bytes()),'members':entries,
        'rows':rows,'node_parameter':n['node_parameter'],'node_triples':n['node_triple_count'],
        'node_norm_multiplicity':n['node_norm_multiplicity'],'packager_sha256':r.digest(Path(__file__).read_bytes()),
        'boundary':'Geometric norm-carrier calculation. A1 results refer to the residual branch polynomials after removing the common fourth power. No rank or Selmer dimension is inferred.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['package','unpack']);args=p.parse_args()
    if args.mode=='package':package()
    else:
        m=r.read(MANIFEST);assert r.digest(BUNDLE.read_bytes())==m['bundle_sha256']
        with zipfile.ZipFile(BUNDLE) as z:
            for row,path in zip(m['members'],PATHS):
                raw=z.read(row['member']);assert r.digest(raw)==row['sha256']
                if path.exists():assert path.read_bytes()==raw
                else:
                    with path.open('xb') as f:f.write(raw)
    print('PASS',args.mode)
