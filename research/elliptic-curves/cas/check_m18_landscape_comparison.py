"""Hash/label replay; --replay-cvp repeats all 960 exact metric CVPs.

Uses only the self-contained frozen input archive for feature recomputation.
Outcome summaries are consulted after validating the feature seal.
"""
import argparse
import json
from pathlib import Path
import tarfile
import tempfile
from m18_landscape_comparison import OUT, ROOT, features, read, sha, require

def check(replay=False):
    manifest=read(OUT/'manifest.json')
    for name,digest in manifest['files'].items():require(sha(OUT/name)==digest,'package hash: '+name)
    protocol=read(OUT/'protocol.json');seal=read(OUT/'feature-seal.json')
    require(sha(OUT/'protocol.json')==seal['protocol_sha256'],'protocol seal')
    for name,digest in protocol['source_hashes'].items():require(sha(ROOT/name)==digest,'feature source hash')
    for name,digest in seal['features'].items():require(sha(OUT/name)==digest,'feature seal')
    with tempfile.TemporaryDirectory() as temp:
        temp=Path(temp)
        with tarfile.open(OUT/'frozen-inputs.tar.gz') as tar:
            for member in tar.getmembers():
                require(member.isfile() and not Path(member.name).is_absolute() and '..' not in Path(member.name).parts,'unsafe archive entry')
            tar.extractall(temp,filter='data')
        for case in protocol['cases']:
            folder=temp/case['case']
            for name,record in case['inputs'].items():require(sha(folder/name)==record['sha256'],'frozen input mismatch')
            saved=read(OUT/'features'/f"{case['case']}.json")
            if replay:
                result=features(read(folder/'seed-input.json'),read(folder/'selection.json'),folder)
                for key,value in result.items():require(value==saved[key],'feature replay: '+key)
                print('REPLAY',case['case'],flush=True)
        # Label verification is deliberately after checking every feature hash.
        outcomes=read(OUT/'outcomes.json')
        require(outcomes['feature_seal_sha256']==sha(OUT/'feature-seal.json'),'label seal binding')
        for row in outcomes['rows']:
            frozen=read(OUT/'features'/f"{row['case']}.json")
            require(row['features']==frozen['features'] and row['j']==frozen['j'],'joined feature mismatch')
            label=row['outcome']
            if label['rank'] is not None:
                source=temp/'outcome-sources'/f"{row['case']}.json"
                require(sha(source)==label['sha256'],'label source hash')
                raw=read(source)
                require(raw['status'] in ['PASS_INDEPENDENT_SEEDED_V3_REPLAY','PASS_INDEPENDENT_DET1092_REPLAY'],'unverified label')
                require(raw['initial_rank']==18 and raw['rank_lower_bound']==label['rank'] and raw['charts']==label['charts'] and raw['stop_reason']==label['stop_reason'],'outcome changed')
    print('PASS',len(protocol['cases']),'states;', '960 CVPs recomputed' if replay else 'hashes, feature seals, archived inputs and outcome bindings')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--replay-cvp',action='store_true');args=p.parse_args();check(args.replay_cvp)
