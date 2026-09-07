#!/usr/bin/env sage-python
"""Finite height-six section domain for three-section elliptic pencils."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,ZZ,matrix,pari
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-all-section-triangles-v1'
def main():
    p=json.loads((D/'protocol.json').read_text())
    for name,h in p['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    h=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    H=matrix(ZZ,4*matrix(QQ,h['seed_height_gram']))
    U=matrix(ZZ,pari(H).qflllgram());assert abs(U.det())==1
    reduced=U.transpose()*H*U
    data=pari(reduced).qfminim(24,100000,2)
    short=matrix(ZZ,data[2]).transpose();assert int(data[0])==2*short.nrows()
    vectors=short*U.transpose()
    result={'schema':'elliptic-curves.mestre-height-six.v1','status':'ENUMERATED_PENDING_INDEPENDENT_COMPLETENESS',
            'scaled_gram':[list(map(int,r)) for r in H.rows()],
            'basis_change':[list(map(int,r)) for r in U.rows()],
            'reduced_gram':[list(map(int,r)) for r in reduced.rows()],
            'reduced_vectors_up_to_sign':[list(map(int,r)) for r in short.rows()],
            'section_words_up_to_sign':[list(map(int,r)) for r in vectors.rows()],
            'signed_vector_count':int(data[0]),'scaled_norm_bound':24,'sources':p['sources']}
    with (D/'height-domain.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS enumerated signed height<=6 words',int(data[0]),flush=True)
if __name__=='__main__':main()
