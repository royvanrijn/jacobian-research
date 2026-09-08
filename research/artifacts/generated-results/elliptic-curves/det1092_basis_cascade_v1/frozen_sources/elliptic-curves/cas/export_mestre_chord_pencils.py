#!/usr/bin/env python3
"""Preserve all nine compiled presentations and the complete bounded selection."""
import json,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/mestre-chord-pencils-v1';OUT=ART/'mestre_chord_pencils_v1.json'
def read(p):return json.loads(p.read_text())
def compute():
 paths=[ART/'mestre_u11_visible_two_neighbor_v1.json']+[D/('pencil'+str(i)+'.json') for i in range(8)]
 records=read(D/'compilation-ledger.json');assert records['terminal'] and len(records['rows'])==8
 baseline=read(ROOT/'artifacts/local/elliptic-curves/mestre-u11-visible-two-neighbor-v1/supervisor.json')
 geometry=read(D/'geometry-supervisor.json');roster=read(D/'roster-supervisor.json')
 for r in records['rows']+[baseline,geometry,roster]:assert r['outcome']=='completed' and r['returncode']==0
 files=paths+[Path(__file__).resolve(),D/'protocol.json',D/'roster.json',D/'admission.json',ART/'mestre_chord_pencil_geometry_v1.json',ART/'mestre_rational_ns_gram_v2.json']
 return {'schema':'elliptic-curves.mestre-chord-pencils.v1','status':'PASS','presentations':[read(p) for p in paths],'ns':read(ART/'mestre_rational_ns_gram_v2.json')['rows'][0],
 'protocol':read(D/'protocol.json'),'roster':read(D/'roster.json'),'admission':read(D/'admission.json'),'geometry':read(ART/'mestre_chord_pencil_geometry_v1.json'),
 'supervision':{'baseline':baseline,'roster':roster,'compilations':records['rows'],'geometry':geometry},
 'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
 'scope':'Nine exact genus-one presentations on the already certified u11 parent, with rational points and Jacobians. The eight screened presentations have exact generic Q ranks8 to10 and are not admitted for new point searches. No claim of pairwise inequivalence or all-fibration exclusion.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();result=compute()
 if args.check:assert result==read(OUT)
 else:
  assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n')
 print('PASS9 PRESERVED PENCILS AND COMPLETE781-WORD ROSTER')
