#!/usr/bin/env python3
"""Use only T!=0 saturation for the parallel-line boundary branch."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=HERE/'weighted-unmarked7/reduced/parallel/input.json';d=json.loads(source.read_text());out=source.parent.parent/'parallel-simple';out.mkdir(exist_ok=True);d['original_reduced_input_sha256']=sha(source);d['source_sha256']=sha(Path(__file__));d['mode']='parallel-simple';d['required_nonzero']='f4';d['saturation']='f4';V=d['variables'];e=[0]*len(V);eq=[[list(e),'1']];e[V.index('f4')]=e[V.index('branch_inv')]=1;eq.append([e,'-1']);d['equations'][-1]=eq;d['scope']='Same complete linear boundary reduction as parallel. Saturate only T!=0; J=1 itself excludes T=1/4 by the critical boundary. This tests a stronger augmented ideal with a lower-degree saturation generator. No inference from a timeout.';(out/'input.json').write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(out)
