#!/usr/bin/env sage -python
"""Export one certified marked edge in the common neighbour-record schema."""
import argparse, hashlib, json
from pathlib import Path
from sage.all import ZZ, matrix, vector
ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser(); ap.add_argument('--marking',type=Path,required=True); ap.add_argument('--certificate',type=Path,required=True); ap.add_argument('--orbit-index',type=int,required=True); ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
def load(p): return matrix(ZZ,[[ZZ(x) for x in l.split()] for l in p.read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')])
def rel(p):
 p=p.resolve()
 try:return str(p.relative_to(ROOT))
 except ValueError:return str(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads(args.marking.resolve().read_text()); c=json.loads(args.certificate.resolve().read_text()); frame=load(ROOT/m['frame_output']); rr=int(m['root_data'][0]); fibre=vector(ZZ,c['source_to_child_basis'][0]); witness=vector(ZZ,fibre[2:]); root=vector(ZZ,witness[:rr]); mw=vector(ZZ,witness[rr:]); cartan=frame[:rr,:rr]; coupling=frame[:rr,rr:]; labels=cartan*root+coupling*mw
record={'q':int(fibre[0]*fibre[1]),'factor_order':[int(fibre[0]),int(fibre[1])],'old_fiber_degree':int(fibre[1]),'orbit_index':args.orbit_index,'mw_projection':[int(x) for x in mw],'dominant_labels':[int(x) for x in labels],'witness':[int(x) for x in witness],'fiber':[int(x) for x in fibre],'child_root_data':c['child']['root_data'],'child_mw_rank':c['child']['mw_rank'],'child_ade':'certified'}
out=args.output.resolve(); out.parent.mkdir(parents=True,exist_ok=True); payload={'status':'PASS_ROOT_ADAPTED_WEYL_NEIGHBORS','frame':m['frame_output'],'input_root_data':m['root_data'],'neighbors':[record],'summaries':[{'q':record['q'],'factor_order':record['factor_order'],'primitive_neighbors':1,'dominant_orbits':1,'dominant_orbits_complete':False}],'inputs':{'paths':[rel(args.marking),rel(args.certificate)],'sha256':{rel(p):sha(p.resolve()) for p in (args.marking,args.certificate)}}}; out.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n'); print(out)
