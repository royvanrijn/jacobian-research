#!/usr/bin/env sage -python
"""Compose exact marked-neighbour certificates and pin a rootless endpoint."""
import argparse, hashlib, json
from pathlib import Path
from sage.all import ZZ, identity_matrix, matrix, pari

ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser(); ap.add_argument('--edge',type=Path,action='append',required=True); ap.add_argument('--source',required=True); ap.add_argument('--status',required=True); ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
PINNED=ROOT/'elkies-k3/data/lattice/rank17_gram.txt'
def load(p): return matrix(ZZ,[[ZZ(x) for x in l.split()] for l in p.read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')])
def rows(m): return [[int(x) for x in r] for r in m.rows()]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):
 p=p.resolve()
 try:return str(p.relative_to(ROOT))
 except ValueError:return str(p)
def qiso(A,B):
 raw=pari(A).qfisom(pari(B)); assert str(raw)!='0'; C=matrix(ZZ,raw); choices=[C,C.transpose()]
 Ci=C.inverse()
 if Ci.change_ring(ZZ)==Ci: Ci=Ci.change_ring(ZZ); choices += [Ci,Ci.transpose()]
 for Q in choices:
  if Q*A*Q.transpose()==B: assert abs(Q.det())==1; return Q
 raise ArithmeticError('no verified row-convention qfisom')
paths=tuple(p.resolve() for p in args.edge); ps=[json.loads(p.read_text()) for p in paths]; route=[]
for p,d in zip(paths,ps):
 assert d['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'; F=matrix(ZZ,d['source_to_child_basis']); I=matrix(ZZ,d['child_to_source_basis']); E=matrix(ZZ,d['equation_A11_to_child_basis']); Ei=matrix(ZZ,d['child_to_equation_A11_basis']); ident=identity_matrix(ZZ,19)
 assert F*I==I*F==ident and E*Ei==Ei*E==ident and abs(F.det())==abs(E.det())==1; assert d['first_edge_exact_horizontal_nef_gate'] and not d['first_edge_exact_negative_horizontal_walls']
 route.append({'certificate':rel(p),'sha256':sha(p),'source_hub':d['source_hub'],'candidate_id':d['candidate_id'],'child_root_data':d['child']['root_data'],'child_mw_rank':d['child']['mw_rank'],'forward_determinant':int(F.det()),'inverse_determinant':int(I.det())})
last=ps[-1]; assert last['child']['root_data']==[0,0,1] and last['child']['mw_rank']==17; childp=ROOT/last['frame_output']; child=load(childp); pinned=load(PINNED); assert child.det()==pinned.det()==948; Q=qiso(child,pinned)
out=args.output.resolve(); out.parent.mkdir(parents=True,exist_ok=True); result={'schema':'elkies-k3.h3-marked-route-to-pinned-r17.v1','status':args.status,'source':args.source,'route':route,'q_sequence':[x['candidate_id']['q'] for x in route],'old_fibre_degrees':[x['candidate_id']['old_fibre_degree'] for x in route],'final_child_frame':rel(childp),'final_child_frame_sha256':sha(childp),'pinned_rank17_frame':rel(PINNED),'pinned_rank17_frame_sha256':sha(PINNED),'rootless_child_to_pinned_r17_isometry':rows(Q),'isometry_relation':'Q * child * Q^t = pinned rank17_gram.txt','isometry_determinant':int(Q.det()),'equation_A11_to_rootless_child_basis':last['equation_A11_to_child_basis'],'rootless_child_to_equation_A11_basis':last['child_to_equation_A11_basis'],'proof_boundary':'Every edge supplies an exact primitive nef isotropic class, marked U, complete finite horizontal-wall audit, full root data, and mutually inverse unimodular NS transports. The terminal rootless positive frame is identified exactly with pinned R17. Equation lifts are separate.','inputs':{'paths':[rel(p) for p in paths+(PINNED,)],'sha256':{rel(p):sha(p) for p in paths+(PINNED,)}}}; out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n'); print(f'MARKEDR17|edges={len(route)}|q={",".join(map(str,result["q_sequence"]))}|det={Q.det()}|status={args.status}|output={out}')
