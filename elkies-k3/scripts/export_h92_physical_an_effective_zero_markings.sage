#!/usr/bin/env sage -python
"""Export effective-zero markings for a certified child of ADE type nA1."""
import argparse, hashlib, json
from pathlib import Path
from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]; U2=matrix(ZZ,((0,1),(1,0)))
ap=argparse.ArgumentParser(); ap.add_argument('--source-marking',type=Path,required=True); ap.add_argument('--certificate',type=Path,required=True); ap.add_argument('--root-rank',type=int,required=True); ap.add_argument('--output-prefix',type=Path,required=True); ap.add_argument('--summary-output',type=Path,required=True); ap.add_argument('--extra-source-curve',action='append',default=[],help='NAME:x0,...,x18 for an additional exact effective curve in source-frame coordinates'); args=ap.parse_args()
SOURCE=args.source_marking.resolve(); CERT=args.certificate.resolve(); PREFIX=args.output_prefix.resolve(); OUTPUT=args.summary_output.resolve(); rr=args.root_rank
def load(p): return matrix(ZZ,[[ZZ(x) for x in l.split()] for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')])
def ent(v): return [int(x) for x in vector(ZZ,v)]
def rows(m): return [[int(x) for x in r] for r in m.rows()]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
s=json.loads(SOURCE.read_text()); c=json.loads(CERT.read_text()); assert c['status']=='PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE'
frame=load(ROOT/s['frame_output']); gram=block_diagonal_matrix(U2,-frame); fibre=vector(ZZ,c['source_to_child_basis'][0]); explicit={n:vector(ZZ,v) for n,v in s['equation_explicit_curves_in_child'].items()}
for item in args.extra_source_curve:
 name,raw=item.split(':',1); value=vector(ZZ,[ZZ(x) for x in raw.split(',')]); assert name and name not in explicit and len(value)==19 and value*gram*value==-2; explicit[name]=value
zeros=sorted(n for n,v in explicit.items() if v*gram*fibre==1); assert zeros
outputs={}
for zn in zeros:
 z=explicit[zn]; u=matrix(ZZ,[list(fibre),list(z+fibre)]); ker=matrix(ZZ,[list(fibre*gram),list((z+fibre)*gram)]).right_kernel_matrix(); pre=u.stack(ker); assert abs(pre.det())==1; pinv=pre.inverse().change_ring(ZZ); comp=-(pre*gram*pre.transpose())[2:,2:]
 half=[vector(ZZ,col) for col in matrix(ZZ,pari(comp).qfminim(2)[2]).columns()]; roots=half+[-r for r in half]; assert len(roots)==2*rr
 known=[]
 for curve in {tuple(v):v for v in explicit.values()}.values():
  if curve*gram*fibre: continue
  inc=ZZ(z*gram*curve); assert inc in (0,1); pos=curve if inc==0 else fibre-curve; tail=vector(ZZ,(pos*pinv)[2:]);
  if tuple(tail) not in {tuple(x) for x in known}: known.append(tail)
 chosen=[]; seen=set(); oldf=vector(ZZ,[1,0]+[0]*17); pinned=vector(ZZ,s['target_fibres_in_root_adapted_hub']['pinned_R17'])
 for root in roots:
  if tuple(root) in seen: continue
  seen|={tuple(root),tuple(-root)}; pair=(root,-root); pick=[r for r in pair if tuple(r) in {tuple(k) for k in known}]
  if not pick:
   pick=[r for r in pair if oldf*gram*(vector(ZZ,[0,0]+list(r))*pre)>0]
  if not pick:
   pick=[r for r in pair if pinned*gram*(vector(ZZ,[0,0]+list(r))*pre)>0]
  if not pick: pick=[min(pair,key=lambda r:tuple(r))]
  assert len(pick)==1; chosen.append(pick[0])
 assert len(chosen)==rr; simple=[vector(ZZ,[0,0]+list(r))*pre for r in chosen]; partial=u.stack(matrix(ZZ,[list(-x) for x in simple])); smith,left,right=partial.smith_form(); n=2+rr; assert smith[:,:n]==matrix.identity(ZZ,n) and not any(smith[:,n:].list())
 basis=block_diagonal_matrix(left.inverse().change_ring(ZZ),matrix.identity(ZZ,19-n))*right.inverse().change_ring(ZZ)
 for i in range(n,19):
  row=basis.row(i); row-=ZZ(row*gram*(z+fibre))*fibre; row-=ZZ(row*gram*fibre)*(z+fibre); basis[i]=row
 assert abs(basis.det())==1; inverse=basis.inverse().change_ring(ZZ); childgram=basis*gram*basis.transpose(); child=-childgram[2:,2:]; assert child[:rr,:rr]==2*matrix.identity(ZZ,rr)
 slug=zn.lower(); fp=Path(str(PREFIX)+f'-{slug}-frame.txt'); mp=Path(str(PREFIX)+f'-{slug}-marking.json'); fp.parent.mkdir(parents=True,exist_ok=True); fp.write_text(f'# physical {rr}A1; zero={zn}\n'+'\n'.join(' '.join(map(str,r)) for r in child.rows())+'\n'); eq=matrix(ZZ,s['equation_A11_to_root_adapted_hub_basis'])
 payload={'schema':'elkies-k3.h3-physical-an-effective-zero-marking.v1','status':'PASS_EXACT_PHYSICAL_AN_EFFECTIVE_ZERO_MARKING','hub':f'physical_{rr}A1_{zn}_zero','zero':zn,'root_data':[rr,2*rr,2**rr],'ade':f'{rr}A1','frame_output':str(fp.relative_to(ROOT)),'frame_sha256':sha(fp),'basis_in_source':rows(basis),'source_in_basis':rows(inverse),'equation_A11_to_root_adapted_hub_basis':rows(basis*eq),'target_fibres_in_root_adapted_hub':{n:ent(vector(ZZ,v)*inverse) for n,v in s['target_fibres_in_root_adapted_hub'].items()},'equation_explicit_curves_in_child':{n:ent(v*inverse) for n,v in explicit.items()},'prefix_operational_score':None,'proof_boundary':'Exact effective zero, nA1 component signs fixed successively by inherited vertical roots, the physical old fibre, and pinned-R17 nef anchor; all transports unimodular.','inputs':{'paths':[str(SOURCE.relative_to(ROOT)),str(CERT.relative_to(ROOT))],'sha256':{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,CERT)}}}; mp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n'); outputs[zn]={'frame':str(fp.relative_to(ROOT)),'marking':str(mp.relative_to(ROOT))}
summary={'schema':'elkies-k3.h3-physical-an-effective-zero-markings.v1','status':'PASS_EXACT_PHYSICAL_AN_ALL_EFFECTIVE_ZERO_MARKINGS','root_rank':rr,'outputs':outputs}; OUTPUT.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n'); print(f'PHYSANZEROS|rank={rr}|count={len(outputs)}|output={OUTPUT}')
