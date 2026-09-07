#!/usr/bin/env sage -python
"""Attach the unique inherited effective zero to the physical A2/MW15 hub."""
import hashlib,json
from pathlib import Path
from sage.all import ZZ,block_diagonal_matrix,matrix,pari,vector
ROOT=Path(__file__).resolve().parents[2]; GEN=ROOT/'artifacts/generated-results'; U=matrix(ZZ,((0,1),(1,0)))
SOURCE=GEN/'elkies-k3-h3-physical-q12-4a1-old_a11_component_0-marking.json'; CERT=GEN/'elkies-k3-h3-physical-4a1-q4o21633-a2-certificate.json'; FRAME=GEN/'elkies-k3-h3-physical-a2-effective-zero-frame.txt'; OUTPUT=GEN/'elkies-k3-h3-physical-a2-effective-zero-marking.json'
def load(p):return matrix(ZZ,[[ZZ(x) for x in l.split()] for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')])
def ent(v):return[int(x) for x in vector(ZZ,v)]
def rows(m):return[[int(x) for x in r]for r in m.rows()]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
s=json.loads(SOURCE.read_text());c=json.loads(CERT.read_text());sf=load(ROOT/s['frame_output']);g=block_diagonal_matrix(U,-sf);fibre=vector(ZZ,c['source_to_child_basis'][0]);explicit={n:vector(ZZ,v)for n,v in s['equation_explicit_curves_in_child'].items()};zeros=[n for n,v in explicit.items()if v*g*fibre==1];assert zeros==['second_old_I6_I4_missing_component'];z=explicit[zeros[0]]
u=matrix(ZZ,[list(fibre),list(z+fibre)]);ker=matrix(ZZ,[list(fibre*g),list((z+fibre)*g)]).right_kernel_matrix();pre=u.stack(ker);assert abs(pre.det())==1;pinv=pre.inverse().change_ring(ZZ);comp=-(pre*g*pre.transpose())[2:,2:];half=[vector(ZZ,col)for col in matrix(ZZ,pari(comp).qfminim(2)[2]).columns()];roots=half+[-r for r in half];assert len(roots)==6
known=explicit['old_A11_component_0'];assert z*g*known==0;known_tail=vector(ZZ,(known*pinv)[2:]);oldf=vector(ZZ,[1,0]+[0]*17);choices=[]
for other in roots:
 if known_tail*comp*other!=-1:continue
 other_source=vector(ZZ,[0,0]+list(other))*pre;affine=fibre-known-other_source
 if oldf*g*other_source>=0 and oldf*g*affine>=0:choices.append(other)
assert len(choices)==1;simple=[known_tail,choices[0]];simple_source=[vector(ZZ,[0,0]+list(r))*pre for r in simple]
partial=u.stack(matrix(ZZ,[list(-x)for x in simple_source]));smith,left,right=partial.smith_form();assert smith[:,:4]==matrix.identity(ZZ,4)and not any(smith[:,4:].list());basis=block_diagonal_matrix(left.inverse().change_ring(ZZ),matrix.identity(ZZ,15))*right.inverse().change_ring(ZZ)
for i in range(4,19):
 row=basis.row(i);row-=ZZ(row*g*(z+fibre))*fibre;row-=ZZ(row*g*fibre)*(z+fibre);basis[i]=row
assert abs(basis.det())==1;inverse=basis.inverse().change_ring(ZZ);childgram=basis*g*basis.transpose();child=-childgram[2:,2:];assert child[:2,:2].det()==3 and ZZ(pari(child[:2,:2]).qfminim(2)[0])==6;FRAME.write_text('# physical A2; zero=second_old_I6_I4_missing_component\n'+'\n'.join(' '.join(map(str,r))for r in child.rows())+'\n');eq=matrix(ZZ,s['equation_A11_to_root_adapted_hub_basis'])
payload={'schema':'elkies-k3.h3-physical-a2-effective-zero-marking.v1','status':'PASS_EXACT_PHYSICAL_A2_EFFECTIVE_ZERO_MARKING','hub':'physical_A2_second_old_I6_I4_missing_component_zero','zero':zeros[0],'root_data':[2,6,3],'ade':'A2','frame_output':str(FRAME.relative_to(ROOT)),'frame_sha256':sha(FRAME),'basis_in_source':rows(basis),'source_in_basis':rows(inverse),'equation_A11_to_root_adapted_hub_basis':rows(basis*eq),'target_fibres_in_root_adapted_hub':{n:ent(vector(ZZ,v)*inverse)for n,v in s['target_fibres_in_root_adapted_hub'].items()},'equation_explicit_curves_in_child':{n:ent(v*inverse)for n,v in explicit.items()},'prefix_operational_score':None,'proof_boundary':'Unique inherited effective zero; the second A2 component sign is the unique root adjacent to the known component for which both it and the resulting affine component have nonnegative old-fibre degree; all transports are unimodular.','inputs':{'paths':[str(SOURCE.relative_to(ROOT)),str(CERT.relative_to(ROOT))],'sha256':{str(p.relative_to(ROOT)):sha(p)for p in(SOURCE,CERT)}}};OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(OUTPUT)
