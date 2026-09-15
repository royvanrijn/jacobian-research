"""Check the actual source-basis rootless U, full nef criterion and MW basis classes."""
from sage.all import *
from pathlib import Path
import json,argparse
P=Path(__file__).resolve().parents[2]/'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
src=json.loads((P/'explicit-source-certificate.json').read_text());F=matrix(ZZ,src['frame_gram']);U=matrix(ZZ,[[0,1],[1,0]]);NS=block_diagonal_matrix(U,-F);z=json.loads((P/'source-rootless-nef20.json').read_text());assert z['status']=='PASS_FULL_NEF_ROOTLESS_SOURCE_TRANSPORT';M=matrix(ZZ,z['transport']);G=matrix(ZZ,z['frame']);h=vector(ZZ,z['ample_class']);D=M.row(0)
assert abs(M.det())==1 and M*NS*M.transpose()==block_diagonal_matrix(U,-G)
assert G.is_positive_definite() and G.det()==1020 and G.__pari__().qfminim(2,10000,2)[0]==0
f=vector(ZZ,[1]+[0]*18);O=vector(ZZ,[-1,1]+[0]*17);root=F[:16,:16];aa=root.inverse()*vector(ZZ,[1]*16)
assert h==vector(ZZ,[22,20]+list(-2*aa)+[0]) and h*NS*h==514 and h*NS*O==2
assert list(-h[2:18])==[2,2,2,2,3,4,3,9,16,21,24,25,24,21,16,9]
curves=[O];offset=0
for n in [3,3,4,10]:
 total=vector(ZZ,19)
 for j in range(n-1):
  r=vector(ZZ,19);r[2+offset+j]=1;curves.append(r);total+=r
 curves.append(f-total);offset+=n-1
assert all(r*NS*r==-2 and h*NS*r>0 and D*NS*r>=0 for r in curves)
# Written ample bound: for an irreducible horizontal C other than O,
# component multiplicities sum to C.f in each I_n, and O.C>=0.
# Hence h.C>= (42-(2+2+4+25))*C.f =9*C.f. All vertical components are positive.
assert 42-(2+2+4+25)==9
assert D*NS*D==0 and D[1]>0
bound=int((D*NS*h-1)//9);b=D[1];w=vector(ZZ,D[2:]);tested=0;enumerated=0
for d in range(1,bound+1):
 vcol=F*w;Q=block_matrix(ZZ,[[b*b*F,(-b*d*vcol).column()],[(-b*d*vcol).row(),matrix(ZZ,1,1,[d*d*(w*F*w)+2*b*b])]])
 L=Q.LLL_gram().transpose();Qred=L*Q*L.transpose();raw=Qred.__pari__().qfminim(4*b*b-1,10000,2);assert raw[0]<20000
 for rr in matrix(ZZ,raw[2]).columns():
  vv=rr*L
  if abs(vv[-1])!=1:continue
  vv*=vv[-1];v=vector(ZZ,vv[:17]);norm=v*F*v
  if (norm-2)%(2*d):continue
  r=vector(ZZ,[(norm-2)//(2*d),d]+list(v));assert r*NS*r==-2
  enumerated+=1;assert D*NS*r>=0
 tested+=1
Z=M.row(1)-D;assert Z*NS*Z==-2 and Z*NS*D==1
sections=[]
for i in range(17):
 S=(G[i,i]//2-1)*D+M.row(1)+M.row(i+2)
 assert S*NS*S==-2 and S*NS*D==1
 sections.append(S)
# Shioda projection with no reducible fibers is exactly the integral frame basis.
projections=[S-Z-(S*NS*Z+2)*D for S in sections]
assert matrix(ZZ,projections)==M[2:,:] and -(matrix(ZZ,projections)*NS*matrix(ZZ,projections).transpose())==G
out={'status':'PASS_PHYSICAL_NEF_MW17_CLASSES','source_claim':'EC-K3-DET1020-EXPLICIT-RATIONAL-SOURCE-20260915','basis':json.loads((P/'physical-source-basis.json').read_text())['basis'],'transport':[[int(c) for c in r] for r in M.rows()],'fiber':list(map(int,D)),'zero':list(map(int,Z)),'section_classes':[[int(c) for c in S] for S in sections],'height_gram':[[int(c) for c in r] for r in G.rows()],'MW_rank':17,'MW_torsion':0,'height_determinant':1020,'source_fiber_degree':int(D[1]),'full_nef_check':{'ample_class':list(map(int,h)),'ample_square':514,'ample_degree':int(D*NS*h),'horizontal_degrees_checked':tested,'upper_degree_bound':bound,'negative_roots_found':enumerated,'method':'Positive-definite augmented-form enumeration of all roots negative on D in the complete ample-bounded degree range, plus O and all vertical components.'},'written_inputs':['Source full rational NS and geometric Picard19 certificate','Ample criterion and K3 Riemann-Roch; an effective irreducible negative wall has square-2','Root reflections preserve the positive cone, bounding h.C<h.D for a negative wall','Primitive nef square-zero K3 divisor gives an elliptic pencil; rootless frame gives irreducible fibers','Degree-one square-2 effective classes on a rootless elliptic K3 are sections; the unimodular U gives integral MW saturation'],'boundary':'Actual nef pencil and saturated rational section DIVISOR CLASSES. Rational functions, Weierstrass equation and section coordinates are not computed.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();target=P/'physical-mw17-certificate.json'
if args.check:assert json.loads(target.read_text())==out
else:target.write_text(json.dumps(out,indent=2)+'\n')
print(out['status'],'degree',out['source_fiber_degree'],'nef degrees checked',tested)
