#!/usr/bin/env python3
"""Construct all six optimal linear/nonlinear output-frame tradeoffs."""
from pathlib import Path
import json,hashlib
import sympy as sp
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def encode(F,X):return [[{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*X).terms())] for p in F]
def main():
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());X=sp.symbols('x1:13');F=[sum(sp.Rational(t['coefficient'])*sp.prod(v**k for v,k in zip(X,t['exponent'])) for t in row) for row in raw['coefficients']];G=list(F);M=sp.eye(12);stages=[]
 def record():
  rows=[sp.Poly(p,*X) for p in G];nonlinear=[sum(sum(e)>=2 for e,c in p.terms()) for p in rows];linear=[sum(sum(e)==1 for e,c in p.terms()) for p in rows]
  stages.append({'stage':len(stages),'linear_matrix':[[str(v) for v in M.row(i)] for i in range(12)],'coefficients':encode(G,X),'linear_counts':linear,'nonlinear_counts':nonlinear,'linear_total':sum(linear),'nonlinear_total':sum(nonlinear),'total':sum(linear)+sum(nonlinear)})
 record();operations=[(1,6,-3),(1,9,-12),(2,7,-8),(8,7,3),(2,8,-4)]
 for i,j,c in operations:
  G[i]=sp.expand(G[i]+c*G[j]);M[i,:]=M[i,:]+c*M[j,:];record()
 assert [(d['linear_total'],d['nonlinear_total']) for d in stages]==[(12+k,36-k) for k in range(6)] and M.det()==1
 points=[[sp.Rational(v) for v in p] for p in raw['points']];images=[[p.subs(dict(zip(X,point))) for p in G] for point in points];assert images[0]==images[1]
 nonlinear_columns=set().union(*(set(e for e,c in sp.Poly(p,*X).terms() if sum(e)>=2) for p in F));assert len(nonlinear_columns)==30
 result={'status':'constructed-independent-replay-required','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'variables':raw['variables'],'dimension':12,'degree':3,'stages':stages,'operations':[{'target':i+1,'from':j+1,'coefficient':str(c)} for i,j,c in operations],'points':raw['points'],'common_image':[str(v) for v in images[0]],'nonlinear_source_columns':30,'output_block_rows':[[1],[2,7,10],[3,8,9],[4,12],[5],[6],[11]],'total_block_lower_bounds':[5,15,16,6,2,2,2],'scope':'Fixed source coordinates; arbitrary invertible linear output frames. Final map has31 nonlinear and17 linear occurrences,48 total, nonidentity unipotent linear part. Not an identity-linear-part36 improvement or global optimum.'}
 Minv=M.inv();sub=dict(zip(X,list(Minv*sp.Matrix(X))));C=[sp.expand(p.subs(sub,simultaneous=True)) for p in G];Ccounts=[sum(sum(e)>=2 for e,c in sp.Poly(p,*X).terms()) for p in C]
 result['identity_normalized_conjugate']={'inverse_linear_matrix':[[str(v) for v in Minv.row(i)] for i in range(12)],'coefficients':encode(C,X),'nonlinear_counts':Ccounts,'nonlinear_total':sum(Ccounts),'points':[[str(v) for v in M*sp.Matrix(point)] for point in points],'common_image':[str(v) for v in images[0]],'scope':'One exact source normalization M*K*M^-1. No global affine-orbit claim.'}
 print('IDENTITY NORMALIZED CONJUGATE',sum(Ccounts),Ccounts,flush=True)
 folder=HERE/'output-frame31';folder.mkdir(exist_ok=True);(folder/'witness.json').write_text(json.dumps(result,indent=2)+'\n');print('CONSTRUCTED',[(d['linear_total'],d['nonlinear_total']) for d in stages],'total48; independent replay required',flush=True)
if __name__=='__main__':main()
