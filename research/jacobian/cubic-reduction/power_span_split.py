#!/usr/bin/env python3
"""Certify the cube-locus classification by rational linear-factor splitting trees."""
from pathlib import Path
import hashlib,json,time
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(B):
 if not B.cols:return B
 R,_=B.T.rref();rows=[list(R.row(i)) for i in range(R.rows) if any(R.row(i))]
 return sp.Matrix(rows).T if rows else sp.zeros(B.rows,0)
def restriction(terms,B,Y):
 lin=[sum(B[i,j]*Y[j] for j in range(B.cols)) for i in range(B.rows)]
 return sp.expand(sum(sp.Rational(c)*sp.prod(f**a for f,a in zip(lin,e)) for e,c in terms))
def main():
 out=HERE/"power-span-split";out.mkdir(exist_ok=True);logs=[];started=time.monotonic()
 for case in range(15):
  source=HERE/"power-single-compress"/("case-%02d.json"%case);data=json.loads(source.read_text());m=len(data["support"]);constraints=data["constraints"]
  nodes=[];memo={};fallbacks=0
  def visit(B):
   nonlocal fallbacks
   B=canonical(B);key=tuple(tuple(str(a) for a in row) for row in B.tolist())
   if key in memo:return memo[key]
   index=len(nodes);memo[key]=index;node=dict(basis=[[str(a) for a in row] for row in B.tolist()]);nodes.append(node)
   if B.cols==0:node["kind"]="leaf";return index
   Y=sp.symbols("v1:"+str(B.cols+1));polys=[restriction(f,B,Y) for f in constraints]
   if not any(polys):node["kind"]="leaf";return index
   choices=[]
   def consider(f,combo):
    if not f:return
    coefficient,factors=sp.factor_list(f,*Y)
    if all(sp.Poly(g,*Y).total_degree()==1 and sp.Poly(g,*Y).coeff_monomial((0,)*len(Y))==0 for g,power in factors):
     choices.append((len(factors),sum(len(sp.Poly(g,*Y).terms()) for g,power in factors),coefficient,factors,combo))
   for i,f in enumerate(polys):consider(f,[[i,"1"]])
   if not choices:
    fallbacks+=1;mons=sorted(set(e for f in polys if f for e,c in sp.Poly(f,*Y).terms()))
    M=sp.Matrix([[sp.Poly(f,*Y).coeff_monomial(e) for e in mons] for f in polys]);R,_=M.row_join(sp.eye(len(polys))).rref()
    for i in range(R.rows):
     f=sum(R[i,j]*sp.prod(v**a for v,a in zip(Y,e)) for j,e in enumerate(mons))
     consider(f,[[j,str(R[i,len(mons)+j])] for j in range(len(polys)) if R[i,len(mons)+j]])
   if not choices:raise RuntimeError('No complete linear factor split in case%d dimension%d'%(case,B.cols))
   _,_,coefficient,factors,combo=min(choices,key=lambda x:(x[0],x[1]));node.update(kind="split",combination=combo,scalar=str(coefficient),factors=[],children=[])
   for f,power in factors:
    row=[sp.expand(f).coeff(y) for y in Y];kernel=sp.Matrix([row]).nullspace();K=sp.Matrix.hstack(*kernel) if kernel else sp.zeros(B.cols,0)
    node["factors"].append(dict(row=list(map(str,row)),power=int(power)));node["children"].append(visit(B*K))
   return index
  root=visit(sp.eye(m));payload=dict(case=case,root=root,source_sha256=sha(Path(__file__)),case_sha256=sha(source),nodes=nodes,
    scope="Exact factor-splitting tree; independent coefficient replay is required before classifying all cube points.")
  cp=out/("case-%02d.json"%case);cp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
  row=dict(case=case,nodes=len(nodes),leaves=sum(n['kind']=='leaf' for n in nodes),fallbacks=fallbacks,
    leaf_dimensions=sorted({len(n['basis'][0]) for n in nodes if n['kind']=='leaf'}),tree_sha256=sha(cp),seconds=round(time.monotonic()-started,3));logs.append(row);print(json.dumps(row),flush=True)
  (out/'checkpoint.json').write_text(json.dumps(dict(status='running',source_sha256=sha(Path(__file__)),cases=logs),indent=2,sort_keys=True)+'\n')
 receipt=dict(status='primary-trees-complete',source_sha256=sha(Path(__file__)),cases=logs,scope='Rational factor trees prove a finite linear-subspace cover after independent replay; no Keller dimension minimum follows directly.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
