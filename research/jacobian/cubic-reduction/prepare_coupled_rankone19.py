#!/usr/bin/env python3
"""Coupled rank6 companion: z=-4/3, freew, with grouped and rotated cube forms."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,Counter
import hashlib,json,random,time
import sympy as S
from flint import fmpq_mat
from prepare_companion_pencil19 import cubes,decompose,primitive,encode
HERE=Path(__file__).resolve().parent
N=19
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def point_score(rows):
 roots=Counter();generic=sum(bool(row) for row in rows)
 for row in rows:
  ids=[j for j in range(19) if row.get(j+19)]
  if not ids:continue
  j=ids[0];root=-row.get(j,Q(0))/row[j+19]
  if all(row.get(j,Q(0))+root*row.get(j+19,Q(0))==0 for j in range(19)):roots[root]+=1
 if not roots:return generic,Q(0)
 root=min(roots,key=lambda r:(-roots[r],abs(r.numerator)+r.denominator,r));return generic-roots[root],root

def main():
 start=time.monotonic();out=HERE/'coupled-rankone19';out.mkdir(exist_ok=True);base=HERE/'cubic12-integer-witness.json';bd=json.loads(base.read_text());X=S.symbols('x1:20');w=S.Symbol('w');z=S.Rational(-4,3);a,b,c,R=X[:4];T=X[9];SS=X[6];U=X[10];tau=X[18];K=[sum(S.Rational(t['coefficient'])*S.prod(x**k for x,k in zip(X,t['exponent'])) for t in row) for row in bd['coefficients']];E=list(X[:12]);E[7]+=w*c*R;E[11]+=w*R**2;E[8]+=z*c*R;E[10]+=z*R*(3*b-SS);G=[S.expand(f.subs(dict(zip(X[:12],E)),simultaneous=True)) for f in K];assert all(S.Poly(f,*X).total_degree()<=3 for f in G);Q2=[];C3=[]
 for f in G:
  terms=S.Poly(f,*X).terms();Q2.append(sum(co*S.prod(x**k for x,k in zip(X,e)) for e,co in terms if sum(e)==2));C3.append(sum(co*S.prod(x**k for x,k in zip(X,e)) for e,co in terms if sum(e)==3))
 Cstar=C3[:5]+[a*c*R];iota=S.zeros(12,6)
 for i in range(5):iota[i,i]=1
 iota[9,5]=w;iota[10,5]=z;assert list(iota*S.Matrix(Cstar))==C3
 H=[tau*Q2[i]-tau**2*sum(iota[i,j]*X[12+j] for j in range(6)) for i in range(12)]+Cstar+[tau**3];H=list(map(S.expand,H));H0=[S.expand(f.subs(w,0)) for f in H];H1=list(map(lambda f:S.expand(S.diff(f,w)),H));assert all(S.expand(f-g-w*h)==0 for f,g,h in zip(H,H0,H1));P0=[{e:Q(c0) for e,c0 in S.Poly(f,*X).terms() if c0} for f in H0];P1=[{e:Q(c0) for e,c0 in S.Poly(f,*X).terms() if c0} for f in H1]
 # Group the rank-one cubic exactly before selecting polarization formulas.
 dd,ee,hh,gg=S.symbols('D E H G');aliases=[3*b-SS,b+T,3*b-2*SS,3*b-2*U];Y=[dd,ee,hh,gg];group=list(H);wanted={12:-a**2*dd+a*c*X[5],13:12*a*dd*ee+12*a*R*X[7]-12*a*c*X[11]+12*w*a*c*R-12*w*R**2*SS,14:2*a*c*gg+4*a*X[8]*dd+8*a*R*X[4]-S.Rational(4,3)*R*hh**2-S.Rational(16,3)*c*R**2,11:a*ee*tau+w*R**2*tau,10:a*X[8]*tau+z*R*dd*tau-z*tau**2*X[17]}
 for i,term in wanted.items():assert S.expand(H[i]-term.subs(dict(zip(Y,aliases)),simultaneous=True))==0;group[i]=term
 vectors=[tuple(int(i==j) for i in range(N)) for j in range(N)]+[tuple(Q(S.expand(f).coeff(x)) for x in X) for f in aliases];terms=[];options=[];pool=set();parents={}
 for tag in ['companion-pencil19','companion-pencil19-mixed','companion-w19']:
  ip=HERE/tag/'input.json';d=json.loads(ip.read_text());parents[tag]=sha(ip);pool.update(map(tuple,d['forms']))
 reference=list(map(tuple,json.loads((HERE/'companion-pencil19/input.json').read_text())['forms'][:89]));L=S.eye(N);L[1,6]=-S.Rational(2,3);L[6,6]=-1;L[9,6]=S.Rational(2,3);L[10,1]=S.Rational(8,3);L[10,6]=-S.Rational(8,9);L[10,10]=-1;assert L*L==S.eye(N);L2=S.Matrix(L);L2[2,2]=-1;L2[8,2]=S.Rational(1,9);L2[11,11]=-1;L2[11,3]=-S.Rational(2,3);assert L2*L2==S.eye(N)
 for M in [L,L2]:
  for f in reference:pool.add(primitive((S.Matrix(1,N,f)*M).tolist()[0])[0])
 for row,f in enumerate(group):
  for exponent,co in S.Poly(f,*(X+tuple(Y))).terms():
   if not co:continue
   factors=[v for v,k in zip(vectors,exponent) for _ in range(k)];assert len(factors)==3;pco=S.Poly(co,w);assert pco.degree()<=1;coef=[Q(pco.nth(0)),Q(pco.nth(1))];terms.append(dict(row=row,coefficient0=str(coef[0]),coefficient1=str(coef[1]),factors=[[str(a) for a in v] for v in factors]));ops=[decompose(factors,k) for k in range(5)] if len(set(factors))==3 else [decompose(factors,0)];options.append(ops)
   for op in ops:pool.update(op)
 rng=random.Random(261123);best=999;bestroot=Q(0);bestchoices=None;bestrows=None
 for restart in range(40):
  choices=[0 if restart==0 else rng.randrange(len(op)) for op in options];used=Counter(v for ops,ch in zip(options,choices) for v in ops[ch])
  for sweep in range(12):
   change=False;order=list(range(len(terms)));rng.shuffle(order)
   for i in order:
    for v in options[i][choices[i]]:used[v]-=1
    costs=[sum(used[v]==0 for v in op) for op in options[i]];which=[j for j,cost in enumerate(costs) if cost==min(costs)];pick=choices[i] if choices[i] in which else rng.choice(which);change|=pick!=choices[i];choices[i]=pick
    for v in options[i][pick]:used[v]+=1
   if not change:break
  rows={}
  for term,ops,ch in zip(terms,options,choices):
   for form,co in ops[ch].items():
    col=rows.setdefault(form,{})
    for part in range(2):
     j=term['row']+19*part;v=Q(term['coefficient'+str(part)])*co
     if v:col[j]=col.get(j,Q(0))+v
  rows={f:{j:v for j,v in col.items() if v} for f,col in rows.items()};rows={f:col for f,col in rows.items() if col};score,root=point_score(list(rows.values()))
  if score<best:best=score;bestroot=root;bestchoices=choices;bestrows=rows;print('TERMWISE',best,'w',root,'restart',restart,flush=True)
 forms=sorted(pool);positions={f:i for i,f in enumerate(forms)};seed=[positions[f] for f,col in bestrows.items() if any(col.get(j,Q(0))+bestroot*col.get(j+19,Q(0)) for j in range(19))];priority=seed+[positions[f] for f in bestrows if positions[f] not in seed]+[i for i in range(len(forms)) if forms[i] not in bestrows];polys=[cubes(f) for f in forms];mons=sorted(set().union(*(set(p) for p in polys+P0+P1)));m=len(forms);C=fmpq_mat([[str(p.get(e,0)) for p in polys] for e in mons]);TT=fmpq_mat([[str(p.get(e,0)) for p in P0+P1] for e in mons]);ordered=fmpq_mat([[str(C[i,j]) for j in priority] for i in range(len(mons))]);rr,r=ordered.rref();basis=[priority[next(j for j in range(m) if rr[i,j])] for i in range(r)];joint=fmpq_mat([[str(C[i,j]) for j in basis]+[str(C[i,j]) for j in range(m)]+[str(TT[i,j]) for j in range(38)] for i in range(len(mons))]);rr,jr=joint.rref();assert jr==r;A=[{str(j):str(rr[i,r+j]) for j in range(m) if rr[i,r+j]} for i in range(r)];B=[{str(j):str(rr[i,r+m+j]) for j in range(38) if rr[i,r+m+j]} for i in range(r)];padding=positions[(0,)*18+(1,)];assert all(Q(row.get('18','0'))==int(basis[i]==padding) and not Q(row.get('37','0')) for i,row in enumerate(B));actual,root=point_score([{int(k):Q(c0) for k,c0 in row.items()} for row in B]);assert actual<=best
 payload=dict(source_sha256=sha(Path(__file__)),base_sha256=sha(base),parent_inputs=parents,variables=list(map(str,X)),dimension=19,parameter='w',fixed_z='-4/3',excluded_parameter_values=[],source_shear=list(map(str,E)),base_family=list(map(str,G)),embedding=[[str(v) for v in row] for row in iota.tolist()],cubic_star=list(map(str,Cstar)),H0=[encode(p) for p in P0],H1=[encode(p) for p in P1],forms=forms,padding_index=padding,rank=r,monomials=len(mons),grouped_terms=terms,grouped_choices=bestchoices,termwise_best=best,termwise_parameter=str(bestroot),initial_basis_best=actual,initial_basis_parameter=str(root),rotation_matrices=[[str(v) for v in row] for row in L.tolist()]+[[str(v) for v in row] for row in L2.tolist()],scope='Coupled w,z rank6 normalized companion, z=-4/3, w arbitrary including0. Extra cubic outputacR, embeddingw*e10+z*e11. Grouped repeated factor and two rotated old89-form pools. Exact primary finite-pool matrices only; no improved Keller witness without independent full cube/determinant/collision verification.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');(out/'basis.json').write_text(json.dumps(dict(input_sha256=sha(out/'input.json'),basis=basis,A=A,B=B),indent=2,sort_keys=True)+'\n');print('PREPARED',m,'RANK',r,'MONS',len(mons),'INITIAL',actual,'w',root,'SECONDS',round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
