#!/usr/bin/env python3
"""Enlarge the raw13 pattern to every weight-compatible degree-seven tail atom."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def main():
 out=HERE/'weighted-tail7';out.mkdir(exist_ok=True);u,w,v,beta,sat=s.symbols('u w v beta sat');f=1-3*u-w
 supports=[]
 for row in [1,2]:
  target=-row;atoms=[]
  for b in range(8):
   for c in range(8):
    a=target+b+2*c
    if a<0 or a+b+c>7 or a+b+c<2:continue
    atoms.append((b,c,(a,b,c)))
  supports.append(sorted(atoms,reverse=True))
 assert list(map(len,supports))==[8,9]
 G=s.symbols('g1:9');H=s.symbols('h1:10');g=u+sum(a*u**b*w**c for a,(b,c,e) in zip(G,supports[0]));h=w+sum(a*u**b*w**c for a,(b,c,e) in zip(H,supports[1]))
 J=s.expand(g*(s.diff(h,u)-3*s.diff(h,w))+2*h*(-s.diff(g,u)+3*s.diff(g,w))+f*(s.diff(g,u)*s.diff(h,w)-s.diff(g,w)*s.diff(h,u))-1)
 eqs=[s.expand(c) for e,c in s.Poly(J,u,w).terms() if c]
 boundaryg=s.Poly(s.expand(g.subs(w,1-3*u)),u);boundaryh=s.Poly(s.expand(h.subs(w,1-3*u)),u)
 boundary=[c for (e,),c in boundaryg.terms() if e>=2]+[c for (e,),c in boundaryh.terms() if e>=3]
 eqs+=boundary+[s.expand(g.subs({u:-beta,w:1+3*beta})+beta),1-sat*beta]
 V=tuple(G)+tuple(H)+(beta,sat);payload=dict(source_sha256=sha(Path(__file__)),variables=list(map(str,V)),tail_atoms=[dict(row=i+1,coefficient=str(a),exponent=list(e),quotient_exponent=[b,c]) for i,(coeff,atoms) in enumerate(zip([G,H],supports)) for a,(b,c,e) in zip(coeff,atoms)],
  equations=[encode(a,V) for a in eqs],boundary_equations=[encode(a,V) for a in boundary],first_output='a-3*a^2*b-a^3*c',scope='All degree-at-most7 tail monomials of weights -1,-2 for source weights(1,-1,-2), with fixed first output. Collision beta!=0. Boundary degree equations follow from G*Hprime-2*Gprime*H=1 on f=0.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(a) for a in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),tail_atom_counts=[8,9],nonlinear_atom_total=19,seconds=round(time.monotonic()-started,3),scope='Enlarged negative-weight tail model; proper basis needs parameter extraction and actual determinant/collision replay. No new sparsity or global minimum claim.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
