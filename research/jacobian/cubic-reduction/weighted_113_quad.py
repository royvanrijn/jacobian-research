#!/usr/bin/env python3
"""Clean parameter units before the remaining quadratic-graph field/ring solve."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 mode=sys.argv[1];assert mode in ['field','clean'];source=HERE/'weighted-113/boundary/both-quad/reduced/input.json';d=json.loads(source.read_text());N=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],N));P=by['c2'];inv=by['branch_inv'];Z=tuple(by[name] for name in d['free_tail_variables']);clean=[];powers=[]
 for terms in d['equations'][:-1]:
  assert all(e[d['variables'].index('branch_inv')]==0 for e,c in terms);k=min(e[d['variables'].index('c2')] for e,c in terms);f=sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(N,e)) for e,c in terms));clean.append(sp.expand(f/P**k));powers.append(k)
 out=source.parent.parent/mode;out.mkdir(exist_ok=True);V=(P,)+Z if mode=='field' else Z+(P,inv);eqs=clean if mode=='field' else clean+[1-inv*P];payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=['c2'] if mode=='field' else [],variables=list(map(str,Z if mode=='field' else V)),removed_parameter_powers=powers,equations=[enc(f,V) for f in eqs],scope='Same complete quadratic-graph boundary reduction onP!=0. Each retained generator is divided only by its recorded monomialP power; saturation uses1-inv*P. This is equivalent on that branch. A field unit still requires every denominator root; proper ideals require actual candidate extraction and full replay.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r='+('(0,c2)' if mode=='field' else '0')+',('+','.join(map(str,Z if mode=='field' else V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-changed-quadratic-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:5000],flush=True)
if __name__=='__main__':main()
