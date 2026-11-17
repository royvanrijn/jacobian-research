#!/usr/bin/env python3
"""Export exact certificates for raw13-pattern rigidity in the specified collision family."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expression(terms,V):return s.Add(*(s.Rational(c)*s.prod(x**a for x,a in zip(V,e)) for e,c in terms))
def execute(out,tag,V,generators,targets):
 lines=['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in generators)+';','ideal G=std(I); matrix M; matrix P; int h;']
 for i,target in enumerate(targets):
  lines+=['if(reduce('+sing(target)+',G)!=0){print("ERROR:MEMBERSHIP");quit;}',
   'M=lift(I,ideal('+sing(target)+'));','P=matrix(I)*M;','if(P[1,1]!='+sing(target)+'){print("ERROR:IDENTITY");quit;}',
   'print("CERT:%d:%d");'%(i,len(generators)),'for(h=1;h<=nrows(M);h++){print(string(M[h,1]));}']
 lines+=['quit;'];script=out/(tag+'.sing');script.write_text('\n'.join(lines)+'\n');r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=30);text=r.stdout+r.stderr;(out/(tag+'.log')).write_text(text);assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000]
 raw=text.splitlines();pos=0;certs=[]
 for i,target in enumerate(targets):
  assert raw[pos]=='CERT:%d:%d'%(i,len(generators));pos+=1;mult=[]
  for line in raw[pos:pos+len(generators)]:
   f=s.Poly(s.sympify(line.replace('^','**'),locals={str(v):v for v in V}),*V);mult.append([[list(e),str(c)] for e,c in f.terms() if c])
  pos+=len(generators);certs.append(dict(target=[[list(e),str(c)] for e,c in s.Poly(target,*V).terms() if c],multipliers=mult))
 assert pos==len(raw)
 payload=dict(source_sha256=sha(Path(__file__)),singular_source_sha256=sha(script),variables=list(map(str,V)),certificates=certs)
 cp=out/(tag+'.json');cp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');return sha(cp)
def main():
 out=HERE/'raw13-rigidity';out.mkdir(exist_ok=True);rawpath=HERE/'raw13-specialization/input.json';normalpath=HERE/'raw13-normalized/input.json';raw=json.loads(rawpath.read_text());normal=json.loads(normalpath.read_text());started=time.monotonic();hashes={}
 V=s.symbols(' '.join(raw['variables']));gens=[expression(f,V) for f in raw['equations']]
 for i in [0,1]:hashes['raw-unit-%d'%i]=execute(out,'raw-unit-%d'%i,V,gens+[V[i]],[s.Integer(1)])
 N=s.symbols(' '.join(normal['variables']));ngens=[expression(f,N) for f in normal['equations']];baseline=[s.Rational(atom['coefficient']) for atom in raw['atoms'][2:]]+[s.Integer(1),s.Integer(1)];assert len(N)==len(baseline)==13
 assert all(s.expand(f.subs(dict(zip(N,baseline))))==0 for f in ngens)
 hashes['normalized-squares']=execute(out,'normalized-squares',N,ngens,[(x-a)**2 for x,a in zip(N,baseline)])
 receipt=dict(status='exact-primary-identities',source_sha256=sha(Path(__file__)),raw_input_sha256=sha(rawpath),normalized_input_sha256=sha(normalpath),certificates=hashes,identities=15,seconds=round(time.monotonic()-started,3),scope='Two unit identities exclude zero first-output coefficients; thirteen squared-coordinate memberships force the normalized13-atom coefficient vector to the known representative over fields. Independent reconstruction and replay required.')
 (out/'certification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
