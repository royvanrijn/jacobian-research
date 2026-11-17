#!/usr/bin/env python3
"""Primary cubic carrier rank-drop queries, retaining every projective lambda chart."""
from pathlib import Path
from collections import defaultdict
import hashlib,json,subprocess,time
import sympy as S
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
W=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(folder):
 ip=folder/'input.json';d=json.loads(ip.read_text());Z=S.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],Z));gens=[sum(S.Rational(c)*S.prod(z**k for z,k in zip(Z,e) if k) for e,c in terms) for terms in d['high_equations']];groups=defaultdict(list)
 for i in range(3,12):groups[W[i]].append(i)
 branches=[];l1,l2=S.symbols('lam1 lam2')
 for weight,rows in groups.items():
  if len(rows)==1:continue
  for a,pivot in enumerate(rows[:-1]):
   later=rows[a+1:];lam={pivot:S.Integer(1)}|dict(zip(later,[l1,l2]));coeff={}
   for atom in d['atoms']:
    if atom['row'] in lam:
     e=tuple(atom['exponent']);coeff[e]=coeff.get(e,0)+lam[atom['row']]*S.Rational(atom['base_coefficient'])*by[atom['parameter']]
   branches.append(dict(weight=weight,pivot=pivot+1,rows=[i+1 for i in rows],later_rows=[i+1 for i in later],equations=[S.expand(v) for v in coeff.values() if v]))
 lines=['ring r=0,('+','.join(map(str,Z))+',lam1,lam2),dp;','short=0;','ideal I='+','.join(map(sing,gens))+';','ideal B=std(I);','print("MODEL_DIM:"+string(dim(B)));']
 for j,b in enumerate(branches):
  lines+=['ideal J'+str(j)+'=B,'+','.join(map(sing,b['equations']))+';','ideal G'+str(j)+'=std(J'+str(j)+');','print("BRANCH:'+str(j)+':PIVOT:'+str(b['pivot'])+':DIM:"+string(dim(G'+str(j)+')));','if(G'+str(j)+'[1]!=1){print("PROPER:'+str(j)+'");print(G'+str(j)+');}']
 lines+=['quit;'];script=folder/'rank-model.sing';script.write_text('\n'.join(lines)+'\n');started=time.monotonic();print('RANK START',folder.name,'BRANCHES',len(branches),flush=True)
 try:
  p=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=p.stdout+p.stderr;status='primary-bases-complete' if p.returncode==0 and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (folder/'rank-model.log').write_text(text);payload=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(ip),branches=[{k:(list(map(str,v)) if k=='equations' else v) for k,v in b.items()} for b in branches],seconds=round(time.monotonic()-started,3),scope='Linear dependence of nonlinear carrier parts within same-weight blocks; charts normalize first nonzero lambda to1. Last-pivot and singleton charts are the separately retained identity-carrier queries. Only an extracted/replayed proper point would give a cubic dimension deletion; primary units are not independent certificates.');(folder/'rank-checkpoint.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('RANK',folder.name,status,text[:4000],flush=True)
def main():
 folders=[HERE/'carrier-reparam-full']+sorted((HERE/'carrier-order-search').glob('*'))
 for folder in folders:
  if (folder/'input.json').exists():run(folder)
if __name__=='__main__':main()
