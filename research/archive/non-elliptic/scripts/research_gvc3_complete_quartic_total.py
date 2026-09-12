#!/usr/bin/env python3
"""Modular Singular probe of the complete weighted-quartic GVC(3) jet.
Research-only: modular emptiness is not promoted automatically to a Q theorem.
"""
from __future__ import annotations
import argparse,hashlib,json,math,subprocess,tempfile,time
from collections import defaultdict
from pathlib import Path

NAMES=("a","u0","u1","u2","u3","u4","r0","r1","r2","r3","r4")
AT=(((1,1,1),None),((2,0,0),"a"),((0,4,0),"u0"),((0,3,1),"u1"),((0,2,2),"u2"),((0,1,3),"u3"),((0,0,4),"u4"))
PT=(((1,2,0),None),((2,0,0),None),((0,4,0),"r0"),((0,3,1),"r1"),((0,2,2),"r2"),((0,1,3),"r3"),((0,0,4),"r4"))
CHARTS={
 "r2":({"r2":1},(),"D(r2)/Gm"),"r3":({"r3":1},(),"D(r3)/Gm"),"r4":({"r4":1},(),"D(r4)/Gm"),
 "u0_r0":({"u0":1},("r0",),"D(u0*r0)/Gm"),"u0_r1":({"u0":1},("r1",),"D(u0*r1)/Gm"),
 "u1_r1":({"r1":1},("u1",),"D(u1*r1)/Gm")}

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def inc(a,i):
 b=list(a);b[i]+=1;return tuple(b)
def terms(ts,vals,idx):
 out=[]
 for e,n in ts:
  if n is None: out.append((e,1,None))
  elif n in vals:
   if vals[n]: out.append((e,int(vals[n]),None))
  else: out.append((e,1,idx[n]))
 return out

def power(ts,m,nv,p):
 z=(0,)*nv;cur={(0,0,0):{z:1}}
 for _ in range(m):
  nxt=defaultdict(lambda:defaultdict(int))
  for e0,q in cur.items():
   for e1,c,i in ts:
    t=nxt[add(e0,e1)]
    for pe,pc in q.items(): t[pe if i is None else inc(pe,i)]=(t[pe if i is None else inc(pe,i)]+pc*c)%p
  cur={e:{x:c for x,c in q.items() if c} for e,q in nxt.items()}
 return cur

def moment(m,vals,p):
 ns=tuple(n for n in NAMES if n not in vals);idx={n:i for i,n in enumerate(ns)}
 A=power(terms(AT,vals,idx),m,len(ns),p);P=power(terms(PT,vals,idx),m,len(ns),p);out=defaultdict(int)
 for e,qa in A.items():
  qp=P.get(e)
  if not qp: continue
  w=math.factorial(e[0])*math.factorial(e[1])*math.factorial(e[2])%p
  for a,ca in qa.items():
   for b,cb in qp.items(): out[add(a,b)]=(out[add(a,b)]+ca*cb*w)%p
 return ns,{e:c for e,c in out.items() if c}

def mul(A,B,p):
 o=defaultdict(int)
 for a,ca in A.items():
  for b,cb in B.items(): o[add(a,b)]=(o[add(a,b)]+ca*cb)%p
 return {e:c for e,c in o.items() if c}
def subst(poly,i,repl,p):
 n=len(next(iter(poly),()));mp=max((e[i] for e in poly),default=0);pw=[{(0,)*(n-1):1}]
 for _ in range(mp): pw.append(mul(pw[-1],repl,p))
 o=defaultdict(int)
 for e,c in poly.items():
  rest=e[:i]+e[i+1:]
  for x,v in mul({rest:c},pw[e[i]],p).items(): o[x]=(o[x]+v)%p
 return {e:c for e,c in o.items() if c}
def reduced(vals,cut,p):
 ns,mu1=moment(1,vals,p);i=ns.index("a");ae=[0]*len(ns);ae[i]=1;ae=tuple(ae);ci=pow(mu1[ae],-1,p);rep={}
 for e,c in mu1.items():
  if e!=ae: rep[e[:i]+e[i+1:]]=(-c*ci)%p
 rn=ns[:i]+ns[i+1:];out=[]
 for m in range(2,cut+1): out.append(subst(moment(m,vals,p)[1],i,rep,p))
 return rn,out,rep

def mon(e,ns):
 f=[]
 for n,k in zip(ns,e):
  if k==1:f.append(n)
  elif k>1:f.append(f"{n}^{k}")
 return "*".join(f) or "1"
def expr(q,ns,p):
 if not q:return "0"
 inv=pow(next(iter(sorted(q.items())))[1],-1,p);a=[]
 for e,c in sorted(q.items(),reverse=True):
  c=c*inv%p;m=mon(e,ns);a.append(str(c) if m=="1" else m if c==1 else f"{c}*{m}")
 return "+".join(a)
def program(chart,cut,p,ns,qs,invert):
 iv=tuple(f"iv{i}" for i in range(len(invert)));vs=ns+iv;eq=[expr(q,ns,p) for q in qs if q]
 eq += [f"{x}*{y}-1" for x,y in zip(iv,invert)]
 return f'''option(redSB);\nring R={p},({','.join(vs)}),dp;\nideal I={',\n'.join(eq)};\ntimer=1; ideal G=slimgb(I); int tm=timer;\nprint("RESULT_BEGIN");\nprint("basis_size="+string(size(G))); print("dimension="+string(dim(G))); print("elapsed_ticks="+string(tm));\nif((size(G)==1)&&(G[1]==1)){{print("unit=1");}}else{{print("unit=0");}}\nprint("BASIS_BEGIN"); print(G); print("BASIS_END"); print("RESULT_END"); quit;\n'''
def parse(s):
 o={}
 for line in s.splitlines():
  line=line.strip()
  for k in ("basis_size","dimension","elapsed_ticks","unit"):
   if line.startswith(k+"="):o[k]=int(line[len(k)+1:])
 if "BASIS_BEGIN" in s and "BASIS_END" in s:
  b=s.split("BASIS_BEGIN",1)[1].split("BASIS_END",1)[0].strip();o["basis_sha256"]=hashlib.sha256(b.encode()).hexdigest();o["basis_preview"]=b[:4000]
 return o

def case(singular,name,cut,p,timeout,wd):
 vals,inv,cov=CHARTS[name];t=time.monotonic();ns,qs,rep=reduced(vals,cut,p);gen=time.monotonic()-t;src=program(name,cut,p,ns,qs,inv);fn=wd/f"{name}_{p}_{cut}.sing";fn.write_text(src)
 o={"chart":name,"coverage":cov,"values":vals,"invert":inv,"prime":p,"cutoff":cut,"variables":ns,"term_counts":list(map(len,qs)),"generation_seconds":gen,"a_replacement_terms":len(rep),"input_sha256":hashlib.sha256(src.encode()).hexdigest()}
 try:
  t=time.monotonic();r=subprocess.run([singular,"-q",str(fn)],text=True,capture_output=True,timeout=timeout);o["singular_seconds"]=time.monotonic()-t;o["returncode"]=r.returncode;o.update(parse(r.stdout));o["stderr"]=r.stderr[-4000:];o["stdout_tail"]=r.stdout[-4000:];o["status"]="completed" if r.returncode==0 else "failed"
 except subprocess.TimeoutExpired as x:o["singular_seconds"]=time.monotonic()-t;o["status"]="timeout";o["stderr"]="";o["stdout_tail"]=""
 return o

def main():
 ap=argparse.ArgumentParser();ap.add_argument("--singular",default="Singular");ap.add_argument("--primes",type=int,nargs="+",default=[101,103,107]);ap.add_argument("--cutoffs",type=int,nargs="+",default=[6,7]);ap.add_argument("--charts",nargs="+",default=list(CHARTS));ap.add_argument("--timeout",type=int,default=1200);ap.add_argument("--output",type=Path,default=Path("artifacts/generated-results/gvc3_complete_quartic_total_probe.json"));a=ap.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
 art={"format":"gvc3-complete-quartic-total-probe-v1","status":"exact modular research computation; not a characteristic-zero theorem","model":{"operator":"dz*dL*dM+a*dz^2+sum u_i*dL^(4-i)*dM^i","polynomial":"z*L^2+z^2+sum r_i*L^(4-i)*M^i","grading":"2z+L+M=4","support_target":"(a,r2,r3,r4,u0*r0,u0*r1,u1*r1)"},"primes":a.primes,"cutoffs":a.cutoffs,"cases":[]}
 with tempfile.TemporaryDirectory() as d:
  for p in a.primes:
   for n in a.charts:
    for c in a.cutoffs:
     print("RUN",n,p,c,flush=True);r=case(a.singular,n,c,p,a.timeout,Path(d));art["cases"].append(r);a.output.write_text(json.dumps(art,indent=2,sort_keys=True)+"\n");print("DONE",n,p,c,r.get("status"),r.get("unit"),r.get("dimension"),r.get("basis_size"),flush=True)
 art["summary"]={str(c):{"unit":sum(x.get("cutoff")==c and x.get("unit")==1 for x in art["cases"]),"nonunit":sum(x.get("cutoff")==c and x.get("unit")==0 for x in art["cases"]),"timeout":sum(x.get("cutoff")==c and x.get("status")=="timeout" for x in art["cases"])} for c in a.cutoffs};a.output.write_text(json.dumps(art,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
