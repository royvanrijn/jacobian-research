#!/usr/bin/env python3
"""Exact controls for ordered one-vertex natural D5 lobe-tree reachability.

Fixed finite N9/N13 controller endpoints are substituted, not unrecovered
literal controller words replayed. The written proof supplies uniformity.
"""
from pathlib import Path
import argparse,hashlib,json,math,platform
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
import verify_uniform_hurwitz_attachment as h
import verify_d5_lobe_core_absorption as u

def power(p,k):
 if k<0:return power(h.inverse(p),-k)
 out=tuple(range(len(p)))
 while k:
  if k&1:out=h.compose(p,out)
  p=h.compose(p,p);k//=2
 return out

def inverse_word(word):return [-a for a in reversed(word)]
def replay(t,word):return u.replay(t,word)
def support(g):return {i for i,a in enumerate(g)if a!=i}

def root_pair(pair,v):
 for j in range(-2,3):
  new=replay(pair,[1 if j>=0 else -1]*abs(j))
  if new[1][v]==v:
   a=new[0][v];b=new[1][a];c=new[0][b];d=new[1][c]
   assert len({v,a,b,c,d})==5
   n=len(new[0])
   assert new==(h.permutation(n,[(v,a),(b,c)]),h.permutation(n,[(a,b),(c,d)]))
   return new,j,(a,b,c,d)
 raise AssertionError("no rooted endpoint")

def local_controls():
 dts=[]
 for vals in __import__("itertools").combinations(range(5),4):
  a,b,c,d=vals
  for edges in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))):
   dts.append(h.permutation(5,edges))
 assert len(dts)==15
 count=0;choices={}
 for a in dts:
  for b in dts:
   if len(h.cycles(h.compose(a,b)))!=1:continue
   pair=(a,b);C=h.product(pair)
   assert replay(pair,[1]*5)==pair
   for v in range(5):
    new,j,path=root_pair(pair,v)
    assert h.product(new)==C
    assert len([k for k in range(-2,3)if replay(pair,[1 if k>=0 else -1]*abs(k))[1][v]==v])==1
    count+=1;choices[j]=choices.get(j,0)+1
 assert count==600
 return {"natural_pairs":120,"shared_sheet_cases":600,"signed_power_frequencies":choices}

def rotate_old(t,r,v):
 old=t[:r];P=h.product(old)
 order=len(P)
 assert len(h.cycles(P))==1+len(P)-(2*r+1)
 k=next(j for j in range(2*r+1)if power(P,j)[0]==v)
 if k>r:k-=2*r+1
 twist=list(range(1,r))*r
 word=(inverse_word(twist) if k>0 else twist)*abs(k)
 moved=replay(t,word)
 assert moved[r:]==t[r:]
 expected=tuple(h.conjugate(power(P,k),g)for g in old)
 assert moved[:r]==expected
 relabel=power(P,-k)
 moved=tuple(h.conjugate(relabel,g)for g in moved)
 assert moved[:r]==old and relabel[v]==0
 return moved,len(word),k

def standardize_pair(t,r):
 n=len(t[0]);pair=t[r:r+2]
 rooted,j,path=root_pair(pair,0)
 localcode=(r+1 if j>=0 else -(r+1))
 t=replay(t,[localcode]*abs(j))
 assert t[r:r+2]==rooted
 rename=list(range(n))
 for a,b in zip(path,range(2*r+1,2*r+5)):rename[a]=b
 assert sorted(rename)==list(range(n))
 t=tuple(h.conjugate(tuple(rename),g)for g in t)
 expected=tuple(g+tuple(range(2*r+5,n))for g in u.initial(r))
 assert t[:r+2]==expected
 return t,abs(j)

def controls_rotation():
 rows=[]
 for r in range(4,9):
  n=2*r+5
  for v in range(2*r+1):
   t=tuple(g+tuple(range(2*r+1,n))for g in u.wr(r))
   fresh=list(range(2*r+1,n))
   pair=(h.permutation(n,[(v,fresh[0]),(fresh[1],fresh[2])]),
         h.permutation(n,[(fresh[0],fresh[1]),(fresh[2],fresh[3])]))
   # Deliberately choose an internal old attachment in many representatives.
   j=(v+r)%5
   pair=replay(pair,[1]*j)
   total=t+pair
   moved,cost,k=rotate_old(total,r,v)
   normalized,local=standardize_pair(moved,r)
   assert normalized==u.initial(r)
   assert cost<=r*r*(r-1) and local<=2
   rows.append({"r":r,"v":v,"k":k,"rotation_moves":cost,"local_moves":local})
 return rows

def tree(t,mode):
 n=4*t+1
 first=(h.permutation(n,[(0,1),(2,3)]),h.permutation(n,[(1,2),(3,4)]))
 first=replay(first,[1]*(mode%5))
 out=first
 for i in range(1,t):
  old=4*i+1
  if mode==0:v=0
  elif mode==1:v=old-1
  elif mode==2:v=(3*i+2)%old
  else:v=(i*i+3*i+1)%old
  a,b,c,d=range(old,old+4)
  pair=(h.permutation(n,[(v,a),(b,c)]),h.permutation(n,[(a,b),(c,d)]))
  pair=replay(pair,[1]*((2*i+mode)%5))
  out+=pair
 return out

def base_normalize(total):
 n=len(total[0]);C=h.product(total[:4]);standard=u.wr(4)
 P=h.product(standard);q=list(range(n));a,b=8,C[0]
 for _ in range(9):
  q[a]=b;a=P[a];b=C[b]
 assert q[0]==0 and sorted(q)==list(range(n))
 target=tuple(h.conjugate(tuple(q),g+tuple(range(9,n)))for g in standard)
 assert h.product(target)==C
 total=target+total[4:]
 total=tuple(h.conjugate(h.inverse(tuple(q)),g)for g in total)
 assert total[:4]==tuple(g+tuple(range(9,n))for g in standard)
 return total

def full_control(t,mode):
 total=tree(t,mode);n=len(total[0])
 assert len(h.cycles(h.product(total)))==1
 total=base_normalize(total);cost=0;rootings=0
 core=tuple(h.permutation(13,cs)for cs in u.CORE_CYCLES)
 target,q=u.canonical_target(core)
 for r in range(4,2*t,2):
  oldN=2*r+1
  assert total[:r]==tuple(g+tuple(range(oldN,n))for g in u.wr(r))
  used=set.union(*(support(g)for g in total[r:r+2]))
  old=used&set(range(oldN))
  assert len(old)==1 and used-old==set(range(oldN,oldN+4))
  v=next(iter(old))
  total,cc,k=rotate_old(total,r,v)
  total,jj=standardize_pair(total,r)
  cost+=cc+jj;rootings+=jj
  total=replay(total,u.prefix_word(r));em=u.core_embedding(r)
  assert total[:6]==tuple(u.embed(g,em,n)for g in core)
  total=tuple(u.embed(g,em,n)for g in target)+total[6:]
  total=replay(total,u.final_word(r))
  active=tuple(g[:oldN+4]for g in total[:r+2])
  rename=u.normalize_star(active,2*r+1)+tuple(range(oldN+4,n))
  total=tuple(h.conjugate(rename,g)for g in total)
  cost+=14*r-38
 assert total==u.wr(2*t)
 bound=sum((2*j)**2*(2*j-1)+14*(2*j)-38+2 for j in range(2,t))
 closed=2*t*t*(t-1)**2-(2*t*(t-1)*(2*t-1))//3+14*t*t-50*t+40
 assert bound==closed and cost<=bound
 return {"lobes":t,"mode":mode,"noncontroller_moves_replayed":cost,
         "noncontroller_upper_bound":bound,"local_pair_normalization_moves":rootings,
         "canonical_endpoint":True}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path,default=HERE/"d5_lobe_tree_reachability_verification.json")
 args=parser.parse_args()
 loc=local_controls();rot=controls_rotation()
 cases=[full_control(t,m)for t in (2,3,4,5,8,12)for m in range(4)]
 imported=json.loads((HERE/"d5_lobe_core_absorption_verification.json").read_text())
 assert imported["existential_literal_controller_word_bound"]==446265984
 assert imported["degree_nine_base_control"]["literal_controller_bound"]==6656
 paths=[Path(__file__).resolve(),HERE/"NATURAL_D5_LOBE_TREE_HURWITZ_REACHABILITY.md",
        HERE/"verify_uniform_hurwitz_attachment.py",
        HERE/"verify_d5_lobe_core_absorption.py",HERE/"D5_LOBE_BOUNDED_CORE_ABSORPTION.md",
        HERE/"d5_lobe_core_absorption_verification.json",
        HERE/"DEGREE_FOURTEEN_HURWITZ_ORBIT.md",
        HERE/"degree_fourteen_hurwitz_verification.json",
        HERE/"IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md",
        HERE/"sources/jones-primitive-cycle-1209.5169v1.pdf"]
 report={"schema_version":1,"result":"PASS","python":platform.python_version(),"local":loc,"old_cycle_rotation_controls":rot,
         "tree_controls":cases,
         "uniform_bound":"6656+(t-2)*446265984+sum_{j=2}^{t-1}[(2j)^2(2j-1)+14(2j)-38+2]",
         "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
         "scope":"Ordered one-vertex naturalD5 pair trees only, with arbitrary existing attachment and arbitrary shared vertex inside the new lobe. Fixed-controller endpoints are substituted, not unrecovered literal words replayed. Uniform reachability uses the written proof. No all-DT connectivity, n16 classification or arbitrary pair-attachment theorem.",
         "formal_verification":False,"independent_external_replay":False,"external_human_review":False}
 args.write_receipt.write_text(json.dumps(report,indent=2)+"\n")
 print("PASS",loc,"rotationcases",len(rot),"wholeconstructioncases",len(cases),flush=True)

if __name__=="__main__":main()
