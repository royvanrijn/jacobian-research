#!/usr/bin/env python3
"""Exact port certificates for the uniform threaded-curl transfer.

Two finite tangle cells prove the unbounded induction in the companion note.
Literal source/target braid wiring controls supplement that proof for 5<=b<=64.
No hyperbolic recognition is used here; the fixed parent's separate interval
certificate is a pinned dependency.
"""
from __future__ import annotations
import argparse
from collections import deque
from copy import deepcopy
from pathlib import Path
import hashlib,json,platform

HERE=Path(__file__).resolve().parent
INPUT=HERE/'threaded_curl_uniform_transfer_inputs.json'


def require(value,message):
 if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def p(node,port):return(node,port%4)


def join(g,a,b):
 require(a!=b and a not in g and b not in g,'fresh paired ports')
 g[a]=b;g[b]=a


def graph(edges):
 g={}
 for a,b in edges:join(g,tuple(a),tuple(b))
 return g


def edges(g):
 return {frozenset((a,b))for a,b in g.items()}


def riii(g,A,a,B,b,C,c):
 require(len({A,B,C})==3,'three distinct RIII crossings')
 require(sum(j%2 for j in(a,b,c))in(1,2),'RIII over/under admissibility')
 tri=[(p(A,a+1),p(B,b)),(p(B,b+1),p(C,c)),(p(C,c+1),p(A,a))]
 require(all(g.get(x)==y for x,y in tri),'RIII triangular-face boundary')
 old=[p(C,c-1),p(C,c-2),p(A,a-1),p(A,a-2),p(B,b-1),p(B,b-2)]
 new=[p(A,a),p(B,b+1),p(B,b),p(C,c+1),p(C,c),p(A,a+1)]
 sub=dict(zip(old,new));oldset=set(old)
 border={frozenset((x,g[x]))for x in old}
 for e in border:
  for x in e:g.pop(x,None)
 for x,y in tri:g.pop(x,None);g.pop(y,None)
 for e in border:
  x,y=tuple(e);join(g,sub.get(x,x),sub.get(y,y))
 for x,y in[(p(A,a-1),p(B,b+2)),(p(B,b-1),p(C,c+2)),(p(C,c-1),p(A,a+2))]:join(g,x,y)


def ri(g,A,corner):
 x,y=p(A,corner),p(A,corner+1)
 require(g.get(x)==y,'RI monogon edge')
 x2,y2=p(A,corner+2),p(A,corner+3);u,v=g[x2],g[y2]
 require(u[0]!=A and v[0]!=A,'RI outside strands')
 for q in(x,y,x2,y2,u,v):g.pop(q,None)
 join(g,u,v)


def rii_at(g,A,B):
 for a in range(4):
  q=g.get(p(A,a));z=g.get(p(A,a+1))
  if q is None or z is None:continue
  node,b=q
  if node==B and z==p(B,b-1)and(a+b)%2==0:
   out=[g[p(A,a+2)],g[p(A,a+3)],g[p(B,b+1)],g[p(B,b+2)]]
   require(all(x[0]not in(A,B)for x in out),'RII four outside strands')
   for i in range(4):g.pop(p(A,i));g.pop(p(B,i))
   for q in out:g.pop(q)
   join(g,out[0],out[3]);join(g,out[1],out[2])
   return True
 return False


def literal_braid(word,n):
 g={};top=[None]*n;wire=[None]*n
 for C,letter in enumerate(word):
  j=abs(letter)-1
  incoming,outgoing=([0,3],[1,2])if letter>0 else([1,0],[2,3])
  for col,port in[(j,incoming[0]),(j+1,incoming[1])]:
   q=p(C,port)
   if wire[col]is None:top[col]=q
   else:join(g,wire[col],q)
  wire[j:j+2]=[p(C,i)for i in outgoing]
 for a,b in zip(wire,top):join(g,a,b)
 require(len(g)==4*len(word),'all literal braid ports paired')
 return g


def compare_fixed_U(A,B,source_U,target_U):
 """Connected oriented rotation-system isomorphism, all U ports fixed."""
 mapping={source_U:(target_U,0)};used={target_U};todo=deque([source_U])
 while todo:
  X=todo.popleft();Y,s=mapping[X]
  for j in range(4):
   C,k=A[p(X,j)];D,l=B[p(Y,j+s)];t=(l-k)%4
   require(t%2==0,'planar isomorphism preserves over/under branches')
   if C in mapping:require(mapping[C]==(D,t),'consistent rotation-system map')
   else:
    require(D not in used,'injective rotation-system map')
    mapping[C]=(D,t);used.add(D);todo.append(C)
 require(len(mapping)*4==len(A)==len(B),'onto rotation-system map')
 for j in range(14):require(mapping[source_U+j]==(target_U+j,0),'each outside U port fixed')
 return mapping


def disk_contract(g,cell):
 nodes={a[0]for a in g};boundary={a[0]for a in g if str(a[0]).startswith('e')}
 reached={next(iter(nodes))};todo=list(reached)
 for node in todo:
  for a,b in g.items():
   if a[0]==node and b[0]not in reached:reached.add(b[0]);todo.append(b[0])
 require(reached==nodes,'connected embedded cell graph')
 require(all(sum(a[0]==node for a in g)==(1 if node in boundary else 4)for node in nodes),
         'crossing valence four and boundary valence one')
 unused=set(g);faces=[]
 while unused:
  start=x=next(iter(unused));face=[]
  while True:
   require(x in unused,'simple rotation-system face partition');unused.remove(x);face.append(x)
   x=g[(x[0],0 if x[0]in boundary else(x[1]+1)%4)]
   if x==start:break
  faces.append(face)
 require(len(nodes)-len(g)//2+len(faces)==2,'cell rotation system is spherical')
 outer=[f for f in faces if any(x[0]in boundary for x in f)]
 require(len(outer)==1 and {x[0]for x in outer[0]if x[0]in boundary}==boundary,
         'all marked boundary ports lie on one outer face')
 order=[x[0]for x in outer[0]if x[0]in boundary];at=order.index('e0');order=order[at:]+order[:at]
 require(order==cell['boundary_cyclic_order'],'same prescribed cyclic boundary order')
 inputs=set(cell['boundary_inputs'])
 def incoming(x):return x[0]not in inputs if x[0]in boundary else x[1]in(0,3)
 require(all(incoming(a)!=incoming(b)for a,b in g.items()),'oriented boundary strands and positive crossings')

def local_cells(data):
 for name,cell in data['local_cells'].items():
  g=graph(cell['before']);disk_contract(g,cell)
  for move in cell['moves']:
   if move[0]=='RIII':riii(g,*move[1:])
   elif move[0]=='RI':ri(g,*move[1:])
   else:raise AssertionError('unknown local move')
  require(edges(g)==edges(graph(cell['after'])),'finite boundary-port cell '+name)
  disk_contract(g,cell)
  # Every named boundary port is still exactly the same outside endpoint.
  require({x for x in g if str(x[0]).startswith('e')}==
          {tuple(x)for x in cell['boundary_ports']},'same external boundary ports')
 return {'clasp_slide_two_RIII':True,'cap_one_RIII_and_one_RI':True,
         'boundary_ports_fixed':True,'over_under_branches_checked':True,
         'embedded_disk_boundaries_and_orientations':True,'cell_count':len(data['local_cells'])}


def source_template(d,Ustart,Y,u,v):
 """Closed-form outside-U graph; boundary ports are actual U ports."""
 es=[]
 def e(a,b):
  pair=frozenset((a,b))
  if pair not in {frozenset(x)for x in es}:es.append((a,b))
 e(p(Y[0],0),p(Ustart+12,1));e(p(Y[0],3),p(Ustart+13,2))
 for i in range(d-1):
  if i:e(p(Y[i],0),p(Y[i-1],2))
  e(p(Y[i],1),p(u[d-1-i],0))
  if i<d-2:e(p(Y[i],2),p(Y[i+1],0))
  else:e(p(Y[i],2),p(u[0],0))
  if i:e(p(Y[i],3),p(v[d-2-i],1))
 for j in range(d):
  if j<d-1:e(p(u[j],1),p(u[j+1],3))
  else:e(p(u[j],1),p(v[j],0))
  e(p(u[j],2),p(v[j],3))
  if j:e(p(u[j],3),p(u[j-1],1))
  else:e(p(u[j],3),p(v[0],2))
  if j<d-1:e(p(v[j],0),p(v[j+1],2))
  if j<=d-3:e(p(v[j],1),p(Y[d-2-j],3))
  elif j==d-2:e(p(v[j],1),p(Ustart+1,3))
  else:e(p(v[j],1),p(Ustart,3))
 return {frozenset(x)for x in es}


def target_cut_certificate(data):
 contract=data['target_cut_contract'];R=data['fixed_six_braid_R'];g=literal_braid(R,6)
 require(contract['literal_cut_edges']==[[[24,2],[0,0]],[[25,1],[0,3]]],
         'actual fourth and fifth closure arcs')
 for a,b in contract['literal_cut_edges']:
  a,b=tuple(a),tuple(b);require(g[a]==b,'literal cut edge pairing');g.pop(a);g.pop(b)
 for name,endpoint in contract['before_boundary_attachments'].items():join(g,(name,0),tuple(endpoint))
 require(contract['target_RII_pairs']==[[6,25],[7,22]],'fixed target RII indices')
 for A,B in contract['target_RII_pairs']:
  require(A!=0 and B!=0 and not(8<=A<22 or 8<=B<22),'target cut and U neighborhoods untouched')
  require(rii_at(g,A,B),'target RII valid with both cut arcs open')
 for name,endpoint in contract['after_boundary_attachments'].items():
  require(g[(name,0)]==tuple(endpoint),'fixed named open target boundary '+name)
 require(len({x[0]for x in g if isinstance(x[0],int)})==22,'normalized target base has22 crossings')
 # Reclose the named pairs, and bind the reduced cut edges to the source ladder.
 for x,y in [('Kin','Kout'),('Cin','Cout')]:
  a,b=g[(x,0)],g[(y,0)]
  for t in[(x,0),(y,0),a,b]:g.pop(t)
  join(g,a,b)
 U=data['prefix_U'];Q=[4]+U+[-4,2,3];P=[5,4,3,2,2,3,4,-5]
 src=literal_braid(Q+P,6);require(rii_at(src,18,25)and rii_at(src,0,15),'source b5 normalization')
 mapping=compare_fixed_U(src,g,1,8)
 require(mapping[17]==(24,0)and mapping[19]==(0,0)and mapping[24]==(5,0),
         'Y1,u0,v0 target cut identification')
 # The cap output is exactly two ordinary positive braid crossings.
 cell=graph(data['local_cells']['cap']['after'])
 require(cell[('u',1)]==('v',0)and cell[('u',2)]==('v',3),
         'positive two-strand braid full-twist cell')
 require(cell[('e0',0)]==('u',0)and cell[('e3',0)]==('u',3)
         and cell[('e1',0)]==('v',1)and cell[('e2',0)]==('v',2),
         'four boundary ports of repeated positive full twist')
 return {'fixed_target_RII_pairs':[[6,25],[7,22]],'both_cut_arcs_kept_open':True,
         'named_boundary_ports_fixed':True,'normalized_base_crossings':22,
         'source_cut_edges':['Y1.2-u0.0','v0.2-u0.3'],
         'inserted_clasp_is_positive_sigma4_squared':True}

def exact_controls(data):
 U=data['prefix_U'];R=data['fixed_six_braid_R']
 require(U==[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3],'literal U')
 require(R==[4,3,2,2,3,4,-5,4]+U+[-4,2,3,5],'literal fixed R')
 checked=[]
 for b in range(5,65):
  d=b-2;r=b-4;q=b-5;rho=list(range(b-1,3,-1))
  Q=rho+U+[-j for j in rho[::-1]]+list(range(2,b-1))
  P=list(range(b,1,-1))+[2]+list(range(3,b))+[-b]
  g=literal_braid(Q+P,b+1)
  Y=[2*b+6+j for j in range(d-1)];u=[3*b+4+j for j in range(d)];v=[5*b-1-j for j in range(d)]
  A,B=len(Q),len(Q+P)-1
  require(g[p(A,2)]==p(B,0)and g[p(A,3)]==p(B,3),'universal outer RII ports')
  require(rii_at(g,A,B),'outer-circle detour RII')
  require(g[p(0,3)]==p(2*b+5,3)and g[p(u[0],3)]==p(v[0],2),'outer RII seeds rho bigon and circle return')
  for i in range(r):
   A,B=i,2*r+13-i
   require(g[p(A,2)]==p(B,0)and g[p(A,3)]==p(B,3),'universal rho RII ports')
   require(g[p(A,0)]==p(v[i],1)and g[p(A,1)]==
           (p(i+1,3)if i<r-1 else p(r+2,0)), 'rho positive outside ports')
   require(g[p(B,1)]==p(B-1,3)and g[p(B,2)]==p(Y[r-i],3),'rho negative outside ports')
   require(rii_at(g,A,B),'rho conjugation detour RII')
  outside={edge for edge in edges(g)if any(not(r<=x[0]<r+14)for x in edge)}
  require(outside==source_template(d,r,Y,u,v),'literal source matches closed-form ladder ports')
  RIIIcount=0
  for s in range(d-3):
   t=Y[d-2-s]
   for j in range(s+1):
    riii(g,u[j],3,t,2,v[s],1);RIIIcount+=1
    if j<s:riii(g,v[j],3,v[s],2,t,1);RIIIcount+=1
   ri(g,t,2)
  require(RIIIcount==q*q,'square RIII count')
  target=literal_braid(R+[4]*(2*q),6)
  # Remove target detours; no U crossing may be removed.
  while True:
   changed=False;nodes=sorted({a[0]for a in target})
   for A in nodes:
    if 8<=A<22:continue
    partners={target[p(A,j)][0]for j in range(4)}
    for B in sorted(partners):
     if B==A or 8<=B<22:continue
     if rii_at(target,A,B):changed=True;break
    if changed:break
   if not changed:break
  compare_fixed_U(g,target,r,8)
  checked.append({'b':b,'RIII':RIIIcount,'RI':q,'rho_and_circle_RII':r+1})
 candidate=HERE/'protected_universal_parent_candidate_inputs.json'
 N=json.loads(candidate.read_text())
 require(N['base_braid_R']==R and N['new_circle_suffix_D']==[6,5,4,4,5,-6],
         'fixed hyperbolic parent input binding')
 return {'bounded_literal_braid_and_port_controls':checked,
  'fixed_U_ports_in_final_spherical_isotopy':14*4,
  'uniform_transfer_proved_by_written_port_induction':
   'M_b unmarked homeomorphic to N(D:(1,-(b-5))) for every b>=5',
  'unbounded_primality_uses_fixed_parent_geometry_and_two_successive_Dehn_filling_arguments':True,
  'fixed_parent_input_sha256':sha(candidate)}


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-receipt',type=Path)
 args=parser.parse_args();data=json.loads(INPUT.read_text())
 report={'schema_version':1,'python':platform.python_version(),'checker_sha256':sha(__file__),
  'input_sha256':sha(INPUT),'local_tangle_cells':local_cells(data),
  'target_cut_contract':target_cut_certificate(data),'literal_controls':exact_controls(data)}
 # A corrupt over/under port and a switched boundary endpoint must fail locally.
 bad=deepcopy(data);bad['local_cells']['cap']['moves'][0][2]=2
 try:local_cells(bad)
 except AssertionError:pass
 else:raise AssertionError('corrupt RIII face accepted')
 bad=deepcopy(data);bad['local_cells']['cap']['after'][0][1][0]='e1'
 try:local_cells(bad)
 except AssertionError:pass
 else:raise AssertionError('switched external boundary port accepted')
 report['negative_controls_rejected']=2
 report['assurance']={'independent_replay':False,'formal_verification':False,'external_review':False}
 report['result']='PASS'
 if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: two finite port cells; literal controls b=5..64; outside U ports fixed; two corrupt controls rejected.')
 print('checker_sha256='+report['checker_sha256']);print('input_sha256='+report['input_sha256'])

if __name__=='__main__':main()