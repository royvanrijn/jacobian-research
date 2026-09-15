#!/usr/bin/env python3
"""Replay symmetry obstruction on both retained complete determinant948 frames."""
import argparse,hashlib,json,resource,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json'
PACKET=ROOT/'artifacts/generated-results/elkies-k3-det948-all-frame-symmetry-v1'
def run():
 start=time.monotonic();data=json.loads(SOURCE.read_text());frames=data['rootless_classes'];assert len(frames)==2
 I=tuple(tuple(int(i==j) for j in range(17)) for i in range(17))
 def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(17)) for j in range(17)) for i in range(17))
 def transpose(a):return tuple(zip(*a))
 def limits():
  resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
 results=[]
 for f in frames:
  H=tuple(map(tuple,f['gram']));literal='['+';'.join(','.join(map(str,r)) for r in H)+']'
  code='default(nbthreads,1);\ndefault(parisize,134217728);\nH='+literal+';a=qfauto(H);print(a[1]);for(k=1,#a[2],print(vector(17,i,vector(17,j,a[2][k][i,j]))));print("DONE");quit;\n'
  p=subprocess.run(['gp','-q','-f'],input=code,text=True,capture_output=True,check=True,timeout=25,preexec_fn=limits)
  lines=p.stdout.strip().splitlines();assert lines[-1]=='DONE',p.stderr
  order=int(lines[0]);assert order==f['automorphism_group_order']
  gens=[tuple(map(tuple,json.loads(line))) for line in lines[1:-1]]
  assert all(mul(mul(transpose(g),H),g)==H for g in gens)
  group={I};todo=[I]
  while todo:
   x=todo.pop()
   for g in gens:
    y=mul(x,g)
    if y not in group:group.add(y);todo.append(y)
   assert len(group)<=order
  assert len(group)==order
  traces=sorted(sum(g[i][i] for i in range(17)) for g in group if g!=I and mul(g,g)==I)
  assert 1 not in traces
  results.append({'class_index':f['class_index'],'group_order':order,'generators':gens,'nonidentity_involution_traces':traces})
 assert [r['group_order'] for r in results]==[2,4]
 assert [r['nonidentity_involution_traces'] for r in results]==[[-17],[-17,-3,3]]
 return {'status':'PASS','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frames':results,'required_trace':1,'compatible_involutions':0,'classification_scope':'All rootless J2 classes of the pinned determinant948 NS, inheriting the complete Niemeier census. No J1 enumeration needed.','independent_group_order_implementation':False,'elapsed_seconds':time.monotonic()-start}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
 if args.write:
  PACKET.mkdir(exist_ok=True);(PACKET/'result.json').write_text(json.dumps(out,indent=2)+'\n')
 else:
  old=json.loads((PACKET/'result.json').read_text());assert json.loads(json.dumps({k:v for k,v in out.items() if k!='elapsed_seconds'}))=={k:v for k,v in old.items() if k!='elapsed_seconds'}
 print('PASS: both determinant948 rootless frame groups closed; involution traces -17 and {-17,-3,3}; none is1')
