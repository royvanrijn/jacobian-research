import json,subprocess,time
from pathlib import Path
x=json.load(open('research/artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json'))['rootless_classes'];out=Path('/tmp/x1092_aut_checkpoints');out.mkdir(exist_ok=True)
for c in x:
 i=c['class_index']
 if c['automorphism_group_order']==2 or i==2:continue
 gram='['+';'.join(','.join(map(str,r)) for r in c['gram'])+']';code='G=qfauto('+gram+');print([G[1],apply(M->vector(17,i,vector(17,j,M[i,j])),G[2])]);\n';started=time.monotonic()
 try:
  p=subprocess.run(['gp','-q','-f'],input=code,text=True,capture_output=True,timeout=5)
  result=json.loads(p.stdout);assert p.returncode==0 and result[0]==c['automorphism_group_order']
  (out/f'class-{i}.json').write_text(json.dumps({'class_index':i,'elapsed':time.monotonic()-started,'group':result,'input_gp':code},indent=2)+'\n');print(i,'PASS',time.monotonic()-started,flush=True)
 except subprocess.TimeoutExpired:
  (out/f'class-{i}-timeout.json').write_text(json.dumps({'class_index':i,'status':'TIMEOUT','wall_seconds':5,'input_gp':code})+'\n');print(i,'TIMEOUT',flush=True);break
