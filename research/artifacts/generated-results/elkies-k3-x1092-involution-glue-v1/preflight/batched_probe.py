import json,subprocess,resource,time
from pathlib import Path
x=json.load(open('research/artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json'))['rootless_classes']
lines=['default(parisizemax,512000000);']
for c in x:
 if c['automorphism_group_order']==2:continue
 gram='['+';'.join(','.join(map(str,row)) for row in c['gram'])+']'
 lines.append('G=qfauto('+gram+');print(['+str(c['class_index'])+',G[1],apply(M->vector(17,i,vector(17,j,M[i,j])),G[2])]);')
code='\n'.join(lines)+'\n';Path('/tmp/x1092_aut_probe.gp').write_text(code)
def limits():resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
r=subprocess.run(['gp','-q','-f'],input=code,text=True,capture_output=True,timeout=25,preexec_fn=limits);assert r.returncode==0,r.stderr
rows=[json.loads(line) for line in r.stdout.splitlines()];Path('/tmp/x1092_aut_probe.json').write_text(json.dumps(rows,indent=2)+'\n');print([(i,o,len(g)) for i,o,g in rows])
