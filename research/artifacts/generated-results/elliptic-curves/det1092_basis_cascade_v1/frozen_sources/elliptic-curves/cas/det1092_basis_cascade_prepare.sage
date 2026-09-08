#!/usr/bin/env sage-python
"""Freeze four generic-only basis arms before any reconstruction outcomes."""
import csv,hashlib,json,os,random,subprocess,sys,datetime
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,identity_matrix,matrix
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1';PKG=ART/'det1092_basis_cascade_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip()
def write(p,d):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump({'starting_commit':commit,**d},f,sort_keys=True,indent=2);f.write('\n')
parent_path=ART/'curve302_recovered_mw17_parent_v1.json';parent=read(parent_path)
R=PolynomialRing(QQ,'t');K=R.fraction_field()
def rf(v):return K(R(v['numerator']))/R(v['denominator'])
aa=list(map(rf,parent['a_invariants']));E=EllipticCurve(K,aa)
P=[E(rf(x),rf(y)) for x,y in parent['basis_weierstrass_coordinates']]
aa0=[v(0) for v in aa];a1,a2,a3,a4,a6=aa0;b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
model=[0,0,0,-27*(b2*b2-24*b4),-54*(-b2**3+36*b2*b4-216*b6)]
short=EllipticCurve(QQ,model)
def specialize(q):
 x,y=q[0](0),q[1](0)
 return short(36*x+3*b2,108*(2*y+a1*x+a3))
reference=list(map(specialize,P))
transforms={'original':identity_matrix(ZZ,17)}
u=identity_matrix(ZZ,17);u[13,15]=-1;transforms['recovered-P14']=u
u=identity_matrix(ZZ,17);u[10]=[0,-1,-1,1,0,0,1,0,-1,0,1,-1,0,0,0,0,1];transforms['recovered-P11']=u
u=identity_matrix(ZZ,17);rng=random.Random(10921731);operations=[]
for k in range(12):
 i,j=rng.sample(range(17),2);sign=rng.choice([-1,1]);u[i]+=sign*u[j];operations.append([i,j,sign])
transforms['random-unimodular']=u
assert not D.exists() and not PKG.exists()
D.mkdir(parents=True);PKG.mkdir()
os.environ['CASCADE_ARM']='original'
worker=SourceFileLoader('basis_worker',str(CAS/'det1092_basis_cascade_v1.sage')).load_module()
import audit_recorded_point_mod2_rank_v3 as mod2
import audit_retained_cloud_modl as modl
sources={**worker.sources(),**mod2.sources(),**modl.sources(),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}
# Preserve all CAS Python/Sage dependencies, even transitive modules not in old source manifests.
for p in CAS.rglob('*'):
 if p.is_file() and p.suffix in ('.py','.sage','.cpp','.h'):
  sources[str(p.relative_to(ROOT))]=sha(p)
for name,digest in sources.items():
 dest=PKG/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((ROOT/name).read_bytes())
worker_sources=worker.sources()
basepolicy={'schema':'det1092-basis-cascade.v1','starting_commit':commit,'anchors_per_shell':16,'canonical_per_shell':25,'exact_cvp_node_limit':2000000,'height':125000,'seconds_per_chart':10,'max_charts':4096,'max_epochs':20,'target_rank':31,'prime_bound':1000,'gp_sha256':sha(Path('/usr/bin/gp')),'sources':worker_sources,
 'metric':'Round canonical heights at10^6 in original generic17 plus own discovered points, then transport exactly by blockdiag(U,I). Exact integer CVP; same underlying metric for equal point extensions.',
 'selection':'Unchanged V3 cheap_shortlist and final_shortlist; full binary extension enumeration; actual-point dedup; stop first certified gain and rebuild; complete no-gain epoch terminates.',
 'resource_limits':{'own_CPU_seconds':7200,'hard_CPU_seconds':7210,'wall_seconds':14400,'address_space_GiB':12,'concurrent_arms':4},
 'allowed':['own transformed generic basis and exact inverse/original generic reference','own transported generic orbit catalogue','own checkpoints and discoveries','pinned CAS algorithms and PARI binary'],
 'forbidden':['other arm outputs','historical exceptional points or saved successful cascade charts','public rank31 points and exceptional complements','network queries or outcome-fitted selectors'],
 'boundary':'Designer and V3 algorithm are retrospectively calibrated. Fresh execution has no exceptional points. This is specialized visibility of the same integral MW17 group, not generic rank raising or a new elliptic curve.',
 'failure_semantics':'Preserve all failures and timeouts. No silent reruns. Wall/CPU censoring does not imply no gain or nonexistence. No different parameters are selected in this experiment.'}
original_orbits=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
roster=[];bindings={str(parent_path.relative_to(ROOT)):sha(parent_path),str(original_orbits.relative_to(ROOT)):sha(original_orbits), 'MATH_STATUS.json':sha(ROOT/'MATH_STATUS.json')}
for aid,U in transforms.items():
 assert abs(U.det())==1;inv=U.inverse();assert inv.denominator()==1
 points=[sum((int(c)*p for c,p in zip(row,reference)),short(0)) for row in U.rows()]
 if aid.startswith('recovered'):
  omitted=14 if aid.endswith('14') else 11;candidate_index=0 if omitted==14 else 2
  source=ART/f'det1092_blind_mw16_v1/arms/arm-{omitted:02d}/certificate.json'
  bindings[str(source.relative_to(ROOT))]=sha(source)
  cr=read(source)['candidates'][candidate_index]['coordinates'];q=E(*map(rf,cr));expected=sum((int(c)*p for c,p in zip(U[omitted-1],P)),E(0))
  assert q==expected and specialize(q)==points[omitted-1]
 arm=D/aid;arm.mkdir()
 fixture={'arm_id':aid,'curve':list(map(str,model)),'points':[list(map(str,q.xy())) for q in points],'reference_points':[list(map(str,q.xy())) for q in reference], 'basis_transform':[list(map(int,row)) for row in U.rows()], 'inverse_transform':[list(map(int,row)) for row in inv.rows()], 'determinant':int(U.det()),'generic_gram':[list(map(str,row)) for row in (U*matrix(QQ,parent['generic_height_gram'])*U.transpose()).rows()], 'random_seed':10921731 if aid=='random-unimodular' else None, 'random_elementary_operations':operations if aid=='random-unimodular' else []}
 write(arm/'fixture.json',fixture);write(PKG/'fixtures'/f'{aid}.json',fixture)
 with original_orbits.open() as src,(arm/'orbits.tsv').open('x') as out:
  rd=csv.DictReader(src,delimiter='\t');wr=csv.DictWriter(out,fieldnames=rd.fieldnames,delimiter='\t',lineterminator='\n');wr.writeheader()
  for row in rd:
   word=matrix(ZZ,1,17,list(map(int,row['parent_MW17_w'].split())))*inv
   row['parent_MW17_w']=' '.join(map(str,word.row(0)));wr.writerow(row)
 policy={**basepolicy,'arm_id':aid,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [arm/'fixture.json',arm/'orbits.tsv']}}
 write(arm/'protocol.json',policy)
 roster.append({'arm_id':aid,'fixture_sha256':sha(arm/'fixture.json'),'orbits_sha256':sha(arm/'orbits.tsv'),'protocol_sha256':sha(arm/'protocol.json'),'transform_determinant':int(U.det())})
write(PKG/'roster.json',{'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'arms':roster})
write(PKG/'protocol.json',{'roster_sha256':sha(PKG/'roster.json'),'policy':basepolicy,'all_source_bindings':sources,'preparation_bindings':bindings,'novelty_policy':'Record every tested pointed quartic with exact birational mapping to302. Classify new equation presentations against historical charts only after all arms terminal. A new elliptic-curve class requires exact Q-isomorphism exclusion against repository catalogue; no global novelty claim from absence.'})
for p in list(PKG.rglob('*'))+list(D.rglob('*')):
 if p.is_file():p.chmod(0o444)
print('FROZEN FOUR ARMS',commit,sha(PKG/'roster.json'),flush=True)
