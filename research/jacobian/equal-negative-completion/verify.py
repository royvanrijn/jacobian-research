#!/usr/bin/env python3
"""Independent AST/Fraction coefficient reconstruction and certificate multiplication.
No CAS or primary solver algebra is imported.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import ast,hashlib,json,time,re
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,q,c=Q(1)):
 r=dict(p)
 for e,v in q.items():r[e]=r.get(e,Q(0))+c*v
 return {e:v for e,v in r.items() if v}
def mul(p,q):
 r={}
 for e,a in p.items():
  for f,b in q.items():k=tuple(x+y for x,y in zip(e,f));r[k]=r.get(k,Q(0))+a*b
 return {e:v for e,v in r.items() if v}
def const(v,n):return {(0,)*n:Q(v)} if v else {}
def var(i,n):e=[0]*n;e[i]=1;return {tuple(e):Q(1)}
def power(p,k,n):
 assert isinstance(k,int) and k>=0;r=const(1,n)
 for _ in range(k):r=mul(r,p)
 return r
def parse(text,names):
 n=len(names)
 def visit(t):
  if isinstance(t,ast.Constant):assert isinstance(t.value,int);return const(t.value,n)
  if isinstance(t,ast.Name):return var(names.index(t.id),n)
  if isinstance(t,ast.UnaryOp) and isinstance(t.op,ast.USub):return {e:-v for e,v in visit(t.operand).items()}
  assert isinstance(t,ast.BinOp)
  if isinstance(t.op,ast.Pow):return power(visit(t.left),ast.literal_eval(t.right),n)
  a,b=visit(t.left),visit(t.right)
  if isinstance(t.op,ast.Add):return add(a,b)
  if isinstance(t.op,ast.Sub):return add(a,b,Q(-1))
  if isinstance(t.op,ast.Mult):return mul(a,b)
  assert isinstance(t.op,ast.Div) and len(b)==1 and (0,)*n in b;return {e:v/b[(0,)*n] for e,v in a.items()}
 return visit(ast.parse(text.replace('^','**'),mode='eval').body)
def derivative(p,i):
 out={}
 for e,v in p.items():
  if e[i]:f=list(e);f[i]-=1;out[tuple(f)]=v*e[i]
 return out
def buckets(p):
 out={}
 for e,v in p.items():key=e[:2];tail=(0,0)+e[2:];out.setdefault(key,{})[tail]=v
 return [out[key] for key in sorted(out,reverse=True)]
def condition(A,B,n,k):
 jac=add(mul(derivative(A,0),derivative(B,1)),mul(derivative(A,1),derivative(B,0)),Q(-1));out=mul(add(const(1,n),var(0,n)),jac);out=add(out,mul(A,derivative(B,1)));out=add(out,mul(B,derivative(A,1)),Q(-1));return add(out,const(-k,n))
def parse_flat(text,names):
 n=len(names);positions={v:i for i,v in enumerate(names)};out={}
 for token in re.findall(r'[+-]?[^+-]+',text.replace(' ','')):
  sign=-1 if token.startswith('-') else 1;token=token.lstrip('+-');es=[0]*n;c=Q(sign)
  for factor in token.split('*'):
   if '^' in factor:base,k=factor.split('^');k=int(k);assert k>=0
   else:base=factor;k=1
   if base in positions:es[positions[base]]+=k
   else:assert k==1;c*=Q(base)
  e=tuple(es);out[e]=out.get(e,Q(0))+c
 return {e:v for e,v in out.items() if v}
def general_condition(f,A,B,n,k):
 jac=lambda a,b:add(mul(derivative(a,0),derivative(b,1)),mul(derivative(a,1),derivative(b,0)),Q(-1));out=mul(f,jac(A,B));out=add(out,mul(A,jac(f,B)));out=add(out,mul(B,jac(A,f)));return add(out,const(-k,n))
def certificate(folder,names,equations):
 n=len(names);log=folder/'run.log';text=log.read_text();assert 'UNIT_IDEAL' in text and text.count('CERT_BEGIN')==text.count('CERT_END')==1;lines=text.split('CERT_BEGIN')[1].split('CERT_END')[0].strip().splitlines();assert len(lines)==len(equations);total={};terms=0;maximum=0
 for i,(line,eq) in enumerate(zip(lines,equations),1):
  index,pol=line.split('|',1);assert int(index)==i;mult=parse_flat(pol,names);assert all(e[0]==e[1]==0 for e in mult);terms+=len(mult);maximum=max(maximum,max((sum(e) for e in mult),default=0));total=add(total,mul(mult,eq))
 assert total==const(1,n);cp=json.loads((folder/'checkpoint.json').read_text());assert cp['status']=='terminal' and cp['exit_code']==0 and cp['script_sha256']==sha(folder/'model.sing') and cp['log_sha256']==sha(log);print('CERTIFICATE PASS',folder.name,len(equations),'generators',terms,'multiplier terms,maxdegree',maximum,flush=True);return dict(chart=folder.name,input_sha256=sha(folder/'input.json'),script_sha256=sha(folder/'model.sing'),certificate_log_sha256=sha(log),generators=len(equations),multiplier_terms=terms,maximum_multiplier_degree=maximum,exact_rational_identity_one=True)
def verify_chart(family,tag):
 folder=HERE/(family+'-'+tag);d=json.loads((folder/'input.json').read_text());names=['Y','Z']+d['variables'];n=len(names);A=parse(('Y' if tag=='I' else 'Z')+'+a0*Y**2+a1*Y*Z+a2*Z**2+a3*Y**3+a4*Y**2*Z+a5*Y*Z**2+a6*Z**3+a7*Y**4',names);linearB=('Z' if tag=='I' else ('Y' if family=='normal' else 'Y+t*Z'));B=parse(linearB+'+b0*Y**2+b1*Y*Z+b2*Z**2+b3*Y**3',names);f=parse({'normal':'1+Y','parabola':'1+Z+Y**2','cstar':'1+Y*Z+Y**3'}[family],names);assert A==parse(d['A'],names) and B==parse(d['B'],names);p=general_condition(f,A,B,n,1 if tag=='I' else -1);assert p==parse(d['condition'],names);equations=buckets(p)+[parse('rho*a7-1',names)];assert equations==[parse(s,names) for s in d['equations']];return certificate(folder,names,equations)
def verify_opposite():
 folder=HERE/'opposite-cubic';d=json.loads((folder/'input.json').read_text());names=['U','V']+d['variables'];n=len(names);f=parse('1+a0*U+a1*V+a2*U**2+a3*U*V+a4*V**2+a5*U**3',names);B=parse('V+b0*U**2+b1*U*V+b2*V**2+b3*U**3',names);g=parse('U+U**2',names);assert f==parse(d['f'],names) and B==parse(d['B'],names) and g==parse(d['g'],names);jac=lambda a,b:add(mul(derivative(a,0),derivative(b,1)),mul(derivative(a,1),derivative(b,0)),Q(-1));p=add(add(mul(f,jac(g,B)),mul(g,jac(f,B))),mul(B,jac(g,f)));p=add(p,const(-1,n));assert p==parse(d['condition'],names);eq=buckets(p)+[parse('rho*a5-1',names)];assert eq==[parse(t,names) for t in d['equations']]
 def lift(poly,shift):
  out={}
  for e,c in poly.items():assert e[0]+e[1]+shift>=0;out[(e[0],e[0]+e[1]+shift,e[1])+e[2:]]=c
  return out
 nn=n+1;xx,yy=var(0,nn),var(1,nn);F=[add(xx,mul(power(xx,2,nn),yy)),lift(f,1),lift(B,-1)];full=det([[derivative(v,j) for j in range(3)] for v in F],nn);assert full==lift(add(p,const(1,n)),0);result=certificate(folder,names,eq);result['full_opposite_3D_jacobian_checked']=True;return result
def det(M,n):
 out={}
 for p in permutations(range(3)):
  value=const((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)),n)
  for i in range(3):value=mul(value,M[i][p[i]])
  out=add(out,value)
 return out
def generic_dictionary():
 names=['Y','Z']+['c'+str(i) for i in range(9)]+['a'+str(i) for i in range(12)]+['b'+str(i) for i in range(12)];n=len(names);f=const(1,n);i=0
 for degree in range(1,4):
  for j in range(degree,-1,-1):mon=mul(power(var(0,n),j,n),power(var(1,n),degree-j,n));f=add(f,mul(var(2+i,n),mon));i+=1
 A=var(0,n);B=var(1,n);i=0
 for degree in range(2,5):
  for j in range(degree,-1,-1):mon=mul(power(var(0,n),j,n),power(var(1,n),degree-j,n));A=add(A,mul(var(11+i,n),mon));B=add(B,mul(var(23+i,n),mon));i+=1
 def lift(p,shift):
  out={}
  for e,c in p.items():assert e[0]+e[1]+shift>=0;out[(e[0]+e[1]+shift,e[0],e[1])+e[2:]]=c
  return out
 F=[lift(f,1),lift(A,-1),lift(B,-1)];source_n=n+1;jac=det([[derivative(p,j) for j in range(3)] for p in F],source_n);expected=lift(add(general_condition(f,A,B,n,1),const(1,n)),0);assert jac==expected;print('GENERIC33-PARAMETER FULL3D JACOBIAN DICTIONARY PASS',flush=True);return dict(parameters=33,all_first_degrees_through7=True,full_source_jacobian_checked=True)
def source_shear_control():
 names=['a','b','c'];n=3;aa,bb,cc=[var(i,n) for i in range(n)];texts=['a-3*a**2*b-a**3*c','3*(1+2*a*b)**2*a*c+(9*(1+2*a*b)**2-6*(1+2*a*b)-2)*b','(1+2*a*b)**3*c+(6*(1+2*a*b)**2+2*(1+2*a*b))*b**2'];F=[parse(t,names) for t in texts];assert det([[derivative(f,j) for j in range(3)] for f in F],n)==const(1,n);E=[aa,add(bb,mul(aa,cc),Q(-1,3)),cc]
 def sub(poly):
  out={}
  for e,c in poly.items():
   v=const(c,n)
   for j,k in enumerate(e):v=mul(v,power(E[j],k,n))
   out=add(out,v)
  return out
 G=[sub(f) for f in F];assert G[0]==parse('a-3*a**2*b',names);assert det([[derivative(f,j) for j in range(3)] for f in G],n)==const(1,n);degree=[max(map(sum,f)) for f in G];assert degree==[3,7,9];counts=[sum(sum(e)>=2 for e in f) for f in G];assert counts==[1,7,11]
 def evaluate(f,p):return sum(c*Q(p[0])**e[0]*Q(p[1])**e[1]*Q(p[2])**e[2] for e,c in f.items())
 points=[(Q(0),Q(1),Q(-8)),(Q(-1),Q(-1,3),Q(4))];images=[[evaluate(f,p) for f in G] for p in points];assert points[0]!=points[1] and images[0]==images[1]==[Q(0),Q(1),Q(0)];out=dict(status='pass',helper_sha256=sha(Path(__file__)),source_shear='b -> b-a*c/3',source_weights=[1,-1,-2],component_degrees=degree,nonlinear_occurrences=counts,total_nonlinear_occurrences=19,dimension=3,components=[[[list(e),str(c)] for e,c in sorted(f.items())] for f in G],collision_points=[[str(v) for v in p] for p in points],common_image=[str(v) for v in images[0]],determinant_one_checked=True,first_component_cubic=True,global_construction_improvement=False,scope='Known map under an elementary source shear; degree9/19 nonlinear occurrences, not a global simplicity improvement. Provides a checked positive control with cubic component; GLDC8 gives only the lower bound6 in that extra class.');(HERE/'source-shear-validation.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');return dict(validation_sha256=sha(HERE/'source-shear-validation.json'),ordinary_degree=9,nonlinear_occurrences=19,cubic_component=True,global_improvement=False)
def main():
 start=time.monotonic();root=HERE.parents[1];ledger={e['id']:e for e in json.loads((root/'MATH_STATUS.json').read_text())['entries']};premises={}
 for key in ['GLDC5','GLDC7','GLDC8','KPAIR1']:
  e=ledger[key];assert e['state']=='proved' and e['artifact_hash']=='sha256:'+sha(root/e['checker']);premises[key]=dict(checker_sha256=sha(root/e['checker']),proof_sha256=sha(root/e['canonical_source']))
 original=json.loads((HERE/'input.json').read_text());names=['Y','Z']+original['variables'];n=len(names);A=var(0,n);B=var(1,n);index=0
 for degree in range(2,5):
  for j in range(degree,-1,-1):mon=mul(power(var(0,n),j,n),power(var(1,n),degree-j,n));A=add(A,mul(var(2+index,n),mon));B=add(B,mul(var(14+index,n),mon));index+=1
 assert index==12 and A==parse(original['A'],names) and B==parse(original['B'],names);p=condition(A,B,n,1);assert p==parse(original['condition'],names) and buckets(p)==[parse(s,names) for s in original['equations']];cp=json.loads((HERE/'checkpoint.json').read_text());assert cp['status']=='timeout-no-theorem' and cp['time_budget_seconds']==40 and not cp['unit_marker'];dictionary=generic_dictionary();shear=source_shear_control();required=[verify_chart('normal',t) for t in ['I','II']]+[verify_chart('parabola',t) for t in ['I','II']]+[verify_opposite()];exploratory=[verify_chart('cstar',t) for t in ['I','II']];cp1=json.loads((HERE/'cstar-general-I/checkpoint.json').read_text());assert cp1['status']=='terminal' and cp1['unit_marker'] and cp1['certificate_complete'] and cp1['log_sha256']==sha(HERE/'cstar-general-I/run.log');cp2=json.loads((HERE/'cstar-general-II/checkpoint.json').read_text());assert cp2['status']=='timeout-no-theorem' and not cp2['unit_marker'] and cp2['log_sha256']==sha(HERE/'cstar-general-II/run.log');result=dict(status='pass',source_sha256=sha(Path(__file__)),proof_sha256=sha(HERE/'PROOF.md'),premises=premises,original_input_sha256=sha(HERE/'input.json'),original_timeout_preserved=True,general_cstar_II_timeout_preserved=True,general_cstar_I_primary_unit_not_replayed_or_used=True,original_equations_reconstructed=True,original_parameters=24,original_equations=34,generic_dictionary=dictionary,source_shear_control=shear,required_certificates=required,exploratory_fixed_cstar_certificates=exploratory,all_required_certificates_multiplied_independently=True,math_coverage='Written quotient/leading-degree/curve-residue and normal-chart proof plus GLDC5/7/8 and Moh. Not formalized. Unchanged timeouts are not UNSAT.',scope='No3D Keller counterexample of weights(1,-1,-1), ordinarydegree<=7 and affine-invariant profile strictly below764. Whole nonlinear first quotient degree<=2 is excluded; for quotient degree3 the strict-profile branch is excluded by the opposite-cubic certificate. Potential degree7 profile775/773 outside this strict range is not completely closed: generalC*chartII timed out and chartI primary unit is unused.',global_witness_bounds_improved=False,global_dense_profile_count_changed=False,seconds=round(time.monotonic()-start,3));(HERE/'validation.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('VERIFIED WEIGHT111 STRICT-PROFILE EXCLUSION; GENERAL OUTSIDE-RANGE CHART REMAINS UNRESOLVED',flush=True)
if __name__=='__main__':main()
