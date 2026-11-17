#!/usr/bin/env python3
"""Exact structural controls for the remaining full111 degree-seven gate; no solver."""
from pathlib import Path
import sympy as s
import hashlib,json,time
B=Path(__file__).resolve().parent;Q=B.parents[1];O=B/'equal-negative-full7'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(poly,variables):
 return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(poly),*variables).terms() if c]
def model(name,names,forms,constants=None):
 vs=s.symbols(' '.join(names));env=dict(zip(names,vs));polys={k:s.sympify(v,locals=env) for k,v in forms.items()}
 z={'variables':names,'polynomials':{k:enc(v,vs) for k,v in polys.items()}}
 if constants:z.update(constants)
 (O/(name+'.json')).write_text(json.dumps(z,indent=2)+'\n');return z
br=lambda A,B,y,z:s.diff(A,y)*s.diff(B,z)-s.diff(A,z)*s.diff(B,y)
def jac(f,A,B,y,z):return s.expand(f*br(A,B,y,z)+A*br(f,B,y,z)+B*br(A,f,y,z))
def main():
 started=time.monotonic();O.mkdir(exist_ok=True)
 prior=Q/'jacobian-followup/equal-negative-completion';E={e['id']:e for e in json.loads((Q/'MATH_STATUS.json').read_text())['entries']}
 old={rel:sha(Q/rel) for rel in E['KCUB23']['software_lock'] if rel.endswith(('.sing','.log','input.json','checkpoint.json'))}
 frozen=O/'prior-evidence.json'
 if frozen.exists():assert json.loads(frozen.read_text())['immutable_old_KCUB23_files']==old
 else:frozen.write_text(json.dumps({'immutable_old_KCUB23_files':old,'old40second_queries_and_unreplayed_primary_not_used':True},indent=2)+'\n')
 # Complete normalized33-parameter model, before any leading-degree reduction.
 Y,Z=s.symbols('Y Z');ps=s.symbols('p0:33');f=s.Integer(1);A=Y;C=Z;native=[];i=0
 for degree in range(1,4):
  for j in range(degree,-1,-1):f+=ps[i]*Y**j*Z**(degree-j);native.append({'row':0,'exponent':[degree+1,j,degree-j],'parameter':str(ps[i])});i+=1
 for row in [1,2]:
  for degree in range(2,5):
   for j in range(degree,-1,-1):
    term=ps[i]*Y**j*Z**(degree-j)
    if row==1:A+=term
    else:C+=term
    native.append({'row':row,'exponent':[degree-1,j,degree-j],'parameter':str(ps[i])});i+=1
 assert i==33
 names=['Y','Z']+list(map(str,ps));J=jac(f,A,C,Y,Z)
 model('whole',names,{'f':str(f),'A':str(A),'B':str(C),'J':str(J)},{'native_atoms':native,'scope':'Complete normalized weight111 ordinarydegree<=7'})
 # Complete post-prime-leading spaces; tangent blocks are free and constants zero.
 low=s.symbols('l0:5');aa=s.symbols('a0:10');bb=s.symbols('b0:6');mons=[Y,Z,Y**2,Y*Z,Z**2,Y**3,Y**2*Z,Y*Z**2,Z**3,Y**4]
 ff=1+sum(c*m for c,m in zip(low,[Y,Z,Y**2,Y*Z,Z**2]))+Y**3
 AA=sum(c*m for c,m in zip(aa,mons));BB=sum(c*m for c,m in zip(bb,[Y,Z,Y**2,Y*Z,Z**2,Y**3]));names=['Y','Z']+list(map(str,low+aa+bb));model('outer',names,{'f':str(ff),'A':str(AA),'B':str(BB),'J':str(jac(ff,AA,BB,Y,Z))})
 # Graph case: complete polynomial lifts after the affine Wronskian restriction.
 W,a,u,v,c,w,beta=s.symbols('W a u v c w beta');Ag=(a+u*W)*Y+v*(W-1);Bg=c*Y+w*(W-1);Jg=s.expand(beta*(W*br(Ag,Bg,Y,W)+Bg*s.diff(Ag,Y)-Ag*s.diff(Bg,Y)))
 model('graph',['Y','W','a','u','v','c','w','beta'],{'A':str(Ag),'B':str(Bg),'J':str(Jg)},{'zero_constant_relation':{'constant':'1','W':'1/2','W2':'1/4'}})
 # Two complete punctured-line lift types. No parameter saturation is used.
 T,U,D,b,d=s.symbols('T U D b d');fh=T*U+D
 Ah=a*T+b*U+fh*(u*T+v);Bh=c*T+d*U+w*fh
 model('hyperbola1',['T','U','D','a','b','c','d','u','v','w'],{'f':str(fh),'A':str(Ah),'B':str(Bh),'J':str(jac(fh,Ah,Bh,T,U))},{'zero_constant_relation':{'constant':'1','TU':'-D/3','T2U2':'D**2/6'}})
 Ah=a*T**2+b*U**2+fh*(u*T+v);Bh=c*T**2+w*fh
 model('hyperbola2',['T','U','D','a','b','c','u','v','w'],{'f':str(fh),'A':str(Ah),'B':str(Bh),'J':str(jac(fh,Ah,Bh,T,U))},{'zero_constant_coefficient':True})
 # Exact normalization of the remaining C* quotient, including beta=0.
 alpha=s.symbols('alpha');q=Y**2-beta*Y+alpha+beta**2;D0=1-alpha*beta-beta**3
 model('hyperbola-coordinate',['Y','Z','alpha','beta'],{'f':'1+alpha*Y+beta*Z+Y*Z+Y**3','T':str(Y+beta),'U':str(Z+q),'D':str(D0),'q':str(q)})
 files={f.name:sha(f) for f in sorted(O.glob('*.json')) if f.name not in ['checkpoint.json','validation.json','integration-validation.json']}
 result={'status':'complete','source_sha256':sha(Path(__file__)),'files':files,'whole_atoms':33,'post_leading_parameters':21,'structural_determinants':3,'new_solver_queries':0,'seconds':time.monotonic()-started}
 (O/'checkpoint.json').write_text(json.dumps(result,indent=2)+'\n');print('COMPLETE full111 structural controls:33 native atoms,21 post-leading coefficients,three small exact coefficient obstructions;no solver',flush=True)
if __name__=='__main__':main()
