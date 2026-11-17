#!/usr/bin/env python3
"""Prepare the full positive-axis quadratic source family and formal obstruction jets."""
from pathlib import Path
import json,hashlib,time
import sympy as s
B=Path(__file__).resolve().parent;O=B/'positive-source-quadratic'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(poly,vs):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(poly),*vs).terms() if c]
def homogeneous(poly,degree,X):return s.Add(*(co*s.prod(x**e for x,e in zip(X,mon)) for mon,co in s.Poly(s.expand(poly),*X).terms() if sum(mon)==degree))
def main():
 start=time.monotonic();O.mkdir(exist_ok=True);raw=B/'raw13-specialization/input.json';X=s.symbols('a b c');a,b,c=X;nu=s.symbols('n1 n2 n3');vs=X+nu;H=list(X)
 for atom in json.loads(raw.read_text())['atoms']:H[atom['row']]+=s.Rational(atom['coefficient'])*s.prod(x**e for x,e in zip(X,atom['exponent']))
 q=nu[0]*b*b+nu[1]*b*c+nu[2]*c*c;E=[a+q,b,c];Ei=[a-q,b,c];F=[s.expand(f.subs({a:a+q},simultaneous=True)) for f in H]
 Q=s.expand(H[0]**2);S=s.expand(H[1]**3+27*H[0]*H[2]**2-9*H[1]*H[2]);T=s.expand(H[1]**2-4*H[0]*S);nonlinear=[Q,S,T]
 obstacles=[homogeneous(q*s.diff(f,a),9,X) for f in nonlinear]+[homogeneous(q*q*s.diff(H[2],a,2)/2,9,X),homogeneous(q**3*s.diff(H[1],a,3)/6,9,X)]
 payload={'source_sha256':sha(Path(__file__)),'raw_input_sha256':sha(raw),'variables':list(map(str,vs)),'source_weights':[1,-1,-2],'parameter_weights':[3,4,5],'q':enc(q,vs),'source_map':[enc(f,vs) for f in E],'source_inverse':[enc(f,vs) for f in Ei],'composed_map':[enc(f,vs) for f in F],'composed_top_degrees':[7,9,10],'composed_top_forms':[enc(homogeneous(f,d,X),vs) for f,d in zip(F,[7,9,10])],'base_nonlinear_functions':[enc(f,X) for f in nonlinear],'degree9_obstruction_heads':[enc(f,vs) for f in obstacles],'obstruction_a_exponents':[5,2,3,1,0],'obstruction_nonzero_constants':[6,-108,144,24,12],'leading_parameter_strata':[{'conditions':'n1!=0','order':3,'leading_q':'n1*b^2'},{'conditions':'n1=0,n2!=0','order':4,'leading_q':'n2*b*c'},{'conditions':'n1=n2=0,n3!=0','order':5,'leading_q':'n3*c^2'}],'scope':'All three quadratic coefficients with the zero case retained separately. Formal weighted degeneration and Taylor coefficients,not a sampled-parameter exclusion.'}
 ip=O/'input.json';ip.write_text(json.dumps(payload,indent=2)+'\n');cp={'status':'complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'raw_input_sha256':sha(raw),'parameters':3,'obstruction_heads':5,'solver_queries':0,'formal_order_bound':15,'scope':'The maximum3d needed in the written formal argument is15;no truncation theorem is inferred from finite coefficient sampling.','seconds':time.monotonic()-start};(O/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n');print('COMPLETE full3-parameter source quadratic model,weighted conjugacy and five degree9 Taylor obstruction heads;no solver',flush=True)
if __name__=='__main__':main()
