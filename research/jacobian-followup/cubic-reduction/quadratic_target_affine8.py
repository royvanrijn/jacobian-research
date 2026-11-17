#!/usr/bin/env python3
"""Complete translated homogeneous matrices for H+(0,L,M)*H1^2;no solver."""
from pathlib import Path
from math import comb
import json,hashlib,time
import sympy as s
B=Path(__file__).resolve().parent;O=B/'quadratic-target-affine8'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(poly,vs):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(poly),*vs).terms() if c]
def main():
 start=time.monotonic();O.mkdir(exist_ok=True);src=B/'raw13-specialization/input.json';raw=json.loads(src.read_text());X=s.symbols('x y z');x,y,z=X;alpha,beta,gamma,L,M=s.symbols('alpha beta gamma L M');vs=(alpha,beta,gamma,L,M);centre=(alpha,beta,gamma);H=list(X)
 for atom in raw['atoms']:H[atom['row']]+=s.Rational(atom['coefficient'])*s.prod(a**e for a,e in zip(X,atom['exponent']))
 HH=[s.expand(f.subs(dict(zip(X,[x+alpha,y+beta,z+gamma])),simultaneous=True)) for f in H];p0=s.expand(HH[0].subs(dict.fromkeys(X,0)));P=s.expand(HH[0]-p0);Q=s.expand(P**2);G=[HH[0],s.expand(HH[1]+L*Q),s.expand(HH[2]+M*Q)];degrees={}
 hp=[s.Poly(f,*X) for f in HH];qp=s.Poly(Q,*X);gp=[s.Poly(f,*X) for f in G]
 for degree in range(2,9):
  mons=[(i,j,degree-i-j) for i in range(degree,-1,-1) for j in range(degree-i,-1,-1)]
  degrees[str(degree)]={'monomials':[list(e) for e in mons],'base_matrix':[[enc(f.coeff_monomial(e),vs) for e in mons] for f in hp],'square_row':[enc(qp.coeff_monomial(e),vs) for e in mons],'matrix':[[enc(f.coeff_monomial(e),vs) for e in mons] for f in gp]}
 payload={'source_sha256':sha(Path(__file__)),'source_input_sha256':sha(src),'centre_and_parameters':list(map(str,vs)),'degrees':degrees,'first_constant':enc(p0,vs),'scope':'Complete affine-output-equivalent centered-square model at every source centre,with all L/M zero branches retained.'};ip=O/'input.json';ip.write_text(json.dumps(payload,indent=2)+'\n')
 # Two symbolic equality obstructions after the top square forces x,z to be coordinate axes.
 eta,zeta=s.symbols('eta zeta');BB=y+eta*x+zeta*z;top=8*BB**3+6*M*x*x*BB+2*M*gamma*x**3+6*M*alpha*x*x*z;vtop=(x,y,z,alpha,gamma,M,eta,zeta)
 cubic=[-gamma*x**3-3*x*x*BB,12*gamma*x*x*BB+24*x*BB**2,6*x*BB*z];v3=(x,y,z,gamma,eta,zeta)
 rp=O/'rigidity.json';rp.write_text(json.dumps({'top_variables':list(map(str,vtop)),'top_cubic':enc(top,vtop),'cubic_variables':list(map(str,v3)),'cubic_outputs':[enc(f,v3) for f in cubic],'scope':'Necessary monomial/equality conditions with source x,z coordinate axes;free third source row y+eta*x+zeta*z.'},indent=2)+'\n')
 cp={'status':'complete','source_sha256':sha(Path(__file__)),'source_input_sha256':sha(src),'input_sha256':sha(ip),'rigidity_sha256':sha(rp),'matrix_entries':483,'homogeneous_degrees':list(range(2,9)),'centre_parameters':3,'shear_parameters':2,'solver_queries':0,'candidate_grids':0,'seconds':time.monotonic()-start};(O/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n');print('COMPLETE483 full translated degree2..8 entries,all5 parameters,two generic equality obstructions;no solver/grid',flush=True)
if __name__=='__main__':main()
