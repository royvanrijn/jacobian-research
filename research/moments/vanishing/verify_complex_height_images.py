#!/usr/bin/env python3
"""Exact Cauchy/geometry controls; normalization and graph proofs are written."""
from pathlib import Path
from fractions import Fraction
import hashlib,json
import sympy as s
import verify_exceptional_bezout_maps as B
import verify_haar_shear_witnesses as H
ROOT=Path(__file__).resolve().parent
z,t,x,y=s.symbols('z t x y',real=True)
def add(a,b):
 out=dict(a)
 for e,c in b.items():
  out[e]=out.get(e,Fraction(0))+c
  if not out[e]:del out[e]
 return out
def E(poly):
 out={}
 for (a,b,c,d),v in poly.items():
  if c:B.add_term(out,(a+1,b,c-1,d),c*v)
  if d:B.add_term(out,(a,b+1,c,d-1),d*v)
 return out
def star(poly):return {(d,c,b,a):v*(-1)**(b+c) for (a,b,c,d),v in poly.items()}
def dec(rows):return {tuple(e):Fraction(c) for e,c in rows}
def main():
 C=s.log((z+1)/(z-1))/2;R=-1/(z*z-1)
 assert s.simplify(s.diff(C,z)-R)==0
 assert s.residue(R,z,-1)==s.Rational(1,2) and s.residue(R,z,1)==-s.Rational(1,2)
 assert s.residue(-1/z**2,z,0)==0
 for n in range(4):
  term=s.Rational(1,2*n+1)/z**(2*n+1)
  assert s.expand(s.series(C,z,s.oo,9).removeO()).coeff(z,-2*n-1)==s.Rational(1,2*n+1)
 # A unit positive interval has exactly the required rational derivative.
 assert s.simplify((1/(z+1)-1/(z-1))/2-R)==0
 # Shifted and rotated segments cannot silently have the same derivative.
 shifted=(1/(z+1-s.I)-1/(z-1-s.I))/2
 assert s.simplify(shifted-R)!=0
 rotated=(1/(z+s.I)-1/(z-s.I))/(2*s.I)
 assert s.simplify(rotated-R)!=0
 # Rank-one polynomial composition versus a genuine two-real-variable map.
 h=x+y;f=h**2+s.I*h**3;g=h**2-s.I*h**3
 assert s.expand(s.diff(f,x)*s.diff(g,y)-s.diff(f,y)*s.diff(g,x))==0
 f2=x+s.I*y;g2=x-s.I*y
 assert s.expand(s.diff(f2,x)*s.diff(g2,y)-s.diff(f2,y)*s.diff(g2,x))==-2*s.I
 # The allowed separating polynomial path is NOT asserted moment matching.
 p=t*t-1+s.I*t
 assert p.subs(t,1)==s.I and p.subs(t,-1)==-s.I
 assert s.solve(s.Eq(y*y-1,1-y*y),y)==[-1,1]
 area=s.integrate(2*(1-y*y),(y,-1,1));assert area==s.Rational(8,3)
 # Test the differential-rank criterion on both actual canonical heights.
 maps=ROOT/'exceptional_bezout_maps.json';data=json.loads(maps.read_text());ranks={}
 for name,row in data.items():
  tau=add(H.mul(dec(row['U']),dec(row['source_g'])),H.mul(dec(row['V']),dec(row['source_f'])))
  adj=star(tau)
  bracket=add(H.mul(B.D(tau),E(adj)),{e:-c for e,c in H.mul(E(tau),B.D(adj)).items()})
  assert bracket
  bidegrees={(a+b,c+d) for a,b,c,d in bracket};assert len(bidegrees)==1
  # A nonzero bihomogeneous polynomial cannot vanish on det=1:rescaling
  # would make it vanish on every invertible matrix,an algebraically dense set.
  ranks[name]={'rank_two':True,'bracket_terms':len(bracket),'bracket_bidegree':list(next(iter(bidegrees)))}
 result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':s.__version__,
  'target_Cauchy_derivative':'-1/(z^2-1)','target_residues':{'-1':'1/2','1':'-1/2'},'pure_Cauchy_derivative_residue':'0',
  'polynomial_path_loop_area':'8/3','loop_example_is_flat_moment_height':False,
  'canonical_height_rank_controls':ranks,'maps_sha256':hashlib.sha256(maps.read_bytes()).hexdigest(),
  'imports_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['verify_exceptional_bezout_maps.py','verify_haar_shear_witnesses.py']},
  'universal_normalization_positive_graph_and_geometry_proofs_machine_verified':False,
  'all_complex_line_valued_equivariant_heights':'EXCLUDED',
  'remaining_height_architectures':['image has interior','unary polynomial of a regular real seed with a separating SYMMETRIZED image graph'],
  'remaining_nonzero_pure_phase_architectures':['image has interior','unary polynomial of a regular real seed with a separating UNSYMMETRIZED image graph'],
  'either_remaining_architecture_sufficient':False,'exceptional_MZ_decisions':'OPEN on both models',
  'independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'complex_height_image_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
