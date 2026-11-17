"""Exact controls for primitive SAFE bifurcations, not the universal proof."""
from pathlib import Path
import hashlib,json
import sympy as S

ROOT=Path(__file__).resolve().parent
F,G,s=S.symbols('F G s')
cases=[]
for group,K,expected,points,generic_genus,infinity,central_genus in [
 ('2O',F**3*G/36-3*G**3,-s**9*(s-16)/6912,[(S.Integer(48),S.Integer(-32))],1,2,0),
 ('2I',12*F**5+G**3/144,s**10*(50000*s**2-1)/331776,[(u,-3*u/5) for u in [S.sqrt(5)/200,-S.sqrt(5)/200]],2,1,1)]:
 p=S.expand(K.subs(G,s-F))
 disc=S.factor(S.discriminant(p,F))
 assert S.expand(disc-expected)==0
 assert S.degree(p,F)==2*generic_genus+2-(infinity==1)
 assert S.diff(S.Poly(p,F).LC(),s)==0
 central=S.cancel(p.subs(s,0)/F**2)
 assert S.degree(central,F)==2*central_genus+2-(infinity==1)
 assert S.discriminant(central,F)!=0
 critical=[]
 for u,v in points:
  sub={F:u,G:v}
  assert S.simplify(K.subs(sub))==0
  assert S.simplify((S.diff(K,F)-S.diff(K,G)).subs(sub))==0
  kg=S.simplify(S.diff(K,G).subs(sub))
  numerator=S.simplify((S.diff(K,F,2)-2*S.diff(K,F,G)+S.diff(K,G,2)).subs(sub))
  assert kg!=0 and numerator!=0
  critical.append({'F':str(u),'G':str(v),'value':str(S.simplify(u+v)),
                   'K_G':str(kg),'Hessian_determinant':str(S.simplify(-2*numerator/kg**2))})
 chi=2-2*generic_genus-infinity
 chi0=2-2*central_genus-infinity-1
 chinode=2-2*(generic_genus-1)-infinity-1
 assert chi0==chinode==chi+1
 assert (1-(1+len(points)))*chi+chi0+len(points)*chinode==0
 cases.append({'group':group,'discriminant':str(disc),'zero_normalization':str(S.factor(central)),
               'critical_points':critical,'generic_genus':generic_genus,'infinity_points':infinity,
               'generic_Euler':chi,'zero_Euler':chi0,'nonzero_critical_Euler':chinode})
out={'date':'2026-11-16','completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'sympy_version':S.__version__,'cases':cases,'exact_discriminant_and_Hessian_controls':'PASS',
     'universal_branch_support_theorem_machine_verified':False,'primitivity_proof_machine_verified':False,
     'pure_SAFE_all_order_proof':'Written weight-selection argument,not finite moment tests',
     'new_holonomic_CAS_or_degree_campaign':False,'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
     'independent_replay':False,'formal_verification':False,'external_review':False}
(ROOT/'haar_branch_support_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
